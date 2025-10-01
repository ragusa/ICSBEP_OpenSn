import os
import sys
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import openmc
import openmc.mgxs


def create_homogeneous_geometry(mat, verbose=False):
    """Build a simple 1×1×1 cube filled with a single material and reflective BCs.

    Parameters
    ----------
    mat : openmc.Material
        The (single) material to fill the cube.
    verbose : bool
        If True, prints the created cell for quick sanity check.

    Returns
    -------
    openmc.Geometry
        Geometry with one region/cell filled by `mat` and reflective boundaries.
    """
    side = 1.0  # half-side length => box spans [-1, 1] in all directions

    # Reflective planes on all faces (infinite medium approximation)
    boundaries = {
        "left":   openmc.XPlane(x0=-side, boundary_type="reflective"),
        "right":  openmc.XPlane(x0= side, boundary_type="reflective"),
        "bottom": openmc.YPlane(y0=-side, boundary_type="reflective"),
        "top":    openmc.YPlane(y0= side, boundary_type="reflective"),
        "front":  openmc.ZPlane(z0=-side, boundary_type="reflective"),
        "back":   openmc.ZPlane(z0= side, boundary_type="reflective"),
    }

    # Region = intersection of the 6 half-spaces
    region = (
        +boundaries["left"]
        & -boundaries["right"]
        & +boundaries["bottom"]
        & -boundaries["top"]
        & +boundaries["front"]
        & -boundaries["back"]
    )

    # Single cell filled with the provided material
    mat_cell = openmc.Cell(name=mat.name, fill=mat, region=region)
    if verbose:
        print(mat_cell)

    root_universe = openmc.Universe(cells=[mat_cell])
    return openmc.Geometry(root_universe)


def setup_simulation(source, my_path, n_batches=50, n_particles=5000, src_type="fixed source"):
    """Create an OpenMC Settings object for a fixed-source run.

    Parameters
    ----------
    source : openmc.Source
        The source distribution (space/energy/angle) for the simulation.
    my_path : pathlib.Path
        Output directory; created if missing.
    n_batches : int
        Number of batches to run (fixed-source batches).
    n_particles : int
        Number of source particles per batch.
    src_type : {"fixed source", "eigenvalue"}
        Run mode string accepted by OpenMC.

    Returns
    -------
    openmc.Settings
    """
    settings = openmc.Settings()
    settings.source = source
    settings.batches = n_batches
    settings.particles = n_particles
    settings.run_mode = src_type

    # Temperature interpolation is a sensible default for tabulated data
    settings.temperature["method"] = "interpolation"

    # Keep tally HDF5 and other files under `my_path`
    settings.output = {"tallies": False, "path": str(my_path)}

    # Ensure output directory exists
    my_path.mkdir(parents=True, exist_ok=True)
    return settings


def setup_mgxs_tallies(model, group_edges="XMAS-172", verbose=False):
    """Attach MGXS tallies to the model on a per-cell basis.

    Parameters
    ----------
    model : openmc.Model
        Model whose geometry will be used to define MGXS domains.
    group_edges : str or array-like
        Named group structure (e.g., "XMAS-172") or explicit energy grid.
    verbose : bool
        If True, prints domain names (cells).

    Returns
    -------
    (openmc.Model, openmc.mgxs.Library)
    """
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(model.geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.scatter_format = "legendre"
    mgxs_lib.legendre_order = 7

    # What to compute
    mgxs_lib.mgxs_types = [
        "total",
        "absorption",
        "reduced absorption",
        "consistent nu-scatter matrix",
        "multiplicity matrix",
    ]

    # Target uncertainty (per-tally trigger)
    mgxs_lib.tally_trigger = openmc.Trigger("std_dev", 1e-4)

    # Macro XS (by_nuclide=False). Flip to True if you want per-nuclide XS.
    mgxs_lib.by_nuclide = False

    # Domains: compute MGXS over "cell" regions (one cell = one XS block)
    mgxs_lib.domain_type = "cell"
    mgxs_lib.domains = model.geometry.get_all_material_cells().values()

    # Build required tallies and attach to the model
    mgxs_lib.build_library()
    tallies = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies, merge=True)
    model.tallies = tallies

    if verbose:
        print("MGXS domains:", [c.name for c in mgxs_lib.domains])

    return model, mgxs_lib


def process_results(sp_file, mgxs_lib, my_path, h5_filename, verbose=False):
    """Load tallies from a statepoint, build an OpenMC MGXS library, and export.

    Parameters
    ----------
    sp_file : str or pathlib.Path
        Path to the statepoint file returned by model.run().
    mgxs_lib : openmc.mgxs.Library
        Configured library used to interpret tallies.
    my_path : pathlib.Path
        Output directory for exports.
    h5_filename : str
        Basename for the HDF5 MGXS library export.
    verbose : bool
        If True, prints cell names and paths.
    """
    if sp_file is None:
        raise RuntimeWarning("sp_file is None")

    cell_names = [cell.name for cell in mgxs_lib.domains]
    if verbose:
        print("all cell_names in order:\n\t", cell_names)

    # Load statepoint and fill MGXS objects with tallied data
    sp = openmc.StatePoint(sp_file)
    mgxs_lib.load_from_statepoint(sp)

    # Create an XS library (macro, by cell) and export it
    mgxs_file = mgxs_lib.create_mg_library(xs_type='macro', xsdata_names=cell_names)
    h5_file_path = str(my_path / f"{h5_filename}.h5")
    if verbose:
        print("Exporting MGXS to:", h5_file_path)
    mgxs_file.export_to_hdf5(filename=h5_file_path)

    # Optional: extract raw arrays and write CSVs (requires `export_results_to_csv`)
    # (Note: `export_results_to_csv` must exist in your environment.)
    for cell in mgxs_lib.domains:
        cell_id = cell.id
        cell_name = cell.name
        xs = {xs_type: mgxs_lib.get_mgxs(cell_id, xs_type).get_xs()
              for xs_type in mgxs_lib.mgxs_types}
        file_path = str(my_path / f"{cell_name}_")
        group_edges = mgxs_lib.energy_groups.group_edges
        # You must provide this helper; leaving call as-is per your template:
        export_results_to_csv(xs, group_edges, file_path)

    sp.close()

def create_source(bbox, energy_dist=(0.0e6, 20.0e6)):
    """Uniform-in-volume, isotropic source within the geometry bounding box.

    Parameters
    ----------
    bbox : openmc.BoundingBox
        Geometry bounding box (min/max in x/y/z).
    energy_dist : (float, float)
        (Emin, Emax) in eV for a simple flat energy distribution.

    Returns
    -------
    openmc.Source
    """
    ll = bbox.lower_left
    ur = bbox.upper_right

    # Uniform spatial distribution in the box
    space = openmc.stats.Box(ll, ur, only_fissionable=False)

    # Flat energy in [Emin, Emax]
    emin, emax = energy_dist
    energy = openmc.stats.Uniform(emin, emax)

    # Isotropic angle
    angle = openmc.stats.Isotropic()

    src = openmc.Source(space=space, angle=angle, energy=energy)
    return src

# -------------- YOU WILL INSERT THE 'mat' LOADER RIGHT HERE (see next code block) --------------

# mat = reload from XML (use the utility below), e.g.:
# from pathlib import Path
# mat = load_single_material_from_xml(Path("materials.xml"), name="Lead")  # or id=10, etc.

geom = create_homogeneous_geometry(mat, verbose=True)
print(geom.bounding_box.lower_left, geom.bounding_box.upper_right)

bbox = geom.bounding_box
energy_dist = (0.0e6, 20.0e6)

# Provide your own helper returning an openmc.Source compatible with `bbox` and `energy_dist`.
source = create_source(bbox, energy_dist)

# Assemble the OpenMC model
model = openmc.Model()
model.materials = openmc.Materials([mat])  # IMPORTANT: only the selected material
model.geometry = geom

my_path = Path(f'./mgxs_{mat.name}_2025/')
print('results path =', my_path)

settings = setup_simulation(source, my_path)
model.settings = settings

model, mgxs_lib = setup_mgxs_tallies(model)

# Small trick to avoid leaving an open StatePoint handle in notebooks
try:
    sp
    print('sp found\n')
    sp.close()
except NameError:
    print('sp NOT found\n')

file_path = my_path / f"{mat.name}.xml"
print('model file path =', file_path)

# Run and then process results
sp_file = model.run(path=file_path)
process_results(sp_file, mgxs_lib, my_path, h5_filename=mgxs_lib.domains[0].name)


