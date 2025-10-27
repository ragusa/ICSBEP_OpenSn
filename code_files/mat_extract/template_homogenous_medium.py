import os
from pathlib import Path

import openmc
import openmc.mgxs

# ---------- helper to load exactly one material from materials.xml ----------
from typing import Optional, Union

def load_single_material_from_xml(xml_path: Union[str, Path],
                                  *,
                                  name: Optional[str] = None,
                                  id: Optional[int] = None,
                                  index: Optional[int] = None) -> openmc.Material:
    xml_path = Path(xml_path)
    if sum(x is not None for x in (name, id, index)) != 1:
        raise ValueError("Specify exactly one of: name=..., id=..., index=...")
    mats = openmc.Materials.from_xml(str(xml_path))
    if name is not None:
        matches = [m for m in mats if (m.name == name)]
    elif id is not None:
        matches = [m for m in mats if (getattr(m, "id", None) == id)]
    else:
        if index < 0 or index >= len(mats):
            raise ValueError(f"index={index} is out of range (0..{len(mats)-1})")
        matches = [mats[index]]
    if not matches:
        listing = "\n".join([f"- idx {i:2d}: id={getattr(m, 'id', None)}, name={repr(m.name)}"
                             for i, m in enumerate(mats)])
        raise ValueError(f"No material found in {xml_path}.\nAvailable materials:\n{listing}")
    return matches[0]

# ---------- geometry ----------
def create_homogeneous_geometry(mat: openmc.Material, verbose: bool = False) -> openmc.Geometry:
    side = 1.0
    boundaries = {
        "left":   openmc.XPlane(x0=-side, boundary_type="reflective"),
        "right":  openmc.XPlane(x0= side, boundary_type="reflective"),
        "bottom": openmc.YPlane(y0=-side, boundary_type="reflective"),
        "top":    openmc.YPlane(y0= side, boundary_type="reflective"),
        "front":  openmc.ZPlane(z0=-side, boundary_type="reflective"),
        "back":   openmc.ZPlane(z0= side, boundary_type="reflective"),
    }
    region = (+boundaries["left"] & -boundaries["right"]
              & +boundaries["bottom"] & -boundaries["top"]
              & +boundaries["front"] & -boundaries["back"])
    mat_cell = openmc.Cell(name=mat.name or "material", fill=mat, region=region)
    if verbose:
        print(mat_cell)
    root_universe = openmc.Universe(cells=[mat_cell])
    return openmc.Geometry(root_universe)

# ---------- settings ----------
def setup_simulation(source: openmc.Source, my_path: Path,
                     n_batches: int = 50, n_particles: int = 5000,
                     src_type: str = "fixed source") -> openmc.Settings:
    settings = openmc.Settings()
    settings.source = source
    settings.batches = n_batches
    settings.particles = n_particles
    settings.run_mode = src_type
    settings.temperature["method"] = "interpolation"
    # keep minimal; files will be created under cwd when we call model.run(cwd=my_path)
    settings.output = {"tallies": False}
    my_path.mkdir(parents=True, exist_ok=True)
    return settings

# ---------- MGXS tallies ----------
def setup_mgxs_tallies(model: openmc.Model, group_edges="XMAS-172", verbose: bool = False):
    groups = openmc.mgxs.EnergyGroups(group_edges)
    mgxs_lib = openmc.mgxs.Library(model.geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.scatter_format = "legendre"
    mgxs_lib.legendre_order = 7
    mgxs_lib.mgxs_types = [
        "total",
        "absorption",
        "reduced absorption",
        "consistent nu-scatter matrix",
        "multiplicity matrix",
    ]
    mgxs_lib.tally_trigger = openmc.Trigger("std_dev", 1e-4)
    mgxs_lib.by_nuclide = False
    mgxs_lib.domain_type = "cell"
    mgxs_lib.domains = model.geometry.get_all_material_cells().values()
    mgxs_lib.build_library()
    tallies = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies, merge=True)
    model.tallies = tallies
    if verbose:
        print("MGXS domains:", [c.name for c in mgxs_lib.domains])
    return model, mgxs_lib

# ---------- post-processing ----------
def process_results(sp_file: Path, mgxs_lib: openmc.mgxs.Library,
                    my_path: Path, h5_filename: str, verbose: bool = False) -> Path:
    if sp_file is None:
        raise RuntimeError("No statepoint path returned from model.run().")
    cell_names = [cell.name for cell in mgxs_lib.domains]
    if verbose:
        print("Domains:", cell_names)
    with openmc.StatePoint(sp_file) as sp:
        mgxs_lib.load_from_statepoint(sp)
    mgxs_file = mgxs_lib.create_mg_library(xs_type="macro", xsdata_names=cell_names)
    out = my_path / f"{h5_filename}.h5"
    mgxs_file.export_to_hdf5(filename=str(out))
    return out

# ---------- source ----------
def create_source(bbox: openmc.BoundingBox, energy_dist=(0.0e6, 20.0e6)) -> openmc.Source:
    ll = bbox.lower_left
    ur = bbox.upper_right
    # only_fissionable is deprecated; use a plain Box for uniform-in-volume sampling
    space = openmc.stats.Box(ll, ur)
    emin, emax = energy_dist
    energy = openmc.stats.Uniform(emin, emax)
    angle = openmc.stats.Isotropic()
    return openmc.Source(space=space, angle=angle, energy=energy)

# ---------- driver ----------
if __name__ == "__main__":
    # choose which material to homogenize the box with (id=1 or id=2 from your XML)
    mat = load_single_material_from_xml("../../spherical_cases/heu-met-fast-028/case-1/materials/materials.xml", id=1)

    geom = create_homogeneous_geometry(mat, verbose=True)
    bbox = geom.bounding_box
    source = create_source(bbox, energy_dist=(0.0e6, 20.0e6))

    model = openmc.Model()
    model.materials = openmc.Materials([mat])
    model.geometry = geom

    my_path = Path(f"./mgxs_{mat.name or 'material'}_2025")
    settings = setup_simulation(source, my_path)
    model.settings = settings

    model, mgxs_lib = setup_mgxs_tallies(model, group_edges="XMAS-172")

    # run inside my_path and get the statepoint path
    sp_path = model.run(cwd=my_path)

    # export a single macro MGXS library named after the cell/material
    h5_out = process_results(sp_path, mgxs_lib, my_path, h5_filename=mgxs_lib.domains[0].name)
    print("MGXS library written to:", h5_out)
