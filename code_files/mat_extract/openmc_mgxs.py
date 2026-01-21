import os
import openmc
import openmc.mgxs as mgxs

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator


# Inputs and global options

MATERIALS_XML = "materials.xml"  # original file
BOX_LENGTH = 10.0  # side length of infinite box [cm]

# Defaults for FISSIONABLE materials (eigenvalue problems)
N_PARTICLES_FISS   = 20000
N_BATCHES_FISS     = 120
N_INACTIVE_FISS    = 20

# Lighter settings for NON-FISSIONABLE materials (fixed source)
N_PARTICLES_NONFISS = 5000
N_BATCHES_NONFISS   = 40   # total batches; all are active

# XMAS-172 group structure
xmas_edges = mgxs.GROUP_STRUCTURES["XMAS-172"]
groups = mgxs.EnergyGroups(xmas_edges)


# Decide if a material is fissionable based on its nuclides

FISSIONABLE_PREFIXES = ("U", "Pu", "Np", "Th", "Am", "Cm", "Cf")

def material_is_fissionable(mat: openmc.Material) -> bool:
    for name, _, _ in mat.nuclides:
        if name.startswith(FISSIONABLE_PREFIXES):
            return True
    return False


# Plotting utilities for MGXS

def plot_total_cross_section(mgxs_filename, material_name=None, save_plot=False):
    # Plot total and absorption cross sections from an MGXS HDF5 file.
    with h5py.File(mgxs_filename, "r") as f:
        # Get group structure (energy boundaries)
        group_edges = f.attrs["group structure"]

        # Get material name from file if not provided
        if material_name is None:
            material_name = list(f.keys())[0]

        # Get temperature group (e.g., '294K')
        temp_key = list(f[material_name].keys())[0]

        # Extract cross section data
        xs_total = f[f"{material_name}/{temp_key}/total"][:]
        xs_absorption = f[f"{material_name}/{temp_key}/absorption"][:]

    # Flip group edges to have highest energies first
    group_edges = np.flip(group_edges)

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot as histogram-style step plot
    ax.step(
        group_edges[:-1],
        xs_total,
        where="post",
        linewidth=2,
        label="Total",
        color="blue",
    )
    ax.step(
        group_edges[:-1],
        xs_absorption,
        where="post",
        linewidth=2,
        label="Absorption",
        color="red",
        linestyle="--",
    )

    # Set log scales for both axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Labels and formatting
    ax.set_xlabel("Energy [eV]", fontsize=12)
    ax.set_ylabel("Cross Section [cm$^{-1}$]", fontsize=12)
    ax.set_title(f"Cross Sections for {material_name}", fontsize=14)
    ax.grid(True, which="both", alpha=0.3, linestyle=":")
    ax.legend(fontsize=10)

    plt.tight_layout()

    # Save if requested
    if save_plot:
        plot_filename = mgxs_filename.replace(".h5", "_plot.png")
        plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {plot_filename}")

    plt.show()

    return fig, ax


def plot_all_materials_in_directory(root_dir="."):
    # Plot cross sections for all MGXS files found under root_dir.
    mgxs_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.startswith("mgxs_") and filename.endswith(".h5"):
                mgxs_files.append(os.path.join(dirpath, filename))

    print(f"Found {len(mgxs_files)} MGXS files")

    for mgxs_file in mgxs_files:
        print(f"\nPlotting {mgxs_file}")
        plot_total_cross_section(mgxs_file, save_plot=True)


# Load materials

materials = openmc.Materials.from_xml(MATERIALS_XML)

# Fix C0 -> C12/C13 and renumber IDs sequentially
for new_id, mat in enumerate(materials, start=1):
    # Fix C0 in this material
    for name, frac, frac_type in list(mat.nuclides):
        if name == "C0":
            mat.remove_nuclide("C0")
            c12_frac = frac * 0.9893   # 98.93% C-12
            c13_frac = frac * 0.0107   # 1.07% C-13
            mat.add_nuclide("C12", c12_frac, frac_type)
            mat.add_nuclide("C13", c13_frac, frac_type)
    # Renumber material ID
    mat.id = new_id

# Overwrite original materials.xml with fixed, renumbered version
materials.export_to_xml(MATERIALS_XML)
print("Fixed C0 to C12/C13 and renumbered material IDs starting at 1 in materials.xml")

all_materials = materials


# Loop over materials, one separate run per material

root_dir = os.getcwd()

for mat in all_materials:

    # Always ensure a non-empty, filesystem-safe name
    base_name = mat.name.strip() if (mat.name is not None and mat.name.strip()) else f"mat_{mat.id}"
    safe_name = base_name.replace(" ", "_")

    run_dir = os.path.join(root_dir, f"material_{mat.id}_{safe_name}")
    print(f"Processing material id={mat.id}, name='{mat.name}' in '{run_dir}'")

    # Make directory and change into it
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)
    
    # Ensure IDs start fresh for each material
    openmc.reset_auto_ids()

    # Write materials.xml containing ONLY this (already fixed) material
    one_mat = openmc.Materials([mat])
    one_mat.export_to_xml("materials.xml")

    # Single reflective box filled with this material to represent infinite medium
    L = BOX_LENGTH
    x0 = openmc.XPlane(x0=0.0, boundary_type="reflective")
    x1 = openmc.XPlane(x0=L,   boundary_type="reflective")
    y0 = openmc.YPlane(y0=0.0, boundary_type="reflective")
    y1 = openmc.YPlane(y0=L,   boundary_type="reflective")
    z0 = openmc.ZPlane(z0=0.0, boundary_type="reflective")
    z1 = openmc.ZPlane(z0=L,   boundary_type="reflective")

    region = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    cell = openmc.Cell(fill=mat, region=region)
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()

    # Settings: eigenvalue (fissionable) vs fixed-source (non-fiss)
    settings = openmc.Settings()

    if material_is_fissionable(mat):
        settings.batches   = N_BATCHES_FISS
        settings.inactive  = N_INACTIVE_FISS
        settings.particles = N_PARTICLES_FISS
        settings.run_mode  = "eigenvalue"
    else:
        settings.batches   = N_BATCHES_NONFISS
        settings.inactive  = 0
        settings.particles = N_PARTICLES_NONFISS
        settings.run_mode  = "fixed source"
        settings.max_particle_events = 200000
        
    lower_left  = (0.0, 0.0, 0.0)
    upper_right = (L,   L,   L)
    space = openmc.stats.Box(lower_left, upper_right)
    settings.source = openmc.IndependentSource(space=space)
    settings.source_rejection_fraction = 0.0
    settings.export_to_xml()
    expected_sp = f"statepoint.{settings.batches}.h5"

    # MGXS library definition for this single material
    mgxs_lib = mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.correction = None
    mgxs_lib.mgxs_types = [
        "total",
        "absorption",
        "nu-fission",
        "fission",
        "chi",
        "nu-scatter matrix",
        "multiplicity matrix",
    ]
    mgxs_lib.domain_type = "material"
    mgxs_lib.domains = [mat]
    mgxs_lib.by_nuclide = False

    mgxs_lib.check_library_for_openmc_mgxs()
    mgxs_lib.build_library()
    tallies = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies, merge=True)
    tallies.export_to_xml()

    # Run OpenMC in this subdirectory and load tallies
    openmc.run(cwd=".")

    if not os.path.exists(expected_sp):
        raise RuntimeError(f"Expected statepoint '{expected_sp}' not found.")
    sp_filename = expected_sp

    sp = openmc.StatePoint(sp_filename)
    summary = openmc.Summary("summary.h5")
    sp.link_with_summary(summary)
    mgxs_lib.load_from_statepoint(sp)

    # Create and export MGXS HDF5 file for this material
    xsdata_name = base_name  # same non-empty base name used above

    xsdata = mgxs_lib.get_xsdata(
        domain=mat,
        xsdata_name=xsdata_name,
        xs_type="macro",
        apply_domain_chi=True,
    )

    mg_file = openmc.MGXSLibrary(groups)
    mg_file.add_xsdata(xsdata)

    mg_filename = f"mgxs_{safe_name}.h5"
    mg_file.export_to_hdf5(mg_filename)

    print(f"Wrote MGXS library: {os.path.join(run_dir, mg_filename)}")

    # Plot the freshly written MGXS file for this material
    try:
        plot_total_cross_section(mg_filename, material_name=xsdata_name, save_plot=True)
    except Exception as e:
        print(f"Warning: plotting failed for {mg_filename}: {e}")

    # Return to the original root directory for the next material
    os.chdir(root_dir)
