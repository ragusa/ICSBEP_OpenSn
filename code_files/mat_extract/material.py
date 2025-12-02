import os
import sys
import openmc
import openmc.mgxs as mgxs

# Inputs and global options

# If a path to materials.xml is provided on the command line, use it.
# Otherwise, default to "materials.xml" in the current working directory.
if len(sys.argv) > 1:
    materials_xml_path = sys.argv[1]
else:
    materials_xml_path = "materials.xml"

# Resolve to an absolute path and change into the directory containing materials.xml
materials_xml_path = os.path.abspath(materials_xml_path)
materials_dir = os.path.dirname(materials_xml_path)
os.chdir(materials_dir)

# From this point on, use just the filename, and cwd is the materials directory
MATERIALS_XML = os.path.basename(materials_xml_path)

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
    # Ensure every material has a non-empty name
    if not (mat.name and mat.name.strip()):
        mat.name = f"mat_{mat.id}"

# Overwrite original materials.xml with fixed, renumbered version
materials.export_to_xml(MATERIALS_XML)
print("Fixed C0 to C12/C13, renumbered material IDs starting at 1, and ensured non-empty material names in materials.xml")

all_materials = materials

def get_material_label(mat: openmc.Material) -> str:
    """
    Return a safe, non-empty label for this material to use in
    directory names and XSdata names.
    """
    raw = (mat.name or "").strip()
    if raw:
        return raw
    return f"mat_{mat.id}"

# Loop over materials, one separate run per material

root_dir = os.getcwd()  # this is now the materials/ directory for the given benchmark

for mat in all_materials:

    label = get_material_label(mat)
    safe_name = label.replace(" ", "_")
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
    xsdata_name = get_material_label(mat)

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

    # Return to the original root directory for the next material
    os.chdir(root_dir)
