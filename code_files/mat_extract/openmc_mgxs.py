import os
import openmc
import openmc.mgxs as mgxs

import h5py
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path


# ── Inputs and global options ─────────────────────────────────────────────────
SHOW_GRAPH = False
MATERIALS_XML = "materials.xml"  # original file
BOX_LENGTH = 10.0  # side length of infinite box [cm]

# Defaults for FISSIONABLE materials (eigenvalue problems)
N_PARTICLES_FISS = 300000
N_BATCHES_FISS = 520
N_INACTIVE_FISS = 120

# Lighter settings for NON-FISSIONABLE materials (fixed source)
N_PARTICLES_NONFISS = 5000
N_BATCHES_NONFISS = 40
N_INACTIVE_NONFISS = 0

# LANL 70 group structure
group_edges = np.loadtxt("../../../../code_files/mat_extract/LANL70g_eV.txt")
groups = mgxs.EnergyGroups(group_edges)

# Heuristic fissionable detection
FISSIONABLE_PREFIXES = ("U", "Pu", "Np", "Th", "Am", "Cm", "Cf")

MATERIALS_XML_IN  = "materials.xml"         # tracked/original
MATERIALS_XML_OUT = "materials_fixed.xml"   # generated

# ── Helper utilities ──────────────────────────────────────────────────────────
def sanitize_name(name: str, replacement: str = "_") -> str:
    # Make a string safe for filesystem paths and HDF5 group keys.
    if name is None:
        name = ""
    s = str(name).strip()
    s = s.replace("/", replacement).replace("\\", replacement)
    s = s.replace(" ", replacement)
    while replacement * 2 in s:
        s = s.replace(replacement * 2, replacement)
    s = s.strip(replacement)
    return s if s else "material"


def material_is_fissionable(mat: openmc.Material) -> bool:
    for name, _, _ in mat.nuclides:
        if name.startswith(FISSIONABLE_PREFIXES):
            return True
    return False


def export_results_to_csv(xs, group_edges, file_path, verbose=False):
    # Export cross-section data to CSV files.
    for xs_type, data in xs.items():
        if verbose:
            print(xs_type)
        filename = file_path + f"xs_{xs_type.replace(' ', '_')}.csv"
        if data.ndim == 3:
            pn = data.shape[-1]
            aux = data[:, :, 0]
            for ip in range(1, pn):
                aux = np.vstack((aux, data[:, :, ip]))
            np.savetxt(filename, aux, delimiter=",")
        else:
            np.savetxt(filename, data, delimiter=",")
    ng = len(group_edges) - 1
    np.savetxt(file_path + f"energy_edges_{ng}g.csv", group_edges, delimiter=",")
    print("Cross-section data saved to CSV files.")


def process_results(sp_filename, mgxs_lib, my_path, h5_filename, verbose=False):
    # Process OpenMC statepoint results and export MGXS HDF5 + CSV data.
    if sp_filename is None:
        raise RuntimeError("sp_filename is None")

    cell_names = [cell.name for cell in mgxs_lib.domains]
    if verbose:
        print("all cell_names in order:\n\t", cell_names)

    sp = openmc.StatePoint(sp_filename)
    summary = openmc.Summary(os.path.join(os.path.dirname(sp_filename), "summary.h5"))
    sp.link_with_summary(summary)
    mgxs_lib.load_from_statepoint(sp)

    mgxs_file = mgxs_lib.create_mg_library(xs_type="macro", xsdata_names=cell_names)

    h5_file_path = str(my_path / f"{h5_filename}.h5")
    print(h5_file_path)
    mgxs_file.export_to_hdf5(filename=h5_file_path)

    for cell in mgxs_lib.domains:
        cell_id = cell.id
        cell_name = cell.name
        xs = {
            xs_type: mgxs_lib.get_mgxs(cell_id, xs_type).get_xs()
            for xs_type in mgxs_lib.mgxs_types
        }
        file_path_csv = str(my_path / f"{cell_name}_")
        group_edges_arr = mgxs_lib.energy_groups.group_edges
        export_results_to_csv(xs, group_edges_arr, file_path_csv)

    sp.close()


# ── Plotting utilities ────────────────────────────────────────────────────────
def plot_total_cross_section(
    mgxs_filename,
    material_name=None,
    material_key=None,
    display_name=None,
    save_plot=False,
):
    # Plot total and absorption cross sections from an MGXS HDF5 file.
    with h5py.File(mgxs_filename, "r") as f:
        group_edges_local = f.attrs["group structure"]

        if material_key is None:
            material_key = material_name
        if material_key is None:
            material_key = list(f.keys())[0]
        if display_name is None:
            display_name = material_key

        temp_key = list(f[material_key].keys())[0]
        xs_total = f[f"{material_key}/{temp_key}/total"][:]
        xs_absorption = f[f"{material_key}/{temp_key}/absorption"][:]

    group_edges_local = np.flip(group_edges_local)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(group_edges_local[:-1], xs_total, where="post", linewidth=2,
            label="Total", color="blue")
    ax.step(group_edges_local[:-1], xs_absorption, where="post", linewidth=2,
            label="Absorption", color="red", linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("linear")
    ax.set_xlabel("Energy [eV]", fontsize=12)
    ax.set_ylabel("Cross Section [cm$^{-1}$]", fontsize=12)
    ax.set_title(f"Cross Sections for {display_name}", fontsize=14)
    ax.grid(True, which="both", alpha=0.3, linestyle=":")
    ax.legend(fontsize=10)
    plt.tight_layout()

    if save_plot:
        plot_filename = mgxs_filename.replace(".h5", "_plot.png")
        plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {plot_filename}")
    if SHOW_GRAPH:
        plt.show()
    plt.close(fig)

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


# ── Load and fix materials ────────────────────────────────────────────────────
materials = openmc.Materials.from_xml(MATERIALS_XML_IN)

your_files = os.getcwd()
if "ragusa" in your_files:
    os.environ["OPENMC_CROSS_SECTIONS"] = (
        "/home/ragusa/xs/endfb-viii.0-hdf5/cross_sections.xml"
    )

# Fix C0 -> C12/C13 and renumber IDs sequentially
for new_id, mat in enumerate(materials, start=1):
    for name, frac, frac_type in list(mat.nuclides):
        if name == "C0":
            mat.remove_nuclide("C0")
            c12_frac = frac * 0.9893   # 98.93% C-12
            c13_frac = frac * 0.0107   # 1.07% C-13
            mat.add_nuclide("C12", c12_frac, frac_type)
            mat.add_nuclide("C13", c13_frac, frac_type)
    mat.id = new_id
    
materials.export_to_xml(MATERIALS_XML_OUT)
all_materials = openmc.Materials.from_xml(MATERIALS_XML_OUT)
print("Fixed C0 to C12/C13 and renumbered material IDs starting at 1.")

all_materials = materials

# ── Main loop: one separate run per material ──────────────────────────────────
root_dir = os.getcwd()

for mat in all_materials:
    # ── Name handling ─────────────────────────────────────────────────────
    base_name = (
        mat.name.strip()
        if (mat.name is not None and mat.name.strip())
        else f"mat_{mat.id}"
    )
    safe_name = sanitize_name(base_name)

    run_dir = os.path.join(root_dir, f"material_{mat.id}_{safe_name}")
    print(f"\nProcessing material id={mat.id}, name='{mat.name}' in '{run_dir}'")

    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)

    # Ensure IDs start fresh for each material
    openmc.reset_auto_ids()

    # ── Materials: write then reload fresh to avoid registry collisions ───
    openmc.Materials([mat]).export_to_xml("materials.xml")
    one_mat = openmc.Materials.from_xml("materials.xml")[0]

    # ── Geometry (reflective box → infinite homogeneous medium) ───────────
    L = BOX_LENGTH
    x0 = openmc.XPlane(x0=0.0, boundary_type="reflective")
    x1 = openmc.XPlane(x0=L, boundary_type="reflective")
    y0 = openmc.YPlane(y0=0.0, boundary_type="reflective")
    y1 = openmc.YPlane(y0=L, boundary_type="reflective")
    z0 = openmc.ZPlane(z0=0.0, boundary_type="reflective")
    z1 = openmc.ZPlane(z0=L, boundary_type="reflective")

    region = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    cell = openmc.Cell(name=safe_name, fill=one_mat, region=region)
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()

    bbox = geometry.bounding_box

    # ── Source: only apply fissionable constraint when it makes sense ─────
    uniform_dist = openmc.stats.Box(bbox.lower_left, bbox.upper_right)

    if material_is_fissionable(one_mat):
        source = openmc.IndependentSource(
            space=uniform_dist,
            constraints={"fissionable": True},
        )
        print("  Using source constraint: fissionable=True")
    else:
        source = openmc.IndependentSource(space=uniform_dist)
        print("  Using unconstrained uniform box source (nonfissionable material)")

    # ── Settings ──────────────────────────────────────────────────────────
    settings = openmc.Settings()
    settings.source = source    

    if material_is_fissionable(mat):
        settings.batches   = N_BATCHES_FISS
        settings.inactive  = N_INACTIVE_FISS
        settings.particles = N_PARTICLES_FISS
        settings.run_mode  = "eigenvalue"
    else:
        settings.batches   = N_BATCHES_NONFISS
        settings.inactive  = N_INACTIVE_NONFISS
        settings.particles = N_PARTICLES_NONFISS
        settings.run_mode  = "fixed source"
        settings.max_particle_events = 200000
    
    settings.temperature["method"] = "interpolation"
    settings.export_to_xml()
    expected_sp = f"statepoint.{settings.batches}.h5"


    # ── MGXS library ─────────────────────────────────────────────────────
    mgxs_lib = mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups

    mgxs_lib.scatter_format = "legendre"
    mgxs_lib.legendre_order = 7

    mgxs_lib.mgxs_types = [
        "total",
        "absorption",
        "fission",
        "nu-fission",
        "chi",
        "reduced absorption",
        "scatter matrix",
        "nu-scatter matrix",
        "consistent nu-scatter matrix",
        "multiplicity matrix",
    ]

    mgxs_lib.by_nuclide = False

    mgxs_lib.domain_type = "cell"
    mgxs_lib.domains = list(geometry.get_all_material_cells().values())

    mgxs_lib.build_library()
    tallies = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies, merge=True)
    tallies.export_to_xml()

    # ── Run OpenMC ───────────────────────────────────────────────────────
    openmc.run(cwd=".")

    if not os.path.exists(expected_sp):
        raise RuntimeError(f"Expected statepoint '{expected_sp}' not found.")

    # ── Process results & export MGXS HDF5 + CSV ─────────────────────────
    my_path = Path(".")
    ng = len(groups.group_edges) - 1
    gr_name = f"_LANL{ng}g"
    h5_filename = f"{safe_name}{gr_name}"

    process_results(expected_sp, mgxs_lib, my_path, h5_filename=h5_filename)

    # ── Plot the freshly written MGXS file ───────────────────────────────
    mgxs_h5_path = str(my_path / f"{h5_filename}.h5")
    try:
        plot_total_cross_section(
            mgxs_h5_path,
            material_key=safe_name,
            display_name=base_name,
            save_plot=True,
        )
    except Exception as e:
        print(f"Warning: plotting failed for {mgxs_h5_path}: {e}")

    # Return to the original root directory for the next material
    os.chdir(root_dir)

print("\n=== All materials processed. ===")
