import os
import time
import openmc
import openmc.mgxs as mgxs

import h5py
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path


# ── Inputs and global options ─────────────────────────────────────────────────
SHOW_GRAPH = False
BOX_LENGTH = 10.0  # side length of infinite box [cm]

# Defaults for FISSIONABLE materials (eigenvalue problems)
N_PARTICLES_FISS = 200000
N_BATCHES_FISS = 360
N_INACTIVE_FISS = 60

# Lighter settings for NON-FISSIONABLE materials (fixed source)
N_PARTICLES_NONFISS = 5000
N_BATCHES_NONFISS = 40
N_INACTIVE_NONFISS = 0

# LANL 70 group structure
group_edges = np.loadtxt("../../../../code_files/mat_extract/LANL70g_eV.txt")
groups = mgxs.EnergyGroups(group_edges)

MATERIALS_XML_IN = "materials.xml"
MATERIALS_XML_OUT = "materials_fixed.xml"


def _fmt_seconds(seconds: float) -> str:
    seconds = float(seconds)
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = seconds % 60.0
    if h > 0:
        return f"{h:d}:{m:02d}:{s:05.2f}"
    return f"{m:d}:{s:05.2f}"


def sanitize_name(name: str, replacement: str = "_") -> str:
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
    """
      - Return True for HEU (e.g., U-235 enrichment >> 5%)
      - Return False for natural/depleted uranium (~0.7% U-235)
      - Return True for typical Pu-bearing fuels even at modest fractions
    """

    # Fissile isotopes that strongly indicate eigenvalue-worthy fuel
    fissile_set = {
        "U233", "U235",
        "Pu239", "Pu241",
        "Am242m"
    }

    # Actinide element prefixes (used to define "heavy metal" inventory)
    actinide_prefixes = ("U", "Pu", "Th", "Np", "Am", "Cm", "Cf")

    # Collect atomic amounts
    ao = {name: amt for (name, amt, _frac_type) in mat.nuclides}

    # Total heavy metal (actinides) present
    heavy_metal_ao = sum(v for k, v in ao.items() if k.startswith(actinide_prefixes))
    if heavy_metal_ao <= 0.0:
        return False

    fissile_ao = sum(ao.get(n, 0.0) for n in fissile_set)
    if fissile_ao <= 0.0:
        return False

    # Special case: uranium-only (or uranium-dominated) → use U-235 enrichment threshold
    u_total = sum(v for k, v in ao.items() if k.startswith("U"))
    u235 = ao.get("U235", 0.0)
    u_is_dominant = (u_total / heavy_metal_ao) >= 0.95  # mostly uranium

    if u_is_dominant and u_total > 0.0:
        enrichment = u235 / u_total  # fraction (not %)
        # Conservative cutoff: treat <= ~5% as not eigenvalue-worthy in this "bare box" context
        return enrichment >= 0.05

    # General case (e.g., Pu-bearing, mixed actinides): use fissile fraction among heavy metal
    fissile_frac = fissile_ao / heavy_metal_ao

    # Pu-bearing fuels can be eigenvalue-worthy at lower fissile fraction than LEU
    pu_present = any(k.startswith("Pu") for k in ao.keys())
    if pu_present:
        return fissile_frac >= 0.01  # 1% fissile among heavy metal

    # Otherwise, require a modest fissile share
    return fissile_frac >= 0.02  # 2% fissile among heavy metal



def export_results_to_csv(xs, group_edges, file_path):
    for xs_type, data in xs.items():
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


def process_results(sp_filename, mgxs_lib, my_path, h5_filename):
    if sp_filename is None:
        raise RuntimeError("sp_filename is None")

    cell_names = [cell.name for cell in mgxs_lib.domains]

    sp = openmc.StatePoint(sp_filename)
    summary = openmc.Summary(os.path.join(os.path.dirname(sp_filename), "summary.h5"))
    sp.link_with_summary(summary)
    mgxs_lib.load_from_statepoint(sp)

    mgxs_file = mgxs_lib.create_mg_library(xs_type="macro", xsdata_names=cell_names)
    mgxs_file.export_to_hdf5(filename=str(my_path / f"{h5_filename}.h5"))

    for cell in mgxs_lib.domains:
        cell_id = cell.id
        cell_name = cell.name
        xs = {
            xs_type: mgxs_lib.get_mgxs(cell_id, xs_type).get_xs()
            for xs_type in mgxs_lib.mgxs_types
        }
        export_results_to_csv(xs, mgxs_lib.energy_groups.group_edges, str(my_path / f"{cell_name}_"))

    sp.close()


def plot_total_cross_section(mgxs_filename, material_key, display_name, save_plot=True):
    with h5py.File(mgxs_filename, "r") as f:
        group_edges_local = np.flip(f.attrs["group structure"])
        temp_key = list(f[material_key].keys())[0]
        xs_total = f[f"{material_key}/{temp_key}/total"][:]
        xs_absorption = f[f"{material_key}/{temp_key}/absorption"][:]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(group_edges_local[:-1], xs_total, where="post", linewidth=2, label="Total", color="blue")
    ax.step(group_edges_local[:-1], xs_absorption, where="post", linewidth=2, label="Absorption", color="red", linestyle="--")
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
    if SHOW_GRAPH:
        plt.show()
    plt.close(fig)


# ── Load and fix materials ────────────────────────────────────────────────────
script_t0 = time.perf_counter()

materials = openmc.Materials.from_xml(MATERIALS_XML_IN)

your_files = os.getcwd()
if "ragusa" in your_files:
    os.environ["OPENMC_CROSS_SECTIONS"] = (
        "/home/ragusa/xs/endfb-viii.0-hdf5/cross_sections.xml"
    )

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

# Keep your original behavior (you overwrite with in-memory materials)
all_materials = materials


# ── Main loop: one separate run per material ──────────────────────────────────
root_dir = os.getcwd()
all_materials_list = list(all_materials)
n_total = len(all_materials_list)

print(f"TIMING: case_start materials={n_total}", flush=True)

for i_mat, mat in enumerate(all_materials_list, start=1):
    mat_t0 = time.perf_counter()

    base_name = mat.name.strip() if (mat.name is not None and mat.name.strip()) else f"mat_{mat.id}"
    safe_name = sanitize_name(base_name)
    run_dir = os.path.join(root_dir, f"material_{mat.id}_{safe_name}")

    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)

    openmc.reset_auto_ids()

    openmc.Materials([mat]).export_to_xml("materials.xml")
    one_mat = openmc.Materials.from_xml("materials.xml")[0]

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
    uniform_dist = openmc.stats.Box(bbox.lower_left, bbox.upper_right)

    if material_is_fissionable(one_mat):
        source = openmc.IndependentSource(space=uniform_dist, constraints={"fissionable": True})
        fiss_flag = 1
    else:
        source = openmc.IndependentSource(space=uniform_dist)
        fiss_flag = 0

    settings = openmc.Settings()
    settings.source = source

    if material_is_fissionable(mat):
        settings.batches = N_BATCHES_FISS
        settings.inactive = N_INACTIVE_FISS
        settings.particles = N_PARTICLES_FISS
        settings.run_mode = "eigenvalue"
    else:
        settings.batches = N_BATCHES_NONFISS
        settings.inactive = N_INACTIVE_NONFISS
        settings.particles = N_PARTICLES_NONFISS
        settings.run_mode = "fixed source"
        settings.max_particle_events = 200000

    settings.temperature["method"] = "interpolation"
    settings.export_to_xml()
    expected_sp = f"statepoint.{settings.batches}.h5"

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

    run_t0 = time.perf_counter()
    openmc.run(cwd=".", output=False)  # suppress OpenMC console output
    run_dt = time.perf_counter() - run_t0

    if not os.path.exists(expected_sp):
        raise RuntimeError(f"Expected statepoint '{expected_sp}' not found.")

    post_t0 = time.perf_counter()
    my_path = Path(".")
    ng = len(groups.group_edges) - 1
    h5_filename = f"{safe_name}_LANL{ng}g"

    process_results(expected_sp, mgxs_lib, my_path, h5_filename=h5_filename)

    mgxs_h5_path = str(my_path / f"{h5_filename}.h5")
    try:
        plot_total_cross_section(mgxs_h5_path, material_key=safe_name, display_name=base_name, save_plot=True)
    except Exception:
        # Keep quiet by default; helper will still see timing.
        pass

    post_dt = time.perf_counter() - post_t0

    os.chdir(root_dir)

    mat_dt = time.perf_counter() - mat_t0
    print(
        "TIMING: material_done "
        f"{i_mat}/{n_total} id={mat.id} fiss={fiss_flag} "
        f"total={_fmt_seconds(mat_dt)} openmc={_fmt_seconds(run_dt)} post={_fmt_seconds(post_dt)} "
        f"name={safe_name}",
        flush=True,
    )

script_dt = time.perf_counter() - script_t0
print(f"TIMING: case_total total={_fmt_seconds(script_dt)}", flush=True)
