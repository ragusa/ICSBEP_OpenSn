#!/usr/bin/env python3
import os
import sys
import time
import pickle
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
from typing import Optional, Union, Tuple, List, Dict, Any

# Third-party
import openmc
import openmc.mgxs

# -------- configuration --------
CASE_STATUS_NAME = "case_status.txt"
DEFAULT_MATERIAL_TIMEOUT_SEC = int(os.getenv("MATERIAL_TIMEOUT_SEC", 2 * 3600))  # 2 hours
DEFAULT_MAX_WORKERS = int(os.getenv("MAX_WORKERS", os.cpu_count() or 1))
DEFAULT_GROUPS = os.getenv("MG_GROUPS", "XMAS-172")
DEFAULT_LEGENDRE_ORDER = int(os.getenv("LEGENDRE_ORDER", 7))
DEFAULT_N_BATCHES = int(os.getenv("N_BATCHES", 50))
DEFAULT_N_PARTICLES = int(os.getenv("N_PARTICLES", 5000))
DEFAULT_RUN_MODE = os.getenv("RUN_MODE", "fixed source")  # "fixed source" or "eigenvalue"

# -------- helpers from your single-material driver (adapted) --------
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

def create_homogeneous_geometry(mat: openmc.Material) -> openmc.Geometry:
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
    root_universe = openmc.Universe(cells=[mat_cell])
    return openmc.Geometry(root_universe)

def create_source(bbox: openmc.BoundingBox, energy_dist=(0.0e6, 20.0e6)) -> openmc.Source:
    ll = bbox.lower_left
    ur = bbox.upper_right
    space = openmc.stats.Box(ll, ur)
    emin, emax = energy_dist
    energy = openmc.stats.Uniform(emin, emax)
    angle = openmc.stats.Isotropic()
    return openmc.Source(space=space, angle=angle, energy=energy)

def setup_simulation(source: openmc.Source,
                     n_batches: int,
                     n_particles: int,
                     run_mode: str) -> openmc.Settings:
    settings = openmc.Settings()
    settings.source = source
    settings.batches = n_batches
    settings.particles = n_particles
    settings.run_mode = run_mode
    settings.temperature["method"] = "interpolation"
    settings.output = {"tallies": False}
    return settings

def setup_mgxs_tallies(geometry: openmc.Geometry,
                       group_edges: str,
                       legendre_order: int) -> Tuple[openmc.Tallies, openmc.mgxs.Library]:
    groups = openmc.mgxs.EnergyGroups(group_edges)
    mgxs_lib = openmc.mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.scatter_format = "legendre"
    mgxs_lib.legendre_order = legendre_order
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
    mgxs_lib.domains = geometry.get_all_material_cells().values()
    mgxs_lib.build_library()
    tallies = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies, merge=True)
    return tallies, mgxs_lib

def find_latest_statepoint(cwd: Path) -> Path:
    candidates = sorted(cwd.glob("statepoint.*.h5"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No statepoint.*.h5 found under {cwd}")
    return candidates[0]

# -------- discovery --------
def discover_material_tasks(spherical_cases_dir: Path) -> List[Dict[str, Any]]:
    """
    Find all spherical_cases/**/materials/materials.pkl and pair with materials.xml.
    For each entry in the pickle dict, create a task with a selector by id (int) if possible,
    otherwise by name (str).
    """
    tasks = []
    for dirpath, _, filenames in os.walk(spherical_cases_dir):
        if "materials.pkl" not in filenames:
            continue
        mat_dir = Path(dirpath)
        pkl_path = mat_dir / "materials.pkl"
        xml_path = mat_dir / "materials.xml"
        if not xml_path.is_file():
            continue
        # spherical_cases/<bench>/<case-#>/materials/
        parts = mat_dir.parts
        try:
            idx = parts.index("spherical_cases")
            bench = parts[idx + 1] if idx + 1 < len(parts) else "unknown-bench"
            case_label = parts[idx + 2] if idx + 2 < len(parts) else "case-1"
        except ValueError:
            bench, case_label = "unknown-bench", "case-1"

        with open(pkl_path, "rb") as f:
            mat_map = pickle.load(f)

        for key in mat_map.keys():
            selector: Tuple[str, Union[int, str]]
            mat_label: str
            try:
                selector = ("id", int(key))
                mat_label = f"id-{int(key)}"
            except Exception:
                selector = ("name", str(key))
                safe = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in str(key))
                mat_label = f"name-{safe or 'unnamed'}"
            out_dir = mat_dir / "mgxs" / mat_label
            tasks.append({
                "bench": bench,
                "case": case_label,
                "xml_path": str(xml_path),
                "selector": selector,
                "out_dir": str(out_dir),
            })
    return tasks

# -------- worker --------
def run_worker(xml_path: str,
               selector: Tuple[str, Union[int, str]],
               out_dir: str,
               group_edges: str,
               legendre_order: int,
               n_batches: int,
               n_particles: int,
               run_mode: str,
               timeout_sec: int) -> str:
    """
    Build a homogeneous box model of the selected material, export XMLs into out_dir,
    run OpenMC with a timeout, then post-process the statepoint to write one MGXS HDF5.
    Returns the full path to the .h5 file.
    """
    start = time.monotonic()
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # select material
    if selector[0] == "id":
        mat = load_single_material_from_xml(xml_path, id=int(selector[1]))
    elif selector[0] == "name":
        mat = load_single_material_from_xml(xml_path, name=str(selector[1]))
    else:
        raise ValueError(f"Unknown selector kind: {selector[0]}")

    # model
    geom = create_homogeneous_geometry(mat)
    bbox = geom.bounding_box
    source = create_source(bbox, energy_dist=(0.0e6, 20.0e6))
    settings = setup_simulation(source, n_batches=n_batches, n_particles=n_particles, run_mode=run_mode)
    tallies, mgxs_lib = setup_mgxs_tallies(geom, group_edges=group_edges, legendre_order=legendre_order)

    model = openmc.Model()
    model.materials = openmc.Materials([mat])
    model.geometry = geom
    model.settings = settings
    model.tallies = tallies

    # export all XMLs into out_dir
    orig_cwd = Path.cwd()
    try:
        os.chdir(out_dir)
        model.export_to_xml()
    finally:
        os.chdir(orig_cwd)

    # run openmc as a subprocess with a hard timeout
    remaining = timeout_sec - (time.monotonic() - start)
    if remaining <= 0:
        raise subprocess.TimeoutExpired(cmd="openmc", timeout=timeout_sec)
    subprocess.run(["openmc"], cwd=out_dir, check=True, timeout=remaining)

    # post-process: load SP, fill MGXS, and export
    sp_path = find_latest_statepoint(Path(out_dir))
    with openmc.StatePoint(sp_path) as sp:
        mgxs_lib.load_from_statepoint(sp)
    cell_names = [c.name for c in mgxs_lib.domains]
    h5_name = (cell_names[0] or (mat.name or "material")).replace(" ", "_")
    out_h5 = Path(out_dir) / f"{h5_name}.h5"
    mgxs_file = mgxs_lib.create_mg_library(xs_type="macro", xsdata_names=cell_names)
    mgxs_file.export_to_hdf5(filename=str(out_h5))
    return str(out_h5)

# -------- main driver (logs both SUCCESS and FAIL to case_status.txt) --------
def main():
    # repo_root/ (assuming this file lives under code_files/mgxs_code/)
    repo_root = Path(__file__).resolve().parents[2]
    mgxs_code_dir = Path(__file__).resolve().parent
    spherical_cases_dir = repo_root / "spherical_cases"

    mgxs_code_dir.mkdir(parents=True, exist_ok=True)
    status_path = mgxs_code_dir / CASE_STATUS_NAME

    # reset status file
    with open(status_path, "w") as f:
        f.write("timestamp,status,benchmark,case,selector,output_or_reason\n")

    # discover tasks
    tasks = discover_material_tasks(spherical_cases_dir)
    total = len(tasks)
    print(f"Discovered {total} material tasks; running with {DEFAULT_MAX_WORKERS} workers")

    def log_status(status: str, bench: str, case: str, selector: str, info: str):
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"{ts},{status},{bench},{case},{selector},{info}"
        with open(status_path, "a") as f:
            f.write(line + "\n")

    completed = 0
    with ProcessPoolExecutor(max_workers=DEFAULT_MAX_WORKERS) as executor:
        futures = {}
        for idx, t in enumerate(tasks, start=1):
            fut = executor.submit(
                run_worker,
                t["xml_path"],
                t["selector"],
                t["out_dir"],
                DEFAULT_GROUPS,
                DEFAULT_LEGENDRE_ORDER,
                DEFAULT_N_BATCHES,
                DEFAULT_N_PARTICLES,
                DEFAULT_RUN_MODE,
                DEFAULT_MATERIAL_TIMEOUT_SEC,
            )
            futures[fut] = (idx, t)

        for fut in as_completed(futures):
            idx, t = futures[fut]
            sel_kind, sel_val = t["selector"]
            sel_str = f"{sel_kind}={sel_val}"
            try:
                out_h5 = fut.result()
                completed += 1
                print(f"-----------{completed}/{total}----------- OK: {t['bench']}/{t['case']} | {sel_str} -> {out_h5}")
                log_status("SUCCESS", t["bench"], t["case"], sel_str, out_h5)
            except subprocess.TimeoutExpired as e:
                completed += 1
                reason = f"TIMEOUT after {e.timeout}s in {e.cmd}"
                print(f"-----------{completed}/{total}----------- FAIL: {t['bench']}/{t['case']} | {sel_str} ({reason})")
                log_status("FAIL", t["bench"], t["case"], sel_str, reason)
            except TimeoutError:
                completed += 1
                reason = f"FUTURE TIMEOUT after {DEFAULT_MATERIAL_TIMEOUT_SEC}s"
                print(f"-----------{completed}/{total}----------- FAIL: {t['bench']}/{t['case']} | {sel_str} ({reason})")
                log_status("FAIL", t["bench"], t["case"], sel_str, reason)
            except subprocess.CalledProcessError as e:
                completed += 1
                reason = f"EXIT {e.returncode}"
                print(f"-----------{completed}/{total}----------- FAIL: {t['bench']}/{t['case']} | {sel_str} ({reason})")
                log_status("FAIL", t["bench"], t["case"], sel_str, reason)
            except Exception as e:
                completed += 1
                reason = f"{type(e).__name__}: {e}"
                print(f"-----------{completed}/{total}----------- ERROR: {t['bench']}/{t['case']} | {sel_str} ({reason})")
                log_status("ERROR", t["bench"], t["case"], sel_str, reason)

    print(f"Case status written to: {status_path}")

if __name__ == "__main__":
    main()
