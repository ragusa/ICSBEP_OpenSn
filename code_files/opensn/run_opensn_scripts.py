#!/usr/bin/env python3

from pathlib import Path
import csv
import subprocess
import re
import time

ROOT = Path("~/research/ICSBEP_OpenSn").expanduser()
OPEN_SN = Path("~/research/opensn/build/python/opensn").expanduser()
MESH_REPORT = ROOT / "code_files" / "gmsh_code" / "msh_header_report.csv"
OUTPUT_CSV = ROOT / "code_files" / "opensn" / "opensn_run_report.csv"

SCRIPT_GLOB = "spherical_cases/*/case-*/*.py"
TIMEOUT_SECONDS = 60 * 60 * 3  # 3 hours

K_PATTERNS = [
    re.compile(r"Computed k-eigenvalue:\s*([0-9Ee+\-.]+)"),
    re.compile(r"Final k-eigenvalue\s*:\s*([0-9Ee+\-.]+)")
]


def load_mesh_report(mesh_report_path):
    mesh_info = {}
    if not mesh_report_path.exists():
        return mesh_info

    with mesh_report_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rel_mesh = row["file"].strip()
            mesh_info[rel_mesh] = {
                "num_shells_csv": int(row["below_$PhysicalNames"]),
                "num_nodes": int(row["below_$Nodes"]),
            }
    return mesh_info


def parse_k_eff(text):
    for pattern in K_PATTERNS:
        m = pattern.search(text)
        if m:
            return float(m.group(1))
    return None


def read_radii(case_dir):
    radii_file = case_dir / "mesh" / "radii.txt"
    if not radii_file.exists():
        return None, None

    radii = []
    with radii_file.open("r") as f:
        for line in f:
            s = line.strip()
            if s:
                radii.append(float(s))
    return len(radii), radii


def find_mesh_for_case(case_dir):
    mesh_dir = case_dir / "mesh"
    msh_files = sorted(mesh_dir.glob("*.msh"))
    if not msh_files:
        return None
    if len(msh_files) == 1:
        return msh_files[0]
    # Prefer the n_shells_sphere file if multiple exist
    for msh in msh_files:
        if "n_shells_sphere_" in msh.name:
            return msh
    return msh_files[0]


def get_case_metadata(script_path):
    rel = script_path.relative_to(ROOT)
    parts = rel.parts
    # spherical_cases/<benchmark>/case-xx/<script.py>
    benchmark = parts[1]
    case_name = parts[2]
    return benchmark, case_name, rel.as_posix()


def run_case(script_path, mesh_info):
    benchmark, case_name, rel_script = get_case_metadata(script_path)
    case_dir = script_path.parent

    msh_path = find_mesh_for_case(case_dir)
    rel_msh = msh_path.relative_to(ROOT).as_posix() if msh_path else None

    num_shells_radii, radii = read_radii(case_dir)

    csv_shells = None
    csv_nodes = None

    if rel_msh and rel_msh in mesh_info:
        csv_shells = mesh_info[rel_msh]["num_shells_csv"]
        csv_nodes = mesh_info[rel_msh]["num_nodes"]
        if csv_nodes == 0:
            return {
                "benchmark": benchmark,
                "case_name": case_name,
                "script": rel_script,
                "mesh_file": rel_msh,
                "num_shells": csv_shells if csv_shells is not None else num_shells_radii,
                "radii": ";".join(map(str, radii)) if radii else "",
                "k_eff": "",
                "status": "skipped_bad_mesh",
                "returncode": "",
                "runtime_sec": "",
                "num_nodes": csv_nodes,
                "stdout_log": "",
                "stderr_log": "",
                "notes": "Skipped because msh_header_report.csv shows 0 nodes",
            }

    cmd = [str(OPEN_SN), "-i", script_path.name]
    start = time.time()

    try:
        result = subprocess.run(
            cmd,
            cwd=case_dir,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS
        )
        runtime_sec = round(time.time() - start, 3)

        full_output = (result.stdout or "") + "\n" + (result.stderr or "")
        k_eff = parse_k_eff(full_output)

        status = "ok" if result.returncode == 0 and k_eff is not None else "failed"

        return {
            "benchmark": benchmark,
            "case_name": case_name,
            "script": rel_script,
            "mesh_file": rel_msh or "",
            "num_shells": (
                csv_shells if csv_shells is not None
                else num_shells_radii if num_shells_radii is not None
                else ""
            ),
            "radii": ";".join(map(str, radii)) if radii else "",
            "k_eff": k_eff if k_eff is not None else "",
            "status": status,
            "returncode": result.returncode,
            "runtime_sec": runtime_sec,
            "num_nodes": csv_nodes if csv_nodes is not None else "",
            "stdout_log": (result.stdout or "").strip(),
            "stderr_log": (result.stderr or "").strip(),
            "notes": ""
        }

    except subprocess.TimeoutExpired as e:
        runtime_sec = round(time.time() - start, 3)
        stdout = e.stdout if isinstance(e.stdout, str) else (
            e.stdout.decode(errors="replace") if e.stdout else ""
        )
        stderr = e.stderr if isinstance(e.stderr, str) else (
            e.stderr.decode(errors="replace") if e.stderr else ""
        )

        return {
            "benchmark": benchmark,
            "case_name": case_name,
            "script": rel_script,
            "mesh_file": rel_msh or "",
            "num_shells": (
                csv_shells if csv_shells is not None
                else num_shells_radii if num_shells_radii is not None
                else ""
            ),
            "radii": ";".join(map(str, radii)) if radii else "",
            "k_eff": "",
            "status": "timeout",
            "returncode": "",
            "runtime_sec": runtime_sec,
            "num_nodes": csv_nodes if csv_nodes is not None else "",
            "stdout_log": stdout.strip(),
            "stderr_log": stderr.strip(),
            "notes": f"Timed out after {TIMEOUT_SECONDS} seconds"
        }

    except Exception as e:
        runtime_sec = round(time.time() - start, 3)
        return {
            "benchmark": benchmark,
            "case_name": case_name,
            "script": rel_script,
            "mesh_file": rel_msh or "",
            "num_shells": (
                csv_shells if csv_shells is not None
                else num_shells_radii if num_shells_radii is not None
                else ""
            ),
            "radii": ";".join(map(str, radii)) if radii else "",
            "k_eff": "",
            "status": "error",
            "returncode": "",
            "runtime_sec": runtime_sec,
            "num_nodes": csv_nodes if csv_nodes is not None else "",
            "stdout_log": "",
            "stderr_log": "",
            "notes": str(e)
        }


def main():
    if not OPEN_SN.exists():
        raise FileNotFoundError(f"OpenSn executable not found: {OPEN_SN}")

    mesh_info = load_mesh_report(MESH_REPORT)
    scripts = sorted(ROOT.glob(SCRIPT_GLOB))

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "benchmark",
        "case_name",
        "script",
        "mesh_file",
        "num_shells",
        "radii",
        "k_eff",
        "status",
        "returncode",
        "runtime_sec",
        "num_nodes",
        "stdout_log",
        "stderr_log",
        "notes",
    ]

    rows = []
    for script in scripts:
        row = run_case(script, mesh_info)
        rows.append(row)
        print(f"{row['benchmark']} {row['case_name']} -> {row['status']} k_eff={row['k_eff']}")

    with OUTPUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote report to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
