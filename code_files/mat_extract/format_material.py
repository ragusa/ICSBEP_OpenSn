import os
import sys
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

FAIL_LOG_NAME = "failed_material_cases.txt"
DEFAULT_CASE_TIMEOUT_SEC = int(os.getenv("CASE_TIMEOUT_SEC", 3600))  # default 1 hr
MAX_WORKERS_ENV = os.getenv("MAX_WORKERS")

# Optional extra args for Material.py if needed, e.g. ["--xml"]
MATERIAL_PY_ARGS = []  # keep [] if it just takes the path

def run_worker(xml_path: str, timeout_sec: int, material_code_dir: Path) -> str:
    """
    Execute Material.py for a single materials.xml.
    Returns the xml_path on success.
    """
    start = time.monotonic()

    remaining = timeout_sec - (time.monotonic() - start)
    if remaining <= 0:
        raise subprocess.TimeoutExpired(cmd="material.py", timeout=timeout_sec)

    material_py = material_code_dir / "material.py"
    cmd = [sys.executable, str(material_py), *MATERIAL_PY_ARGS, xml_path]

    # Run from the Material.py folder to satisfy any relative imports/paths
    subprocess.run(cmd, check=True, timeout=remaining, cwd=str(material_code_dir))
    return xml_path

def discover_material_xmls(spherical_cases_dir: Path):
    """
    Find all materials/materials.xml files under spherical_cases.
    """
    xml_paths = []
    for dirpath, _, filenames in os.walk(spherical_cases_dir):
        if "materials.xml" in filenames:
            p = Path(dirpath) / "materials.xml"
            # Narrow to paths inside a 'materials' directory
            if Path(dirpath).name == "materials":
                xml_paths.append(p)
    xml_paths.sort()
    return xml_paths

def main():
    # repo_root/
    repo_root = Path(__file__).resolve().parents[2]
    # code_files/material_code/ (folder containing this script and Material.py)
    material_code_dir = Path(__file__).resolve().parent
    # spherical_cases/ (top-level alongside code_files, icsbep_original, etc.)
    spherical_cases_dir = repo_root / "spherical_cases"

    material_code_dir.mkdir(parents=True, exist_ok=True)
    fail_log_path = material_code_dir / FAIL_LOG_NAME
    with open(fail_log_path, "w") as f:
        f.write("")

    xml_paths = discover_material_xmls(spherical_cases_dir)
    total = len(xml_paths)
    max_workers = int(MAX_WORKERS_ENV) if MAX_WORKERS_ENV else (os.cpu_count() or 1)
    print(f"Discovered {total} materials.xml files; running with {max_workers} workers")

    futures = {}
    completed = 0

    def log_failure(line: str):
        with open(fail_log_path, "a") as f:
            f.write(line + "\n")

    max_workers = 1
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for idx, xml in enumerate(xml_paths, start=1):
            fut = executor.submit(run_worker, str(xml), DEFAULT_CASE_TIMEOUT_SEC, material_code_dir)
            futures[fut] = (idx, xml)

        for fut in as_completed(futures):
            idx, xml = futures[fut]
            try:
                data_path = fut.result()
                completed += 1
                print(f"-----------{completed}/{total}----------- OK: {data_path}")
            except subprocess.TimeoutExpired as e:
                completed += 1
                reason = f"TIMEOUT after {DEFAULT_CASE_TIMEOUT_SEC}s in {e.cmd}"
                print(f"-----------{completed}/{total}----------- FAIL: {xml} ({reason})")
                log_failure(f"Task {idx}/{total} | FAIL | {xml.parent} | {reason}")
            except subprocess.CalledProcessError as e:
                completed += 1
                reason = f"EXIT {e.returncode}"
                print(f"-----------{completed}/{total}----------- FAIL: {xml} ({reason})")
                log_failure(f"Task {idx}/{total} | FAIL | {xml.parent} | {reason}")
            except Exception as e:
                completed += 1
                reason = f"{type(e).__name__}: {e}"
                print(f"-----------{completed}/{total}----------- ERROR: {xml} ({reason})")
                log_failure(f"Task {idx}/{total} | ERROR | {xml.parent} | {reason}")

    print(f"Failure log written to: {fail_log_path}")

if __name__ == "__main__":
    main()
