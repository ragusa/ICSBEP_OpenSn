import os
import sys
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

FAIL_LOG_NAME = "failed_cases.txt"
DEFAULT_CASE_TIMEOUT_SEC = int(os.getenv("CASE_TIMEOUT_SEC", 3*3600))  # 1 hr

def run_worker(radii_list, radii_file_path, timeout_sec, gmsh_code_dir: Path) -> str:
    """
    radii_list: list[str] parsed from radii.txt
    radii_file_path: full path to the radii.txt used for logging
    gmsh_code_dir: folder containing nShell.py
    """
    start = time.monotonic()

    # Run nShell.py from gmsh_code so relative imports/paths work
    remaining = timeout_sec - (time.monotonic() - start)
    if remaining <= 0:
        raise subprocess.TimeoutExpired(cmd="nShell.py", timeout=timeout_sec)

    nshell = gmsh_code_dir / "nShell.py"
    # Ensure we execute from gmsh_code to match script expectations
    subprocess.run(
        [sys.executable, str(nshell), ",".join(radii_list), str(radii_file_path)],
        check=True,
        timeout=remaining,
        cwd=str(gmsh_code_dir),
    )

    return str(radii_file_path)

def discover_cases(spherical_cases_dir: Path):
    """
    Walk spherical_cases/**/mesh/radii.txt and collect:
    - filepaths: list[Path] to each radii.txt
    - all_radii: list[list[str]] parsed values per file
    """
    filepaths, all_radii = [], []
    for dirpath, _, filenames in os.walk(spherical_cases_dir):
        for filename in filenames:
            if filename == "radii.txt":
                filepath = Path(dirpath) / filename
                filepaths.append(filepath)
                with open(filepath, "r") as f:
                    radii = [line.strip() for line in f if line.strip()]
                all_radii.append(radii)
    return filepaths, all_radii

def main():
    # repo_root/
    repo_root = Path(__file__).resolve().parents[2]
    # code_files/gmsh_code/
    gmsh_code_dir = Path(__file__).resolve().parent
    # spherical_cases/ (top-level alongside code_files, icsbep_original, mesh, etc.)
    spherical_cases_dir = repo_root / "spherical_cases"

    gmsh_code_dir.mkdir(parents=True, exist_ok=True)
    fail_log_path = gmsh_code_dir / FAIL_LOG_NAME

    # Truncate/create the failure log at the start of each run
    with open(fail_log_path, "w") as f:
        f.write("")

    # Find all radii.txt under spherical_cases/**/mesh/
    filepaths, all_radii = discover_cases(spherical_cases_dir)

    # Optional: narrow to a specific case (keep or remove as desired)
    #filepaths = filepaths[2:3]
    #all_radii = all_radii[2:3]

    # Concurrency
    max_workers_env = os.getenv("MAX_WORKERS")
    max_workers = int(max_workers_env) if max_workers_env else (os.cpu_count() or 1)
    print(f"Discovered {len(filepaths)} cases; running with {max_workers} workers")

    futures = {}
    completed = 0
    total = len(filepaths)

    def log_failure(line: str):
        with open(fail_log_path, "a") as f:
            f.write(line + "\n")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for idx, (fp, radii_list) in enumerate(zip(filepaths, all_radii), start=1):
            fut = executor.submit(
                run_worker, radii_list, str(fp), DEFAULT_CASE_TIMEOUT_SEC, gmsh_code_dir
            )
            futures[fut] = (idx, fp, radii_list)

        for fut in as_completed(futures):
            idx, fp, radii_list = futures[fut]
            try:
                radii_str = ",".join(radii_list)
            except TypeError:
                radii_str = ",".join(str(x) for x in radii_list)

            try:
                data_path = fut.result()
                completed += 1
                print(f"-----------{completed}/{total}----------- OK: {data_path}")
            except subprocess.TimeoutExpired as e:
                completed += 1
                case_dir = Path(fp).parent
                reason = f"TIMEOUT after {DEFAULT_CASE_TIMEOUT_SEC}s in {e.cmd}"
                print(f"-----------{completed}/{total}----------- FAIL: {fp} ({reason})")
                log_failure(f"Task {idx}/{total} | FAIL | {case_dir} | radii=[{radii_str}] | {reason}")
            except subprocess.CalledProcessError as e:
                completed += 1
                case_dir = Path(fp).parent
                reason = f"EXIT {e.returncode}"
                print(f"-----------{completed}/{total}----------- FAIL: {fp} ({reason})")
                log_failure(f"Task {idx}/{total} | FAIL | {case_dir} | radii=[{radii_str}] | {reason}")
            except Exception as e:
                completed += 1
                case_dir = Path(fp).parent
                reason = f"{type(e).__name__}: {e}"
                print(f"-----------{completed}/{total}----------- ERROR: {fp} ({reason})")
                log_failure(f"Task {idx}/{total} | ERROR | {case_dir} | radii=[{radii_str}] | {reason}")

    print(f"Failure log written to: {fail_log_path}")

if __name__ == "__main__":
    main()
