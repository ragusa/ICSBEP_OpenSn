from __future__ import annotations

import os
import subprocess
import sys
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, List, Optional, Tuple


# -------------------- User-editable range (1-based, inclusive) --------------------
START_CASE_NUM = 1
END_CASE_NUM = None  # None means "last found"

# Only forward these child lines to the console (everything else is hidden)
FORWARD_PREFIXES = ("TIMING:", "ERROR:", "WARN:")


# -------------------- Data model --------------------
@dataclass(frozen=True)
class MatCase:
    num: int
    case_name: str
    case_id: str
    materials_dir: Path
    xml_path: Path


# -------------------- Helpers --------------------
def find_spherical_cases_root(script_dir: Path) -> Path:
    case_loc = (script_dir / ".." / ".." / "spherical_cases").resolve()
    if case_loc.is_dir():
        return case_loc
    raise FileNotFoundError("Could not locate 'spherical_cases' directory relative to this script.")


def parse_case_name_and_id(spherical_root: Path, xml_path: Path) -> Tuple[str, str]:
    rel = xml_path.resolve().relative_to(spherical_root.resolve())
    parts = rel.parts
    if len(parts) < 4:
        return ("UNKNOWN_CASE", "UNKNOWN_ID")
    return (parts[0], parts[1])


def discover_materials_xml(spherical_root: Path) -> List[MatCase]:
    xmls = sorted(spherical_root.glob("**/materials/materials.xml"))
    cases: List[MatCase] = []
    for i, xml_path in enumerate(xmls, start=1):
        case_name, case_id = parse_case_name_and_id(spherical_root, xml_path)
        cases.append(
            MatCase(
                num=i,
                case_name=case_name,
                case_id=case_id,
                materials_dir=xml_path.parent,
                xml_path=xml_path,
            )
        )
    return cases


def clamp_range(start_1based: int, end_1based: Optional[int], total: int) -> Tuple[int, int]:
    if total == 0:
        return (1, 0)
    start = max(1, min(start_1based, total))
    end = total if end_1based is None else max(1, min(end_1based, total))
    if end < start:
        start, end = end, start
    return start, end


def format_seconds(seconds: float) -> str:
    seconds = float(seconds)
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = seconds % 60.0
    if h > 0:
        return f"{h:d}:{m:02d}:{s:05.2f}"
    return f"{m:d}:{s:05.2f}"


def run_code_in_dir_filtered(
    code_path: Path,
    workdir: Path,
    *,
    forward_prefixes: Tuple[str, ...] = FORWARD_PREFIXES,
    tail_lines: int = 200,
) -> Tuple[int, str]:
    """
    Run the child script in workdir.

    - Does NOT show OpenMC (or other) output.
    - Only forwards child lines that start with forward_prefixes.
    - Returns (returncode, combined_output_tail) for failure logging.
    """
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    cmd = [sys.executable, "-u", str(code_path)]
    proc = subprocess.Popen(
        cmd,
        cwd=str(workdir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )

    assert proc.stdout is not None
    tail: Deque[str] = deque(maxlen=tail_lines)

    for line in proc.stdout:
        tail.append(line.rstrip("\n"))
        s = line.lstrip()
        if s.startswith(forward_prefixes):
            sys.stdout.write(line)
            sys.stdout.flush()

    proc.wait()
    return proc.returncode, "\n".join(tail)


# -------------------- Main --------------------
def main() -> int:
    overall_t0 = time.perf_counter()

    script_dir = Path(__file__).resolve().parent
    code_path = script_dir / "openmc_mgxs.py"
    if not code_path.is_file():
        raise FileNotFoundError(f"Expected openmc_mgxs.py next to this script: {code_path}")

    spherical_root = find_spherical_cases_root(script_dir)
    cases = discover_materials_xml(spherical_root)
    total = len(cases)
    start, end = clamp_range(START_CASE_NUM, END_CASE_NUM, total)

    print(f"Found {total} materials.xml files under: {spherical_root}", flush=True)
    if total == 0:
        return 0

    print(f"Running cases in range: {start} through {end}", flush=True)

    failed_path = script_dir / "failed.txt"
    failed_path.write_text("", encoding="utf-8")

    failures = 0
    ran = 0

    for c in cases:
        if not (start <= c.num <= end):
            continue

        ran += 1
        display_name = f"{c.case_name}/{c.case_id}"
        print("\n" + "=" * 88, flush=True)
        print(f"[{c.num}/{total}] Running {display_name} in {c.materials_dir}", flush=True)

        case_t0 = time.perf_counter()
        try:
            returncode, tail = run_code_in_dir_filtered(code_path, c.materials_dir)
        except Exception as e:
            returncode = 999
            tail = f"Helper exception while launching child:\n{type(e).__name__}: {e}"
        case_dt = time.perf_counter() - case_t0

        print(f"[{c.num}/{total}] Case time: {format_seconds(case_dt)} (rc={returncode})", flush=True)

        if returncode != 0:
            failures += 1
            with failed_path.open("a", encoding="utf-8") as ff:
                ff.write(f"CASE #{c.num}: {display_name}\n")
                ff.write(f"  xml_path: {c.xml_path}\n")
                ff.write(f"  workdir:  {c.materials_dir}\n")
                ff.write(f"  returncode: {returncode}\n")
                ff.write(f"  case_time_seconds: {case_dt:.6f}\n")
                ff.write("  combined_output_tail:\n")
                ff.write(tail.strip() + "\n\n")

    overall_dt = time.perf_counter() - overall_t0
    print("\n" + "=" * 88, flush=True)
    print(f"Helper finished. Ran {ran} case(s) in {format_seconds(overall_dt)} total.", flush=True)

    if failures:
        print(f"Done with {failures} failures. See: {failed_path}", flush=True)
        return 2

    print("Done with no failures.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
