from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


# -------------------- User-editable range (1-based, inclusive) --------------------
START_CASE_NUM = 1
END_CASE_NUM = None  # None means "last found"


# -------------------- Data model --------------------
@dataclass(frozen=True)
class MatCase:
    num: int                 # 1-based index in discovery list
    case_name: str           # e.g., "heu-met-fast-002"
    case_id: str             # e.g., "case-1"
    materials_dir: Path      # folder that contains materials.xml
    xml_path: Path           # full path to materials.xml


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


def run_code_in_dir(code_path: Path, workdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(code_path)],
        cwd=str(workdir),
        capture_output=True,
        text=True,
        check=False,
    )


# -------------------- Main --------------------
def main() -> int:
    script_dir = Path(__file__).resolve().parent
    code_path = script_dir / "openmc_mgxs.py"
    if not code_path.is_file():
        raise FileNotFoundError(f"Expected openmc_mgxs.py next to this script: {code_path}")

    spherical_root = find_spherical_cases_root(script_dir)
    cases = discover_materials_xml(spherical_root)
    total = len(cases)
    start, end = clamp_range(START_CASE_NUM, END_CASE_NUM, total)

    print(f"Found {total} materials.xml files under: {spherical_root}")
    if total == 0:
        return 0

    print(f"Running cases in range: {start} through {end}")

    failed_path = script_dir / "failed.txt"
    # Overwrite each run so it reflects the latest attempt
    failed_path.write_text("", encoding="utf-8")

    failures = 0
    for c in cases:
        if not (start <= c.num <= end):
            continue

        display_name = f"{c.case_name}/{c.case_id}"
        print(f"[{c.num}/{total}] Running {display_name} in {c.materials_dir}")

        try:
            proc = run_code_in_dir(code_path, c.materials_dir)
            if proc.returncode != 0:
                failures += 1
                with failed_path.open("a", encoding="utf-8") as ff:
                    ff.write(f"CASE #{c.num}: {display_name}\n")
                    ff.write(f"  xml_path: {c.xml_path}\n")
                    ff.write(f"  workdir:  {c.materials_dir}\n")
                    ff.write(f"  returncode: {proc.returncode}\n")
                    # Keep stderr; include a little stdout if helpful
                    if proc.stderr:
                        ff.write("  stderr:\n")
                        ff.write(proc.stderr.strip() + "\n")
                    if proc.stdout:
                        ff.write("  stdout (tail):\n")
                        ff.write("\n".join(proc.stdout.splitlines()[-40:]) + "\n")
                    ff.write("\n")
        except Exception as e:
            failures += 1
            with failed_path.open("a", encoding="utf-8") as ff:
                ff.write(f"CASE #{c.num}: {display_name}\n")
                ff.write(f"  xml_path: {c.xml_path}\n")
                ff.write(f"  workdir:  {c.materials_dir}\n")
                ff.write(f"  exception: {type(e).__name__}: {e}\n\n")

    if failures:
        print(f"Done with {failures} failures. See: {failed_path}")
        return 2

    print("Done with no failures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
