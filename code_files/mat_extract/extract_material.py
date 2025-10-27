#!/usr/bin/env python3
from pathlib import Path
import shutil

def main():
    # Assume this script lives in mat_extract; repo root is one level up via "../"
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir / "../.."

    spherical_cases_dir = repo_root / "spherical_cases"
    icsbep_original_dir = repo_root / "icsbep_original"

    if not spherical_cases_dir.is_dir():
        raise FileNotFoundError(f"Missing folder: {spherical_cases_dir}")
    if not icsbep_original_dir.is_dir():
        raise FileNotFoundError(f"Missing folder: {icsbep_original_dir}")

    # Collect immediate child folders in spherical_cases as the destination case names
    case_names = sorted([p.name for p in spherical_cases_dir.iterdir() if p.is_dir()])

    print(f"Found {len(case_names)} case folders in spherical_cases")

    copied = 0
    skipped_missing_source = []
    skipped_missing_materials = []

    for name in case_names:
        # Source layout expected:
        # icsbep_original/<name>/openmc/(case-#?)/materials.xml
        src_case_dir = icsbep_original_dir / name
        if not src_case_dir.is_dir():
            skipped_missing_source.append(name)
            continue

        openmc_dir = src_case_dir / "openmc"
        if not openmc_dir.is_dir():
            skipped_missing_source.append(name)
            continue

        # Prefer explicit case-* subfolders; if none, treat openmc itself as "case-1"
        case_dirs = sorted([d for d in openmc_dir.iterdir() if d.is_dir() and d.name.lower().startswith("case-")])
        if not case_dirs:
            case_dirs = [openmc_dir]  # will map to "case-1" on destination

        any_copied_for_name = False

        for cd in case_dirs:
            # Resolve source materials.xml path
            src_materials = cd / "materials.xml"

            # Determine destination case label
            if cd == openmc_dir:
                dest_case_label = "case-1"
            else:
                dest_case_label = cd.name

            if not src_materials.is_file():
                skipped_missing_materials.append(f"{name}/{dest_case_label}")
                continue

            # Destination: spherical_cases/<name>/<case-#>/materials/materials.xml
            dest_case_dir = spherical_cases_dir / name / dest_case_label / "materials"
            dest_case_dir.mkdir(parents=True, exist_ok=True)

            shutil.copy2(src_materials, dest_case_dir / "materials.xml")
            copied += 1
            any_copied_for_name = True
            print(f"Copied: {name}/openmc/{dest_case_label}/materials.xml -> spherical_cases/{name}/{dest_case_label}/materials/materials.xml")

        # If there were explicit case-* folders but none had materials.xml, track at top level too
        if not any_copied_for_name and case_dirs and case_dirs != [openmc_dir]:
            skipped_missing_materials.append(f"{name}/(all cases)")

    print(f"\nCompleted. Files copied: {copied}")
    if skipped_missing_source:
        print("No matching folder in icsbep_original for:", ", ".join(skipped_missing_source))
    if skipped_missing_materials:
        print("Missing materials.xml under openmc for:", ", ".join(skipped_missing_materials))

if __name__ == "__main__":
    main()
