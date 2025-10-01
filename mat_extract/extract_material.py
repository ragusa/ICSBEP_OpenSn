#!/usr/bin/env python3
from pathlib import Path
import shutil

def main():
    # Assume this script lives in mat_extract, and sibling folders are spherical_cases and icsbep_original [web:16].
    script_dir = Path(__file__).resolve().parent  # mat_extract [web:16]
    repo_root = script_dir.parent  # folder containing mat_extract, spherical_cases, icsbep_original [web:16]

    spherical_cases_dir = repo_root / "spherical_cases"  # target cases root [web:16]
    icsbep_original_dir = repo_root / "icsbep_original"  # source cases root [web:16]

    if not spherical_cases_dir.is_dir():
        raise FileNotFoundError(f"Missing folder: {spherical_cases_dir}")  # sanity check [web:7]
    if not icsbep_original_dir.is_dir():
        raise FileNotFoundError(f"Missing folder: {icsbep_original_dir}")  # sanity check [web:7]

    # 1) Collect folder names directly under spherical_cases [web:12][web:18].
    case_names = sorted([p.name for p in spherical_cases_dir.iterdir() if p.is_dir()])  # flat list of dirs [web:12]
    print(f"Found {len(case_names)} case folders in spherical_cases")  # progress [web:12]

    copied = 0
    skipped_missing_source = []
    skipped_missing_materials = []

    for name in case_names:
        # 2) In icsbep_original, find matching folder, then its 'openmc/materials.xml' [web:16].
        src_case_dir = icsbep_original_dir / name  # match by folder name [web:16]
        if not src_case_dir.is_dir():
            skipped_missing_source.append(name)  # no matching source case [web:7]
            continue

        openmc_dir = src_case_dir / "openmc"  # expected structure [web:16]
        materials_xml = openmc_dir / "materials.xml"  # file to copy [web:16]

        if not materials_xml.is_file():
            skipped_missing_materials.append(name)  # no materials.xml to copy [web:7]
            continue

        # 3) In spherical_cases/<name>, create 'materials' and copy materials.xml there [web:2][web:13].
        dest_case_dir = spherical_cases_dir / name  # matched destination [web:16]
        dest_materials_dir = dest_case_dir / "materials"  # new folder [web:16]
        dest_materials_dir.mkdir(parents=True, exist_ok=True)  # ensure exists [web:13][web:10]

        # Use shutil.copy2 to preserve metadata where possible [web:2][web:5].
        shutil.copy2(materials_xml, dest_materials_dir / "materials.xml")  # copy file [web:2][web:3]
        copied += 1
        print(f"Copied: {name}/openmc/materials.xml -> spherical_cases/{name}/materials/materials.xml")  # progress [web:2]

    print(f"\nCompleted. Files copied: {copied}")  # summary [web:2]
    if skipped_missing_source:
        print("No matching folder in icsbep_original for:", ", ".join(skipped_missing_source))  # info [web:7]
    if skipped_missing_materials:
        print("Missing materials.xml under openmc for:", ", ".join(skipped_missing_materials))  # info [web:7]

if __name__ == "__main__":
    main()