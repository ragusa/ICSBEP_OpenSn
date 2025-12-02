# Material Extraction

## Purpose

This folder extracts materials for spherical cases from ICSBEP Criticaly Benchmark cases. Material files are saved in each case folder as separate subfolders for each layer. Failed cases are listed in [failed_material_cases.txt](./failed_material_cases.txt).

## How to Run

From this folder, run:
  - python3 [extract_material.py](./extract_material.py)
  - python3 [format_material.py](./format_material.py)

## Inputs

- icsbep_original: source tree containing {case}/openmc/(case-*|openmc)/materials.xml files expected per case.
- spherical_cases: destination tree with case folders that will receive materials/materials.xml per case-N.

## Outputs

- materials/materials.xml created under spherical_cases/{name}/{case-N}/ for each discovered case with a source materials.xml.
- materials/material_{n}\_{material name}/mgxs\_{material name}.h5 files to be used in OpenSn
- [failed_material_cases.txt](./failed_material_cases.txt) file within [mat_extract](./) folder consisting of cases that did not complete material extraction.

## Rules and Behavior

- Repository root is inferred one level above this folder; paths are resolved relative to mat_extract.
- If openmc has explicit case-* subfolders, each is checked; otherwise openmc itself is treated as case-1.
- Destination directories are created as spherical_cases/{name}/{case-N}/materials before copying with metadata preserved.
- Cases lacking a matching icsbep_original/{name}/openmc layout are listed as missing sources; cases without materials.xml are reported per case or “(all cases)” if none had the file.

### Notes

- The script expects immediate child folders in spherical_cases to define the set of destination case names to process.
- Run from mat_extract to ensure paths resolve correctly; otherwise adjust repo_root logic if relocating the script.
