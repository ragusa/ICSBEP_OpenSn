# Material Extraction

## Purpose

This folder extracts materials for spherical cases from ICSBEP Criticaly Benchmark cases. Material files are saved in each case folder as separate subfolders for each layer. Failed cases are listed in [failed_material_cases.txt](./failed_material_cases.txt).

## How to Run

From this folder, run:
  - python3 [generate_material_mgxs.py](./generate_material_mgxs.py)

## Inputs

- icsbep_original: source tree containing {case}/openmc/(case-*|openmc)/materials.xml files expected per case.
- spherical_cases: destination tree with case folders that will receive materials/materials.xml per case-N.

## Outputs

- materials/material_{n}\_{material name}/mgxs\_{material name}.h5 files to be used in OpenSn
- [failed.txt](./failed.txt) file within [mat_extract](./) folder consisting of cases that did not complete material extraction.

## Rules and Behavior

- Repository root is inferred one level above this folder; paths are resolved relative to mat_extract.
- If openmc has explicit case-* subfolders, each is checked; otherwise openmc itself is treated as case-1.
- Destination directories are created as spherical_cases/{name}/{case-N}/materials before copying with metadata preserved.
- Cases lacking a matching icsbep_original/{name}/openmc layout are listed as missing sources; cases without materials.xml are reported per case or “(all cases)” if none had the file.

### Notes

- The script expects immediate child folders in spherical_cases to define the set of destination case names to process.
- Run from mat_extract to ensure paths resolve correctly; otherwise adjust repo_root logic if relocating the script.
