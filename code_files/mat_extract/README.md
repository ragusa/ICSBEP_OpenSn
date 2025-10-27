# Material Extraction

## Extraction Process

Use the [material extraction script](./extract_material.py) to copy OpenMC materials.xml files from ICSBEP originals into the normalized spherical case tree under spherical_cases, creating a materials subfolder per case as needed.

## How to Run

- Run the script from this folder:
  python3 [extract_material.py](./extract_material.py)

## Inputs

- icsbep_original: source tree containing <case>/openmc/(case-*|openmc)/materials.xml files expected per case.
- spherical_cases: destination tree with case folders that will receive materials/materials.xml per case-N.

## Outputs

- materials/materials.xml created under spherical_cases/<name>/<case-N>/ for each discovered case with a source materials.xml.
- Console summary of copied counts, missing source case folders, and missing materials.xml paths.

## Rules and Behavior

- Repository root is inferred one level above this folder; paths are resolved relative to mat_extract.
- If openmc has explicit case-* subfolders, each is checked; otherwise openmc itself is treated as case-1.
- Destination directories are created as spherical_cases/<name>/<case-N>/materials before copying with metadata preserved (copy2).
- Cases lacking a matching icsbep_original/<name>/openmc layout are listed as missing sources; cases without materials.xml are reported per case or “(all cases)” if none had the file.

### Notes

- The script expects immediate child folders in spherical_cases to define the set of destination case names to process.
- Run from mat_extract to ensure paths resolve correctly; otherwise adjust repo_root logic if relocating the script.
