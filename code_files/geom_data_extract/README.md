# Geometry Extraction

## Extraction Process

Use the [spherical case extraction script](./duplicate_spherical_cases.py) to scan your OpenMC dataset and create a mirrored mesh tree containing only cases with exclusively spherical geometry. A [spherical_results.pkl](./spherical_results.pkl) file will be created with a list of result dictionaries.

## How to Run

- Run the script from this folder:
  python3 [duplicate_spherical_cases.py](./duplicate_spherical_cases.py)

## Outputs

- spherical_results.pkl in this folder with scan results.
- Mirrored mesh tree called [spherical_cases](../../spherical_cases/) with case folders.
- radii.txt within a mesh folder created in each case folder containing sorted sphere radii.

## Notes

- Only cases where every surface in geometry.xml is type "sphere" are included.
- The script processes directories at or below an "openmc" folder.
- Case folders are normalized to "case-N" based on the first case-like segment after "openmc".
- No third-party dependencies (Python standard library only).
