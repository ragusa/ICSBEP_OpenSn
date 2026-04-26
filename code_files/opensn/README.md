# OpenSn Script Generation

## Purpose

This folder auto-generates OpenSn Python run scripts for ICSBEP spherical benchmark cases using each case’s mesh and material data. Each generated script is ready to run a k-eigenvalue calculation and export scalar flux fields.

## How to Run

From this folder, run:
  - python3 [opensn_gen.py](./opensn_gen.py)
  - python3 [run_opensn_scripts.py](./orun_opensn_scripts.py)

## Inputs

- ../../spherical_cases tree discovered automatically from this folder.
- For each case directory:
- mesh/radii.txt listing shell radii (one per line).
- mesh/geometry.xml defining spherical cells and material IDs.
- materials/material_{ID}_{NAME}/ containing OpenMC HDF5 cross sections ({NAME}_LANL70g.h5).

## Outputs

- One OpenSn run script per case named {BENCHMARK}_{case-*}.py written into each case folder.
- Scripts load per-shell multi-group XS, apply volume-correction scaling, solve the k-eigenvalue problem, and export group-wise scalar flux to PVTU.
- result files names similarly to the created OpenSn scripts but with .pvtu and .vtu extensions
- A report of data from cases called [opensn_run_report.csv](./opensn_run_report.csv)

### Notes

- generate_opensn_scripts.py walks all benchmark subfolders under spherical_cases and processes any case-* directories it finds.
- Each spherical shell gets its own MultiGroupXS object with an independent scaling factor derived from exact versus mesh-computed shell volumes.
- Material folders must follow material_{number}_{name} naming and contain matching XS files, or the case is skipped with a warning.
- OPEN_SN needs to be changed according to your system's Open Sn installation