# Gmsh Mesh Generation

## Purpose

This folder generates Gmsh 3D meshes for spherical (n-shell) cases from the project’s mesh tree. Mesh files are saved in each case folder beside radii.txt. Failed cases are listed in [failed_cases.txt](./failed_cases.txt).

## How to Run

- From this folder, run:
  python3 [Create_ICSBEP_Meshes.py](./Create_ICSBEP_Meshes.py)

## Outputs

- n_shells_sphere_{N}_shells.msh files beside each radii.txt.
- [failed_cases.txt](./failed_cases.txt) listing any failed or timed-out cases.

## Notes

- [Create_ICSBEP_Meshes.py](./Create_ICSBEP_Meshes.py) discovers all cases and manages parallel processing of [nShell.py](./nShell.py).
- [nShell.py](./nShell.py) builds and meshes spherical geometries.
- Environment variables can adjust runtime, threading, and meshing parameters.
- If gmsh is missing, install the Python gmsh module.
- Invalid or non-increasing radii will cause validation errors and appear in failed_cases.txt.
