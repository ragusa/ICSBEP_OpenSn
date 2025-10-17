# Gmsh code

## Overview

This folder generates 3D Gmsh meshes for spherical (n‑shell) cases found under the project’s mesh tree, writing binary .msh v4.1 files with physical groups.
It runs as‑is with no edits required and places failures in [failed\_cases.txt](./failed_cases.txt).

## Prerequisites

* Python 3 with the gmsh Python module available.
* Upstream geometry extraction can be done with the [geometry extraction script](../geom_data_extract/duplicate_folder_struct_for_mesh.py), which populates per‑case radii.txt files.

## Folder contents

* [Create\_ICSBEP\_Meshes.py](./Create_ICSBEP_Meshes.py)  orchestrates case discovery, parallel execution, timeouts, and reporting.
* [nShell.py](./nShell.py)  builds concentric spheres, assigns physical groups, configures mesh‑size fields, preflights, meshes, and writes the .msh.
* [failed\_cases.txt](./failed_cases.txt)  created/overwritten each run with any case that timed out or errored.

## Quick start

* From this folder, run:
  python3 [Create\_ICSBEP\_Meshes.py](./Create_ICSBEP_Meshes.py)
* Mesh files are written next to each per‑case radii.txt as n\_shells\_sphere\_{N}\_shells.msh

## Details

### [Create\_ICSBEP\_Meshes.py](./Create_ICSBEP_Meshes.py)

* Walks the mesh/ tree to collect every radii.txt, then launches tasks via ProcessPoolExecutor.
* Each task invokes [nShell.py](./nShell.py) with two arguments: the radii list (as a Python literal string) and the path to that case’s radii.txt.
* Per‑case timeout is set by CASE\_TIMEOUT\_SEC (defaults to 3600s), and parallelism defaults to CPU count or MAX\_WORKERS.
* Results print to the console; any failure is appended to [failed\_cases.txt](./failed_cases.txt) with case path, radii, and reason.
* Note: the current code intentionally limits to the first discovered case via slicing; keep it as‑is if that is the desired behavior.

### [nShell.py](./nShell.py)

* Reads radii, verifies positivity and strict increase, and derives the output directory from the provided data path.
* Builds geometry with occ.addSphere for each radius and uses Boolean cuts to form disjoint shells, then heals/deduplicates the model.
* Creates physical volume groups named Inner, Shell1, Shell2, … to label regions.
* Configures robust meshing defaults and a background field assembled from MathEval, Threshold, Min, and Max to control element sizes near interfaces.
* Runs an optional preflight that generates 1D/2D meshes to estimate node counts and adapts Mesh.MeshSizeFactor before 3D.
* Generates 1D/2D/3D, falls back to an alternative 3D algorithm on error, and writes n\_shells\_sphere\_{N}\_shells.msh (binary v4.1) to the case directory.

## Environment variables (optional)

* Orchestration: MAX\_WORKERS, CASE\_TIMEOUT\_SEC.
* Gmsh threading: GMESH\_THREADS, GMESH\_MAX\_THREADS\_2D, GMESH\_MAX\_THREADS\_3D.
* Geometry robustness: GEOMETRY\_TOL.
* Meshing detail: GMESH\_MIN\_CIRCLE\_POINTS, GMESH\_RANDOM\_FACTOR, GMESH\_N\_THICK, GMESH\_N\_CIRC\_FAR, GMESH\_N\_CIRC\_NEAR, GMESH\_ULTRA\_THIN\_RATIO, ABS\_THIN.
* Preflight behavior: GMESH\_SKIP\_PREFLIGHT, GMESH\_PREFLIGHT\_SEC, GMESH\_BUDGET\_2D\_NODES, GMESH\_SCALE\_INIT, GMESH\_SCALE\_MAX, GMESH\_SCALE\_ITERS, GMESH\_PANIC\_TRIG.

## Outputs

* One .msh per processed case located beside that case’s radii.txt, named n\_shells\_sphere\_{N}\_shells.msh.
* A fresh [failed\_cases.txt](./failed_cases.txt) in this folder on every run (empty if no failures).

## Troubleshooting

* If the gmsh module is missing, ensure it is installed in the active Python environment and visible to this script.
* If a case appears in [failed\_cases.txt](./failed_cases.txt), review the reason and re‑run with adjusted environment variables or fewer workers.
* Ensure radii are strictly increasing; malformed radii.txt will raise a validation error.
