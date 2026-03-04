#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auto-generate OpenSn run scripts for ICSBEP spherical benchmark cases.

This script lives in the code_files/OpenSn/ folder and automatically
finds ../../spherical_cases. It walks the directory tree, finds each
case-* folder, reads its radii.txt, geometry.xml, and material sub-folders,
and writes a ready-to-run OpenSn Python script into each case folder.

Usage:
    python generate_opensn_scripts.py
"""

import os
import sys
import re
import xml.etree.ElementTree as ET
from pathlib import Path


def parse_radii(radii_path: str) -> list[str]:
    """Read radii.txt and return a list of radii strings (one per line), preserving exact precision."""
    radii = []
    with open(radii_path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                radii.append(line)
    return radii


def parse_geometry_xml(geometry_path: str) -> list[int]:
    """
    Parse geometry.xml and return an ordered list of material IDs,
    one per spherical shell (cell), in the order the cells appear.

    Each <cell> element has a 'material' attribute giving the material ID.
    The cells are ordered from the innermost shell outward.
    """
    tree = ET.parse(geometry_path)
    root = tree.getroot()
    cells = root.findall("cell")
    # Sort cells by their id attribute to guarantee inside-out order
    cells_sorted = sorted(cells, key=lambda c: int(c.attrib["id"]))
    material_ids = [int(c.attrib["material"]) for c in cells_sorted]
    return material_ids


def discover_materials(materials_dir: str) -> dict[int, str]:
    """
    Scan the materials/ directory for sub-folders named
    'material_{#}_{name}' and return a dict  {mat_number: mat_name}.
    """
    mat_map = {}
    for entry in sorted(os.listdir(materials_dir)):
        entry_path = os.path.join(materials_dir, entry)
        if not os.path.isdir(entry_path):
            continue
        # Pattern: material_{number}_{name}
        m = re.match(r"^material_(\d+)_(.+)$", entry)
        if m:
            mat_num = int(m.group(1))
            mat_name = m.group(2)
            mat_map[mat_num] = mat_name
    return mat_map


def generate_script(
    benchmark_name: str,
    case_name: str,
    radii: list[str],
    cell_material_ids: list[int],
    mat_map: dict[int, str],
) -> str:
    """
    Build the full OpenSn Python script as a string.

    Each spherical shell gets its own MultiGroupXS object so that every
    shell carries its own per-shell volume-correction scaling factor.
    """

    n_shells = len(radii)
    lines = []

    # ---- header --------------------------------------------------------
    lines.append('#!/usr/bin/env python3')
    lines.append('# -*- coding: utf-8 -*-')
    lines.append('"""')
    lines.append(f'{benchmark_name} {case_name} benchmark')
    lines.append('"""')
    lines.append('')

    # ---- imports -------------------------------------------------------
    lines.append('import sys')
    lines.append('import os')
    lines.append('import numpy as np')
    lines.append('')
    lines.append('if "opensn_console" not in globals():')
    lines.append('    from mpi4py import MPI')
    lines.append('    size = MPI.COMM_WORLD.size')
    lines.append('    rank = MPI.COMM_WORLD.rank')
    lines.append('    # Append parent directory to locate the pyopensn modules')
    lines.append('    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))')
    lines.append('    from pyopensn.mesh import FromFileMeshGenerator, PETScGraphPartitioner')
    lines.append('    from pyopensn.xs import MultiGroupXS')
    lines.append('    from pyopensn.aquad import GLCProductQuadrature3DXYZ')
    lines.append('    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver')
    lines.append('    from pyopensn.fieldfunc import FieldFunctionGridBased')
    lines.append('')

    # ---- main block ----------------------------------------------------
    lines.append('if __name__ == "__main__":')
    lines.append('')
    lines.append('    meshgen = FromFileMeshGenerator(')
    lines.append(f'        filename="./mesh/n_shells_sphere_{n_shells}_shells.msh",')
    lines.append("        partitioner=PETScGraphPartitioner(type='parmetis'),")
    lines.append('    )')
    lines.append('    grid = meshgen.Execute()')
    lines.append('')

    # ---- volumes -------------------------------------------------------
    lines.append('    # get "measured" volumes')
    lines.append('    volumes_per_block = grid.ComputeVolumePerBlockID()')
    lines.append('')

    # ---- radii dict ----------------------------------------------------
    radii_entries = ", ".join(f"{i+1}: {r}" for i, r in enumerate(radii))
    lines.append('    # define radii per block-ID')
    lines.append(f'    radii = {{{radii_entries}}}')
    lines.append('')

    # ---- volume correction ---------------------------------------------
    lines.append('    # sort block IDs by increasing radius')
    lines.append('    block_ids = sorted(volumes_per_block.keys())')
    lines.append('')
    lines.append('    # check:')
    lines.append('    missing = [blk for blk in block_ids if blk not in radii]')
    lines.append('    if missing:')
    lines.append('        raise KeyError(f"No radius provided for block IDs: {missing}")')
    lines.append('')
    lines.append('    # compute exact volumes')
    lines.append('    exact_volumes_per_block = {}')
    lines.append('    prev_R = 0.0')
    lines.append('    for blk in block_ids:')
    lines.append('        R = radii[blk]')
    lines.append('        # shell from prev_R to R')
    lines.append('        exact_volumes_per_block[blk] = (4.0 / 3.0) * np.pi * (R**3 - prev_R**3)')
    lines.append('        prev_R = R')
    lines.append('')
    lines.append('    # build the ratios array')
    lines.append('    ratios = np.array([')
    lines.append('        exact_volumes_per_block[blk] / volumes_per_block[blk]')
    lines.append('        for blk in block_ids')
    lines.append('    ])')
    lines.append('')
    lines.append('    # only rank 0 prints')
    lines.append('    if rank == 0:')
    lines.append('        for idx, blk in enumerate(block_ids):')
    lines.append('            print(f"Block {blk}: measured = {volumes_per_block[blk]:.11e}, "')
    lines.append('                  f" exact = {exact_volumes_per_block[blk]:.11e}, "')
    lines.append('                  f" ratio = {ratios[idx]:.6f}")')
    lines.append('')

    # ---- cross-section loading: one xs object PER SHELL ----------------
    lines.append('    # load XS  one object per shell for independent scaling')
    for shell_idx, mat_id in enumerate(cell_material_ids):
        shell_num = shell_idx + 1
        mat_name = mat_map[mat_id]
        lines.append(f'    xs_shell{shell_num} = MultiGroupXS()')
        lines.append(f'    xs_shell{shell_num}.LoadFromOpenMC("./materials/material_{mat_id}_{mat_name}/{mat_name}_LANL70g.h5", "{mat_name}", 294.0)')
        lines.append(f'    xs_shell{shell_num}.SetScalingFactor(ratios[{shell_idx}])')
        lines.append('')

    # ---- num_groups ----------------------------------------------------
    lines.append('    num_groups = xs_shell1.num_groups')
    lines.append('    if rank == 0:')
    lines.append('        print("num groups =", num_groups)')
    lines.append('')

    # ---- solver setup --------------------------------------------------
    lines.append('    # Solver')
    lines.append('    phys = DiscreteOrdinatesProblem(')
    lines.append('        mesh=grid,')
    lines.append('        num_groups=num_groups,')
    lines.append('        groupsets=[')
    lines.append('            {')
    lines.append('                "groups_from_to": (0, num_groups - 1),')
    lines.append('                "angular_quadrature": GLCProductQuadrature3DXYZ(')
    lines.append('                    n_polar=8,')
    lines.append('                    n_azimuthal=16,')
    lines.append('                    scattering_order=3')
    lines.append('                ),')
    lines.append('                "inner_linear_method": "petsc_gmres",')
    lines.append('                "angle_aggregation_type": "single",')
    lines.append('                "angle_aggregation_num_subsets": 1,')
    lines.append('                "l_max_its": 500,')
    lines.append('                "l_abs_tol": 1.0e-6,')
    lines.append('            },')
    lines.append('        ],')
    lines.append('        xs_map=[')
    for shell_idx in range(n_shells):
        block_id = shell_idx + 1
        lines.append(f'            {{"block_ids": [{block_id}], "xs": xs_shell{block_id}}},')
    lines.append('        ],')
    lines.append('        options={')
    lines.append('            "use_precursors": False,')
    lines.append('            "verbose_inner_iterations": True,')
    lines.append('            "verbose_outer_iterations": True,')
    lines.append('        },')
    lines.append('    )')
    lines.append('')

    # ---- k-eigenvalue solver -------------------------------------------
    vtk_name = f"{benchmark_name}_{case_name}"
    lines.append('    k_solver = NonLinearKEigenSolver(')
    lines.append('        problem=phys,')
    lines.append('        nl_max_its=500,')
    lines.append('        nl_abs_tol=1.0e-10,')
    lines.append('    )')
    lines.append('    k_solver.Initialize()')
    lines.append('    k_solver.Execute()')
    lines.append('    k = k_solver.GetEigenvalue()')
    lines.append('    # only rank 0 prints')
    lines.append('    if rank == 0:')
    lines.append('        print(f"Computed k-eigenvalue: {k}")')
    lines.append('')

    # ---- export --------------------------------------------------------
    lines.append('    # export')
    lines.append('    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)')
    lines.append(f'    vtk_basename = "{vtk_name}"')
    lines.append('    # export only the flux of group g (first []), moment 0 (second [])')
    lines.append('    FieldFunctionGridBased.ExportMultipleToPVTU(')
    lines.append('        [fflist[g][0] for g in range(num_groups)],')
    lines.append('        vtk_basename')
    lines.append('    )')
    lines.append('')

    return "\n".join(lines)


def process_case(benchmark_dir: str, case_dir: str):
    """
    Process a single case directory and write the OpenSn script.
    """
    benchmark_folder_name = os.path.basename(benchmark_dir)
    case_folder_name = os.path.basename(case_dir)

    # Convert folder name to benchmark name (e.g. heu-met-fast-001 -> HEU_MET_FAST_001)
    benchmark_name = benchmark_folder_name.upper().replace("-", "_")
    case_name = case_folder_name

    # Paths
    mesh_dir = os.path.join(case_dir, "mesh")
    materials_dir = os.path.join(case_dir, "materials")
    radii_path = os.path.join(mesh_dir, "radii.txt")
    geometry_path = os.path.join(mesh_dir, "geometry.xml")

    # Validate required files/dirs exist
    for p, desc in [
        (mesh_dir, "mesh directory"),
        (materials_dir, "materials directory"),
        (radii_path, "radii.txt"),
        (geometry_path, "geometry.xml"),
    ]:
        if not os.path.exists(p):
            print(f"  WARNING: {desc} not found at {p}, skipping {case_dir}")
            return

    # Parse inputs
    radii = parse_radii(radii_path)
    cell_material_ids = parse_geometry_xml(geometry_path)
    mat_map = discover_materials(materials_dir)

    # Validate that all materials referenced in geometry exist as folders
    unique_referenced = set(cell_material_ids)
    missing_mats = unique_referenced - set(mat_map.keys())
    if missing_mats:
        print(f"  WARNING: Materials {missing_mats} referenced in geometry.xml "
              f"but no matching folder found in {materials_dir}. Skipping.")
        return

    # Sanity check: number of shells == number of radii
    if len(cell_material_ids) != len(radii):
        print(f"  WARNING: {len(cell_material_ids)} cells in geometry.xml but "
              f"{len(radii)} radii in radii.txt for {case_dir}. Skipping.")
        return

    # Generate script
    script_content = generate_script(
        benchmark_name=benchmark_name,
        case_name=case_name,
        radii=radii,
        cell_material_ids=cell_material_ids,
        mat_map=mat_map,
    )

    # Write script to case directory
    script_filename = f"{benchmark_name}_{case_name}.py"
    script_path = os.path.join(case_dir, script_filename)
    with open(script_path, "w") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    print(f"  Written: {script_path}")


def main():
    # Resolve spherical_cases relative to this script's location (code_files/OpenSn/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    spherical_cases_dir = os.path.normpath(os.path.join(script_dir, "../../spherical_cases"))

    if not os.path.isdir(spherical_cases_dir):
        print(f"ERROR: spherical_cases directory not found at {spherical_cases_dir}")
        print("  Expected this script to be in code_files/OpenSn/")
        sys.exit(1)

    print(f"Scanning: {spherical_cases_dir}")

    # Walk: spherical_cases / benchmark / case-*
    benchmark_dirs = sorted([
        os.path.join(spherical_cases_dir, d)
        for d in os.listdir(spherical_cases_dir)
        if os.path.isdir(os.path.join(spherical_cases_dir, d))
    ])

    total_generated = 0
    for bench_dir in benchmark_dirs:
        bench_name = os.path.basename(bench_dir)
        case_dirs = sorted([
            os.path.join(bench_dir, d)
            for d in os.listdir(bench_dir)
            if os.path.isdir(os.path.join(bench_dir, d)) and d.startswith("case")
        ])
        if not case_dirs:
            continue
        print(f"\nBenchmark: {bench_name}")
        for case_dir in case_dirs:
            print(f"  Processing: {os.path.basename(case_dir)}")
            process_case(bench_dir, case_dir)
            total_generated += 1

    print(f"\nDone. Processed {total_generated} case(s).")


if __name__ == "__main__":
    main()