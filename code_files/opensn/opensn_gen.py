#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auto-generate OpenSn run scripts for ICSBEP spherical benchmark cases.

This script lives in the code_files/opensn/ folder and automatically
finds ../../spherical_cases. It walks the directory tree, finds each
case-* folder, reads its radii.txt, geometry.xml, and material sub-folders,
and writes a ready-to-run OpenSn Python script into each case folder.

Special handling:
- If a cell in geometry.xml has material="void", the generated OpenSn script
  will load the void multigroup XS from:
      ../../../code_files/mat_extract/void.xs

Usage:
python generate_opensn_scripts.py
"""

import os
import sys
import re
import xml.etree.ElementTree as ET


def parse_radii(radii_path: str) -> list[str]:
    """Read radii.txt and return a list of radii strings, preserving exact precision."""
    radii = []
    with open(radii_path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                radii.append(line)
    return radii


def parse_geometry_xml(geometry_path: str) -> list[str]:
    """
    Parse geometry.xml and return an ordered list of material references,
    one per spherical shell (cell), in the order the cells appear.

    Each <cell> element has a 'material' attribute. This may be a numeric
    material ID like "1" or a token like "void".
    """
    tree = ET.parse(geometry_path)
    root = tree.getroot()
    cells = root.findall("cell")
    cells_sorted = sorted(cells, key=lambda c: int(c.attrib["id"]))

    material_refs = []
    for cell in cells_sorted:
        if "material" not in cell.attrib:
            raise KeyError(f"Cell id={cell.attrib.get('id', '?')} missing 'material' attribute")
        material_refs.append(cell.attrib["material"].strip())

    return material_refs


def discover_materials(materials_dir: str) -> dict[int, str]:
    """
    Scan the materials/ directory for sub-folders named
    'material_{#}_{name}' and return a dict {mat_number: mat_name}.
    """
    mat_map = {}
    for entry in sorted(os.listdir(materials_dir)):
        entry_path = os.path.join(materials_dir, entry)
        if not os.path.isdir(entry_path):
            continue

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
    cell_material_refs: list[str],
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
    lines.append('rank = 0')
    lines.append('size = 1')
    lines.append('')
    lines.append('if "opensn_console" not in globals():')
    lines.append('    from mpi4py import MPI')
    lines.append('    size = MPI.COMM_WORLD.size')
    lines.append('    rank = MPI.COMM_WORLD.rank')
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
    lines.append('        partitioner=PETScGraphPartitioner(type="parmetis"),')
    lines.append('    )')
    lines.append('    grid = meshgen.Execute()')
    lines.append('')

    # ---- volumes -------------------------------------------------------
    lines.append('    # get "measured" volumes')
    lines.append('    volumes_per_block = grid.ComputeVolumePerBlockID()')
    lines.append('')

    # ---- radii dict ----------------------------------------------------
    radii_entries = ", ".join(f"{i + 1}: {r}" for i, r in enumerate(radii))
    lines.append('    # define radii per block-ID')
    lines.append(f'    radii = {{{radii_entries}}}')
    lines.append('')

    # ---- volume correction ---------------------------------------------
    lines.append('    # sort block IDs by increasing radius')
    lines.append('    block_ids = sorted(volumes_per_block.keys())')
    lines.append('')
    lines.append('    missing = [blk for blk in block_ids if blk not in radii]')
    lines.append('    if missing:')
    lines.append('        raise KeyError(f"No radius provided for block IDs: {missing}")')
    lines.append('')
    lines.append('    # compute exact shell volumes')
    lines.append('    exact_volumes_per_block = {}')
    lines.append('    prev_R = 0.0')
    lines.append('    for blk in block_ids:')
    lines.append('        R = radii[blk]')
    lines.append('        exact_volumes_per_block[blk] = (4.0 / 3.0) * np.pi * (R**3 - prev_R**3)')
    lines.append('        prev_R = R')
    lines.append('')
    lines.append('    # build scaling ratios')
    lines.append('    ratios = np.array([')
    lines.append('        exact_volumes_per_block[blk] / volumes_per_block[blk]')
    lines.append('        for blk in block_ids')
    lines.append('    ])')
    lines.append('')
    lines.append('    if rank == 0:')
    lines.append('        for idx, blk in enumerate(block_ids):')
    lines.append('            print(f"Block {blk}: measured = {volumes_per_block[blk]:.11e}, "')
    lines.append('                  f" exact = {exact_volumes_per_block[blk]:.11e}, "')
    lines.append('                  f" ratio = {ratios[idx]:.6f}")')
    lines.append('')

    # ---- void xs path --------------------------------------------------
    lines.append('    void_xs_path = os.path.abspath(')
    lines.append('        os.path.join(os.path.dirname(__file__), "../../../code_files/mat_extract/void.xs")')
    lines.append('    )')
    lines.append('')

    # ---- cross-section loading: one xs object PER SHELL ----------------
    lines.append('    # load XS one object per shell for independent scaling')
    for shell_idx, mat_ref in enumerate(cell_material_refs):
        shell_num = shell_idx + 1
        mat_ref_clean = mat_ref.strip()

        lines.append(f'    xs_shell{shell_num} = MultiGroupXS()')

        if mat_ref_clean.lower() == "void":
            lines.append(f'    xs_shell{shell_num}.LoadFromOpenSn(void_xs_path)')
            lines.append(f'    xs_shell{shell_num}.SetScalingFactor(ratios[{shell_idx}])')
            lines.append(f'    if rank == 0:')
            lines.append(f'        print("Shell {shell_num}: using void XS ->", void_xs_path)')
            lines.append('')
        else:
            try:
                mat_id = int(mat_ref_clean)
            except ValueError:
                raise ValueError(
                    f"Unsupported material reference '{mat_ref_clean}' in geometry.xml. "
                    "Expected integer material ID or 'void'."
                )

            if mat_id not in mat_map:
                raise KeyError(
                    f"Material ID {mat_id} referenced in geometry.xml "
                    "but no matching material folder was found."
                )

            mat_name = mat_map[mat_id]
            lines.append(
                f'    xs_shell{shell_num}.LoadFromOpenMC("./materials/material_{mat_id}_{mat_name}/{mat_name}_LANL70g.h5", "{mat_name}", 294.0)'
            )
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
    lines.append('    if rank == 0:')
    lines.append('        print(f"Computed k-eigenvalue: {k}")')
    lines.append('')

    # ---- export --------------------------------------------------------
    lines.append('    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)')
    lines.append(f'    vtk_basename = "{vtk_name}"')
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

    benchmark_name = benchmark_folder_name.upper().replace("-", "_")
    case_name = case_folder_name

    mesh_dir = os.path.join(case_dir, "mesh")
    materials_dir = os.path.join(case_dir, "materials")
    radii_path = os.path.join(mesh_dir, "radii.txt")
    geometry_path = os.path.join(mesh_dir, "geometry.xml")

    for p, desc in [
        (mesh_dir, "mesh directory"),
        (materials_dir, "materials directory"),
        (radii_path, "radii.txt"),
        (geometry_path, "geometry.xml"),
    ]:
        if not os.path.exists(p):
            print(f"WARNING: {desc} not found at {p}, skipping {case_dir}")
            return

    radii = parse_radii(radii_path)
    cell_material_refs = parse_geometry_xml(geometry_path)
    mat_map = discover_materials(materials_dir)

    # Validate material references in geometry.xml
    missing_mats = set()
    invalid_tokens = set()

    for ref in set(cell_material_refs):
        ref_clean = ref.strip()
        if ref_clean.lower() == "void":
            continue
        try:
            mat_id = int(ref_clean)
        except ValueError:
            invalid_tokens.add(ref_clean)
            continue

        if mat_id not in mat_map:
            missing_mats.add(mat_id)

    if invalid_tokens:
        print(
            f"WARNING: Unsupported material token(s) {sorted(invalid_tokens)} "
            f"in {geometry_path}. Expected integer IDs or 'void'. Skipping."
        )
        return

    if missing_mats:
        print(
            f"WARNING: Materials {sorted(missing_mats)} referenced in geometry.xml "
            f"but no matching folder found in {materials_dir}. Skipping."
        )
        return

    if len(cell_material_refs) != len(radii):
        print(
            f"WARNING: {len(cell_material_refs)} cells in geometry.xml but "
            f"{len(radii)} radii in radii.txt for {case_dir}. Skipping."
        )
        return

    script_content = generate_script(
        benchmark_name=benchmark_name,
        case_name=case_name,
        radii=radii,
        cell_material_refs=cell_material_refs,
        mat_map=mat_map,
    )

    script_filename = f"{benchmark_name}_{case_name}.py"
    script_path = os.path.join(case_dir, script_filename)
    with open(script_path, "w") as f:
        f.write(script_content)

    os.chmod(script_path, 0o755)
    print(f"Written: {script_path}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    spherical_cases_dir = os.path.normpath(os.path.join(script_dir, "../../spherical_cases"))

    if not os.path.isdir(spherical_cases_dir):
        print(f"ERROR: spherical_cases directory not found at {spherical_cases_dir}")
        print("Expected this script to be in code_files/opensn/")
        sys.exit(1)

    print(f"Scanning: {spherical_cases_dir}")

    benchmark_dirs = sorted(
        [
            os.path.join(spherical_cases_dir, d)
            for d in os.listdir(spherical_cases_dir)
            if os.path.isdir(os.path.join(spherical_cases_dir, d))
        ]
    )

    total_generated = 0
    for bench_dir in benchmark_dirs:
        bench_name = os.path.basename(bench_dir)
        case_dirs = sorted(
            [
                os.path.join(bench_dir, d)
                for d in os.listdir(bench_dir)
                if os.path.isdir(os.path.join(bench_dir, d)) and d.startswith("case")
            ]
        )
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