#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HEU_SOL_THERM_004 case-4 benchmark
"""

import sys
import os
import numpy as np

rank = 0
size = 1

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, PETScGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(
        filename="./mesh/n_shells_sphere_4_shells.msh",
        partitioner=PETScGraphPartitioner(type="parmetis"),
    )
    grid = meshgen.Execute()

    # get "measured" volumes
    volumes_per_block = grid.ComputeVolumePerBlockID()

    # define radii per block-ID
    radii = {1: 20.873999999999999, 2: 20.975999999999999, 3: 44.344000000000001, 4: 44.597999999999999}

    # sort block IDs by increasing radius
    block_ids = sorted(volumes_per_block.keys())

    missing = [blk for blk in block_ids if blk not in radii]
    if missing:
        raise KeyError(f"No radius provided for block IDs: {missing}")

    # compute exact shell volumes
    exact_volumes_per_block = {}
    prev_R = 0.0
    for blk in block_ids:
        R = radii[blk]
        exact_volumes_per_block[blk] = (4.0 / 3.0) * np.pi * (R**3 - prev_R**3)
        prev_R = R

    # build scaling ratios
    ratios = np.array([
        exact_volumes_per_block[blk] / volumes_per_block[blk]
        for blk in block_ids
    ])

    if rank == 0:
        for idx, blk in enumerate(block_ids):
            print(f"Block {blk}: measured = {volumes_per_block[blk]:.11e}, "
                  f" exact = {exact_volumes_per_block[blk]:.11e}, "
                  f" ratio = {ratios[idx]:.6f}")

    void_xs_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../../code_files/mat_extract/void.xs")
    )

    # load XS one object per shell for independent scaling
    xs_shell1 = MultiGroupXS()
    xs_shell1.LoadFromOpenMC("./materials/material_1_UO2F2_D2O_Solution/UO2F2_D2O_Solution_LANL70g.h5", "UO2F2_D2O_Solution", 294.0)
    xs_shell1.SetScalingFactor(ratios[0])

    xs_shell2 = MultiGroupXS()
    xs_shell2.LoadFromOpenMC("./materials/material_2_321_Stainless_Steel/321_Stainless_Steel_LANL70g.h5", "321_Stainless_Steel", 294.0)
    xs_shell2.SetScalingFactor(ratios[1])

    xs_shell3 = MultiGroupXS()
    xs_shell3.LoadFromOpenMC("./materials/material_3_Heavy_Water/Heavy_Water_LANL70g.h5", "Heavy_Water", 294.0)
    xs_shell3.SetScalingFactor(ratios[2])

    xs_shell4 = MultiGroupXS()
    xs_shell4.LoadFromOpenMC("./materials/material_2_321_Stainless_Steel/321_Stainless_Steel_LANL70g.h5", "321_Stainless_Steel", 294.0)
    xs_shell4.SetScalingFactor(ratios[3])

    num_groups = xs_shell1.num_groups
    if rank == 0:
        print("num groups =", num_groups)

    # Solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=8,
                    n_azimuthal=16,
                    scattering_order=3
                ),
                "inner_linear_method": "petsc_gmres",
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "l_max_its": 500,
                "l_abs_tol": 1.0e-6,
            },
        ],
        xs_map=[
            {"block_ids": [1], "xs": xs_shell1},
            {"block_ids": [2], "xs": xs_shell2},
            {"block_ids": [3], "xs": xs_shell3},
            {"block_ids": [4], "xs": xs_shell4},
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": True,
            "verbose_outer_iterations": True,
        },
    )

    k_solver = NonLinearKEigenSolver(
        problem=phys,
        nl_max_its=500,
        nl_abs_tol=1.0e-10,
    )
    k_solver.Initialize()
    k_solver.Execute()
    k = k_solver.GetEigenvalue()
    if rank == 0:
        print(f"Computed k-eigenvalue: {k}")

    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)
    vtk_basename = "HEU_SOL_THERM_004_case-4"
    FieldFunctionGridBased.ExportMultipleToPVTU(
        [fflist[g][0] for g in range(num_groups)],
        vtk_basename
    )
