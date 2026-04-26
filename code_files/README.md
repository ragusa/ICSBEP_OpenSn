# Python Code Files

Folder neatly containing all code required for modeling criticality benchmarks

## End-to-end workflow

1) Discover and compile spherical cases  
   Run [duplicate_spherical_cases.py](./geom_data_extract/duplicate_spherical_cases.py) from [geom_data_extract](./geom_data_extract/README.md) to create a organized folder of only purely spherical cases and their needed data.

2) Mesh the cases  
   Run [Create_ICSBEP_Meshes.py](./gmsh_code/Create_ICSBEP_Meshes.py) from [gmsh_code](./gmsh_code/README.md) to generate a mesh for each case.

3) Generate material cross sections
   Run [generate_material_mgxs.py](./mat_extract/generate_material_mgxs.py) from [mat_extract](./mat_extract/README.md) to generate multigroup cross section data for each case.

4) Produce OpenSn scripts 
   Run [opensn_gen.py](./opensn/opensn_gen.py) from [opensn](./opensn/README.md) to generate python scripts that can be run with OpenSn to get results for each case.
   Run [run_opensn_scripts.py](./opensn/orun_opensn_scripts.py) from [opensn](./opensn/README.md) to run OpenSn using the generated scripts.

## Troubleshooting

- Ensure radii in radii.txt are positive and strictly increasing; malformed inputs will be rejected by the meshing stage.  
- If a case appears in failed_cases.txt, review the reason and rerun with adjusted environment variables or fewer workers from [gmsh_code](./gmsh_code/README.md).  
- If materials are missing, verify openmc installation and cross section libraries are installed.