# code_files

Folder neatly containing all code required for modeling criticality benchmarks

## End-to-end workflow

1) Discover spherical cases  
   Run [duplicate_folder_struct_for_mesh.py](./geom_data_extract/duplicate_folder_struct_for_mesh.py) from [geom_data_extract](./geom_data_extract/README.md) to mirror only purely spherical cases and write radii.txt per case.

2) Mesh the cases  
   Run [Create_ICSBEP_Meshes.py](./gmsh_code/Create_ICSBEP_Meshes.py) from [gmsh_code](./gmsh_code/README.md) to generate a mesh for each found case.

3) Stage materials  
   Run [extract_material.py](./mat_extract/extract_material.py) from [mat_extract](./mat_extract/README.md) to generate a material file for each mesh.


## Troubleshooting

- Ensure radii in radii.txt are positive and strictly increasing; malformed inputs will be rejected by the meshing stage.  
- If a case appears in failed_cases.txt, review the reason and rerun with adjusted environment variables or fewer workers from [gmsh_code](./gmsh_code/README.md).  
- If materials are missing, verify openmc installation and cross section libraries are installed.
