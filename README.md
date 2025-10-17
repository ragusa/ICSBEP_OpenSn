# ICSBEP_OpenSn

Tools to extract purely spherical ICSBEP/OpenMC cases, mirror them into a normalized case tree, generate 3D Gmsh meshes, and stage OpenMC materials for downstream OpenSn workflows.

### Repository layout

- [geom_data_extract](./code_files/geom_data_extract/README.md)  Scans an OpenMC dataset, keeps only cases where every surface in geometry.xml is a sphere, and creates a mirrored mesh tree with normalized name/case-N folders. Produces one sorted, unique radii.txt per case and a summary spherical_results.pkl.
- [gmsh_code](./code_files/gmsh_code/README.md)  Builds concentric-sphere geometries and writes binary .msh v4.1 meshes with physical groups for each spherical case discovered from the mesh tree. Orchestrates parallel meshing and writes any failures to failed_cases.txt.
- [mat_extract](./code_files/mat_extract/README.md)  Copies OpenMC materials.xml from ICSBEP originals into spherical_cases under name/case-N/materials/materials.xml, creating the destination folders as needed.
- [icsbep_original](./icsbep_original/)  Source tree containing the ICSBEP/OpenMC criticality benchmark inputs (e.g., geometry.xml, materials.xml) used by the extractors.
- [spherical_cases](./spherical_cases/)  Destination tree containing the normalized spherical cases; after running the pipelines it holds per-case radii.txt, generated meshes, and staged materials.

### End-to-end workflow

1) Discover spherical cases  
   Run duplicate_folder_struct_for_mesh.py from [geom_data_extract](./code_files/geom_data_extract/README.md) to mirror only purely spherical cases and write radii.txt per case; a spherical_results.pkl summary is stored in the extractor folder.

2) Mesh the cases  
   Run Create_ICSBEP_Meshes.py from [gmsh_code](./code_files/gmsh_code/README.md) to generate one nshells_sphere_Nshells.msh per case (binary v4.1) next to its radii.txt; failures or timeouts are recorded in failed_cases.txt.

3) Stage materials  
   Run extract_material.py from [mat_extract](./code_files/mat_extract/README.md) to copy each case’s materials.xml from icsbep_original into spherical_cases/name/case-N/materials/materials.xml.

### What each code folder contains

- [geom_data_extract](./code_files/geom_data_extract/README.md)
  - duplicate_folder_struct_for_mesh.py: walks the OpenMC dataset, filters to sphere-only geometries, normalizes output case folders, and writes per-case radii.txt and top-level spherical_results.pkl.

- [gmsh_code](./code_files/gmsh_code/README.md)
  - Create_ICSBEP_Meshes.py: discovers case radii.txt files, runs meshing in parallel with per-case timeouts, and logs any failure lines to failed_cases.txt.
  - nShell.py: constructs concentric shells, assigns physical volume groups (Inner, Shell1, …), configures meshing fields, and writes .msh v4.1 next to each case’s radii.txt.

- [mat_extract](./code_files/mat_extract/README.md)
  - extract_material.py: resolves repo-relative paths, discovers destination case names from spherical_cases, and copies source materials.xml from icsbep_original into each case’s materials/ folder (preserving metadata). Helper scripts like load_single_material_from_xml.py and template_homogenous_medium.py are included for material handling patterns.

### Requirements

- Python 3 for all scripts.  
- Gmsh Python module required to run the meshing pipeline in [gmsh_code](./code_files/gmsh_code/README.md).  
- The geometry and material extractors rely on the Python standard library.

### Outputs

- spherical_results.pkl in the geometry extractor folder summarizing discovered spherical cases.  
- radii.txt in each mirrored case folder containing positive, strictly increasing radii for meshing.  
- One nshells_sphere_Nshells.msh per processed case (binary v4.1 with physical groups), written alongside that case’s radii.txt.  
- failed_cases.txt in gmsh_code recording any case that timed out or errored.  
- materials/materials.xml created per case under spherical_cases after running the material extractor.

### Troubleshooting

- Ensure radii in radii.txt are positive and strictly increasing; malformed inputs will be rejected by the meshing stage.  
- If a case appears in failed_cases.txt, review the reason and rerun with adjusted environment variables or fewer workers from [gmsh_code](./code_files/gmsh_code/README.md).  
- If materials are missing, verify the expected icsbep_original layout (name/openmc/case-*/materials.xml) and rerun from [mat_extract](./code_files/mat_extract/README.md).

### Quick links

- [geom_data_extract](./code_files/geom_data_extract/README.md)  
- [gmsh_code](./code_files/gmsh_code/README.md)  
- [mat_extract](./code_files/mat_extract/README.md)  
- [icsbep_original](./icsbep_original/)  
- [spherical_cases](./spherical_cases/)
