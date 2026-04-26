# ICSBEP_OpenSn

Tools to extract purely spherical ICSBEP cases, mirror them into a normalized case tree, generate 3D Gmsh meshes, and stage OpenMC materials for downstream OpenSn workflows.

### Repository layout

- [icsbep_original](./icsbep_original/)  Source tree containing the all ICSBEP criticality benchmark inputs (e.g., geometry.xml, materials.xml) to be extracted.
- [spherical_cases](./spherical_cases/)  Destination tree containing only spherical cases; Each case folder contains the generated mesh and materials.
- [code_files](./code_files/README.md) Folder containing all code needed for analyzing criticality benchmark cases.

    - [geom_data_extract](./code_files/geom_data_extract/README.md)  Scans [icsbep_original](./icsbep_original/) folder. For cases where every surface in geometry.xml is a sphere, creates a mirrored mesh tree within [spherical_cases](./spherical_cases/) with all necessary files separated into "case-#" folders.
    - [gmsh_code](./code_files/gmsh_code/README.md)  Builds gmsh geometries for each spherical case discovered from the mesh tree.
    - [mat_extract](./code_files/mat_extract/README.md)  Utilizes OpenMC on each materials.xml file in [spherical_cases](./spherical_cases/) to gather multigroup cross section data.
    - [opensn](./code_files/opensn/README.md)  For each completed case in [spherical_cases](./spherical_cases/), generates a python script to run OpenSn on previously generated mesh and geometry data.

