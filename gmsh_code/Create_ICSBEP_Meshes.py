import os
from pathlib import Path
import numpy as np
import gmsh 
import sys, subprocess
from pathlib import Path

def run_worker(radii: str,data_path: str) -> None:
    subprocess.run([sys.executable, "InputVerification.py", radii, data_path], check=True)


# Navigate up one directory
base_dir = Path(__file__).parent.parent
# Go into "mesh" folder
mesh_dir = base_dir / "mesh"
filepaths = []
all_radii = []
# Iterate through all folders inside "mesh"
for dirpath, _, filenames in os.walk(mesh_dir):
    for filename in filenames:
        if filename == "radii.txt":
            filepath = Path(dirpath) / filename
            filepaths.append(filepath)
            with open(filepath, 'r') as file:
                radii = [line.strip() for line in file if line.strip()]
                all_radii.append(radii)
                


# Example: print filepaths and corresponding radii

# #k = 1
filepaths = filepaths[2:3]
all_radii = all_radii[2:3]
for fp, radii_list in zip(filepaths, all_radii):
    run_worker(str(radii_list),str(fp))
#     print(k)
#     print(fp)
#     print(radii_list)
#     #if len(radii_list) > 2:
#     #    k+=1
#     #    continue
#     # Initalize gmsh
#     gmsh.initialize(sys.argv)
#     gmsh.model.add("Sphere" + str(k))
#     k += 1
#     radii = [float(x) for x in radii_list]
#     radii.sort()
#     i = 1
#     sphere = [0]*len(radii)
#     for r in radii:
#         # Sphere
#         origin = 3*[0.0]
#         sphere[i-1] = gmsh.model.occ.addSphere(origin[0], origin[1], origin[2],
#                                     r*2)
#         # By default removeObject and removeTool are True, 
#         # we do not want to lose the inner box
#         if (i>1):
#             gmsh.model.occ.cut([(3,sphere[i-1])],
#                                [(3,sphere[i-2])],
#                                removeTool=False)
#         gmsh.model.occ.synchronize()
#         i+=1
    
#     for j in range(i-1):
#          gmsh.model.addPhysicalGroup(dim=3,tags=[sphere[j]],
#                                      name="Sphere" + str(j+1))
#     gmsh.model.occ.synchronize()
#     gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
#     gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
#     gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
#     gmsh.option.setNumber("Mesh.MeshSizeMax", 10)
#     # Export Mesh to specified version for OpenSn
#     gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    
#     gmsh.option.setNumber("Mesh.Algorithm", 6)
#     # Fast algorithm
#     # gmsh.option.setNumber("Mesh.Algorithm3D", 10)
#     # gmsh.model.mesh.generate(2)
#     gmsh.model.mesh.generate(3)
#     gmsh.write(str(fp.with_name("Spheres.msh")))
    
#     #if '-nopopup' not in sys.argv:
#     #    gmsh.fltk.run()
#     gmsh.finalize()
#     print("finish")