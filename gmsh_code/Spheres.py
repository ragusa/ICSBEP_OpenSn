import numpy as np
import gmsh 
import sys

# Initalize gmsh
gmsh.initialize(sys.argv)
gmsh.model.add("Spheres")

# External Sphere
origin = 3*[0.0]
dx1 = 2.5
sphere1 = gmsh.model.occ.addSphere(origin[0], origin[1], origin[2],
                                dx1, tag=1)

# Internal Sphere
point2 = 3*[0.25]
dx2 = 2.0
sphere2 = gmsh.model.occ.addSphere(origin[0], origin[1], origin[2],
                                dx2, tag=2)

# By default removeObject and removeTool are True, 
# we do not want to lose the inner box
gmsh.model.occ.cut([(3,sphere1)],
                   [(3,sphere2)],
                   removeTool=False)

gmsh.model.occ.synchronize()

gmsh.model.addPhysicalGroup(dim=3,
                            tags=[1],
                            name="Sphere1")
gmsh.model.addPhysicalGroup(dim=3,
                            tags=[2],
                            name="Sphere2")

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)

# Export Mesh to specified version for OpenSn
gmsh.option.setNumber("Mesh.MshFileVersion",2.2)

gmsh.option.setNumber("Mesh.Algorithm", 6)
# Fast algorithm
# gmsh.option.setNumber("Mesh.Algorithm3D", 10)

# gmsh.model.mesh.generate(2)
gmsh.model.mesh.generate(3)
gmsh.write("Spheres.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
