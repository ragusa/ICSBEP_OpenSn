# -*- coding: utf-8 -*-

"""

Generalized N-shell concentric spheres with Gmsh (OCC).

- Radii must be strictly increasing: r[0] < r[1] < ... < r[N]

- Core is the smallest sphere (r[0]); Shell k is between r[k-1] and r[k].

- Boolean cuts: for k = N..1, cut sphere r[k] by sphere r[k-1] with removeTool=False.

- Physical groups: tag=1 -> "Inner", then tag=2.."N+1" -> "Shell1".."ShellN".

"""

from ast import literal_eval

import gmsh

import sys

import os

 

# --------------------- User parameters ---------------------

# Edit this list to change the number and sizes of shells:

#radii = [1.0/3.0, 2.0/3.0, 1.0]   # => 2 shells (Shell1: r0->r1, Shell2: r1->r2)

#radii = [0.20, 0.40, 0.60, 0.80, 1.00, 2]  # example with 4 shells

radlist = sys.argv[1]
data_path = sys.argv[2]
#data_path = "/home/rwszolek3/research/ICSBEP_OpenSn/mesh/heu-met-fast-029/"
lst = literal_eval(radlist)
radii = [float(x) for x in lst] 
#radii = [1.91, 3.15, 4.01, 4.66, 5.35, 6.0, 6.75, 7.55, 8.35, 9.15, 11.0, 12.25]
uniform_size = 0.15        # global target mesh size

model_name   = "n_shells_sphere"  # Gmsh model name

out_dir      = data_path[:-9]      # output directory

# -----------------------------------------------------------

 

# ---------------------- Basic checks -----------------------

if len(radii) < 2:

    raise ValueError("Need at least two radii to form shells.")

if any(r <= 0.0 for r in radii):

    raise ValueError("All radii must be positive.")

if any(r2 <= r1 for r1, r2 in zip(radii, radii[1:])):

    raise ValueError("Radii must be strictly increasing.")

 

N = len(radii) - 1  # number of shells

# -----------------------------------------------------------

 

gmsh.initialize()

gmsh.model.add(model_name)

 

# Optional: print messages in terminal (useful when running without GUI)

gmsh.option.setNumber("General.Terminal", 1)

 

# ------------------- Create concentric spheres -------------------

# Keep track of the Gmsh volume tags for the primitive spheres (by radius index)

sphere_tags = []

for r in radii:

    tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, r)

    sphere_tags.append(tag)

 

# ------------------- Boolean cuts to form shells -----------------

# We reproduce the original logic:

#   For k = N down to 1:

#     Shell k = (sphere r[k]) CUT (sphere r[k-1]), with removeTool=False

#   - This removes the object (outer sphere) and keeps the tool (inner sphere)

#   - The smallest sphere r[0] remains as the "Inner" volume at the end

shell_tags = [None] * (N + 1)  # 1..N will be used; index 0 unused

 

for k in range(N, 0, -1):

    # Object: outer sphere r[k]; Tool: inner neighbor r[k-1]

    out_dimtags, _ = gmsh.model.occ.cut(

        [(3, sphere_tags[k])],

        [(3, sphere_tags[k - 1])],

        removeTool=False  # keep the inner sphere for the next cut

        # removeObject defaults to True, which we want (remove the outer sphere)

    )

    gmsh.model.occ.synchronize()

 

    # Expect exactly one new volume (the shell). Filter for 3D volumes just in case.

    new_vols = [(d, t) for (d, t) in out_dimtags if d == 3]

    if not new_vols:

        gmsh.finalize()

        raise RuntimeError(f"Cut failed to produce a shell for k={k} (r[{k-1}] -> r[{k}]).")

 

    # Take the first (and only) shell volume tag

    shell_tags[k] = new_vols[0][1]

 

# After all cuts:

# - sphere_tags[0] is the remaining inner solid (the core),

# - shell_tags[1..N] are the concentric shells outward.

 

# ----------------------- Physical groups -------------------------

# tag=1 -> Inner (core)

gmsh.model.addPhysicalGroup(3, [sphere_tags[0]], tag=1)

gmsh.model.setPhysicalName(3, 1, "Inner")

 

# tag=2..N+1 -> Shell1..ShellN

for k in range(1, N + 1):

    phys_tag = k + 1

    gmsh.model.addPhysicalGroup(3, [shell_tags[k]], tag=phys_tag)

    gmsh.model.setPhysicalName(3, phys_tag, f"Shell{k}")

 

# --------------------- Proportional mesh sizing ------------------
# Scale mesh size with radius: h(r) = (uniform_size / r_min) * r, clamped to [h_min, h_max]
r_min = radii[0]
scale = uniform_size / r_min     # ensures h(r_min) = uniform_size
h_min = 0.5 * uniform_size       # slightly finer allowed near core
h_max = 8.0 * uniform_size       # coarser allowed outward

# Define background field using only MathEval + Min/Max (no Constant fields)
fld = gmsh.model.mesh.field
# h(r) = scale * sqrt(x^2 + y^2 + z^2)
fld.add("MathEval", 1)
fld.setString(1, "F", f"{scale}*sqrt(x*x + y*y + z*z)")

# Lower bound as constant MathEval: h_min
fld.add("MathEval", 2)
fld.setString(2, "F", f"{h_min}")

# max(h(r), h_min)
fld.add("Max", 3)
fld.setNumbers(3, "FieldsList", [1, 2])

# Upper bound as constant MathEval: h_max
fld.add("MathEval", 4)
fld.setString(4, "F", f"{h_max}")

# min(max(h(r), h_min), h_max)
fld.add("Min", 5)
fld.setNumbers(5, "FieldsList", [3, 4])

# Use as background mesh field
fld.setAsBackgroundMesh(5)

# Make background field authoritative (recommended with background fields)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# (Optional) keep hard clamps consistent with the background field
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h_min)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h_max)


 

# ----------------------- Mesh generation -------------------------

gmsh.model.mesh.generate(3)

 

# -------------------------- Save mesh ----------------------------

os.makedirs(out_dir, exist_ok=True)

outfile = os.path.join(out_dir, f"{model_name}_{N}_shells.msh")

gmsh.write(outfile)

print(f"Mesh written to: {outfile}")

 

# --------------------------- Optional GUI ------------------------

if "-nopopup" not in sys.argv:

    gmsh.fltk.run()

 

gmsh.finalize()