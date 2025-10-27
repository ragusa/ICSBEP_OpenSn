# Import necessary modules
from ast import literal_eval       # Safely converts string representations of lists to actual Python lists
import gmsh, sys, os, math, time   # Gmsh for mesh operations, and standard libraries for system and math utilities

# --------------------- Parse Command-Line Arguments ---------------------
radlist = sys.argv[1]   # First command-line argument: list of radii (string or comma-separated)
data_path = sys.argv[2] # Second argument: path to data/output file

# Convert radius input string into a list of floats
lst = literal_eval(radlist) if radlist.strip().startswith('[') else [float(x) for x in radlist.split(',')]
radii = [float(x) for x in lst]   # Ensure all elements are floats

# Derive output directory by trimming known file suffix (9 characters from end)
out_dir = data_path[:-9]
print(out_dir)  # Print output directory location
print(radii)    # Print parsed radii list for confirmation

# --------------------- User Parameters ---------------------
uniform_size = float(os.getenv("GMESH_UNIFORM_SIZE", "0.4")) # Global mesh element size scaling factor
model_name = "n_shells_sphere"   # Model name for Gmsh
SHOW_POPUP = False               # Disable GUI popup for automated runs

# ---------------------- Basic Checks -----------------------
# Validate radii input to ensure mesh generation consistency
if len(radii) < 1:
    raise ValueError("Need at least one radius to form a sphere.")
if any(r <= 0.0 for r in radii):
    raise ValueError("All radii must be positive.")
if len(radii) > 1 and any(r2 <= r1 for r1, r2 in zip(radii, radii[1:])):
    raise ValueError("Radii must be strictly increasing.")

N = len(radii) - 1        # Number of shell layers between spheres
r_min, r_max = radii[0], radii[-1]  # Minimum and maximum sphere radii
eps = 1e-12               # Small tolerance for floating-point stability
min_tk = min((radii[i] - radii[i-1]) for i in range(1, len(radii))) if len(radii) > 1 else float("inf")

# Initialize Gmsh and the geometric model
gmsh.initialize()
gmsh.model.add(model_name)

# ---------------- Threads and Terminal Options ----------------
num_threads = int(os.getenv("GMESH_THREADS", "1"))  # Number of parallel threads
gmsh.option.setNumber("General.NumThreads", max(1, num_threads))              # Global threading
gmsh.option.setNumber("Mesh.MaxNumThreads2D", max(1, int(os.getenv("GMESH_MAX_THREADS_2D", str(num_threads)))))  # 2D threads
gmsh.option.setNumber("Mesh.MaxNumThreads3D", max(1, int(os.getenv("GMESH_MAX_THREADS_3D", str(num_threads)))))  # 3D threads
gmsh.option.setNumber("General.Terminal", 1)  # Enable terminal output for real-time feedback

# ---------------- OCC Healing and Sewing Options ----------------
# Enable corrections for small edges/faces and geometry stitching to ensure watertight meshes
gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)
gmsh.option.setNumber("Geometry.OCCSewFaces", 1)
gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
gmsh.option.setNumber("Geometry.Tolerance", float(os.getenv("GEOMETRY_TOL", "1e-8")))

# ---------------- Circle Discretization and Randomization ----------------
# Control circle facet density and slight randomness for avoiding degenerate meshes
gmsh.option.setNumber("Mesh.MinimumCirclePoints", int(os.getenv("GMESH_MIN_CIRCLE_POINTS", "8")))
gmsh.option.setNumber("Mesh.RandomFactor", float(os.getenv("GMESH_RANDOM_FACTOR", "5e-5")))

# ------------------- Create Concentric Spheres -------------------
# Add multiple concentric spheres with specified radii
sphere_tags = [gmsh.model.occ.addSphere(0.0, 0.0, 0.0, r) for r in radii]
gmsh.model.occ.synchronize()  # Apply geometric changes to the Gmsh model

# ------------------- Boolean Cuts to Form Shells -----------------
# Create distinct shell volumes by subtracting inner spheres from the outer ones
shell_tags = [None] * (N + 1)
for k in range(N, 0, -1):  
    out_dimtags, _ = gmsh.model.occ.cut([(3, sphere_tags[k])], [(3, sphere_tags[k - 1])], removeTool=False)
    gmsh.model.occ.synchronize()
    new_vols = [(d, t) for (d, t) in out_dimtags if d == 3]  # Extract 3D volume entities
    if not new_vols:
        gmsh.finalize()
        raise RuntimeError(f"Cut failed for k={k} (r[{k-1}] -> r[{k}]).")
    shell_tags[k] = new_vols[0][1]

# Remove duplicate entities left over from boolean operations
try:
    gmsh.model.occ.removeAllDuplicates()
except Exception:
    pass
gmsh.model.occ.synchronize()

# ----------------------- Physical Groups -------------------------
# Define physical volume groups for each shell for post-processing and boundary conditions
if N == 0:
    gmsh.model.addPhysicalGroup(3, [sphere_tags[0]], tag=1)
    gmsh.model.setPhysicalName(3, 1, "Inner")
else:
    gmsh.model.addPhysicalGroup(3, [sphere_tags[0]], tag=1)
    gmsh.model.setPhysicalName(3, 1, "Inner")
    for k in range(1, N + 1):
        phys_tag = k + 1
        gmsh.model.addPhysicalGroup(3, [shell_tags[k]], tag=phys_tag)
        gmsh.model.setPhysicalName(3, phys_tag, f"Shell{k}")

# ---------------- Mesh Algorithm and Global Options ----------------
# Specify mesh generation algorithms and global smoothing/optimization settings
gmsh.option.setNumber("Mesh.Algorithm", 6)        # Use Frontal-Delaunay for 2D meshing
gmsh.option.setNumber("Mesh.Algorithm3D", 10)     # Use HXT for fast parallel 3D meshing
gmsh.option.setNumber("Mesh.ElementOrder", 1)     # First-order elements for simpler meshes
gmsh.option.setNumber("Mesh.Smoothing", 0)        # Disable post-smoothing to retain control
gmsh.option.setNumber("Mesh.Optimize", 0)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)

# Disable curvature, point, and boundary-based local mesh refinements to rely solely on background field
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# Output mesh format configuration
gmsh.option.setNumber("Mesh.Binary", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)

# --------------------- Field builder ------------------
# This function defines mesh size fields for adaptive refinement around sphere interfaces.
# These fields control the target element size at each point in space.
# Smaller elements are placed near thin shells or sharp gradients, improving accuracy.

last_field_ids = []  # Keeps track of existing fields to clear them when rebuilding

def build_fields(n_thick_loc, n_circ_near_loc, n_circ_far_loc, band_coeff, band_extra_mult):
    global last_field_ids
    fld = gmsh.model.mesh.field  # Shortcut for the field API
    
    # Remove old mesh fields before creating new ones
    for fid in last_field_ids:
        try:
            fld.remove(fid)
        except Exception:
            pass
    last_field_ids = []

    r_expr = "sqrt(x*x + y*y + z*z)"  # Expression for current radius in spherical coordinates

    # Compute absolute mesh size limits based on shell thickness
    tk_min = min_tk if math.isfinite(min_tk) else uniform_size
    h_min_abs = max(uniform_size, min(0.35 * tk_min, 0.03 * r_max))  # Fine size near thin shells
    h_max_abs = (2.0 * math.pi * r_max) / max(n_circ_far_loc, 6)     # Coarser size far away

    # Smoothly increase element size radially (quadratic scaling)
    alpha = 2.0
    slope = 0.0 if (r_max - r_min) < eps else (h_max_abs - h_min_abs) / ((r_max - r_min) ** alpha)

    # Field for gradual radial growth
    g_field = fld.add("MathEval"); last_field_ids.append(g_field)
    fld.setString(g_field, "F", f"{h_min_abs} + {slope}*(({r_expr} - {r_min})*({r_expr} - {r_min}))")

    # Define limits on minimum and maximum element sizes
    cmin_field = fld.add("MathEval"); last_field_ids.append(cmin_field)
    fld.setString(cmin_field, "F", f"{h_min_abs}")
    cmax_field = fld.add("MathEval"); last_field_ids.append(cmax_field)
    fld.setString(cmax_field, "F", f"{h_max_abs}")
    g_clamped_max = fld.add("Max"); last_field_ids.append(g_clamped_max)
    fld.setNumbers(g_clamped_max, "FieldsList", [g_field, cmin_field])
    g_final = fld.add("Min"); last_field_ids.append(g_final)
    fld.setNumbers(g_final, "FieldsList", [g_clamped_max, cmax_field])
    to_min = [g_final]  # Store all refinement fields for later combination

    # Helper function: returns target mesh size around a given radius rr
    def h_circ_near_fun(rr, ncirc=None):
        nn = n_circ_near_loc if ncirc is None else ncirc
        return (2.0 * math.pi * max(rr, eps)) / max(nn, 6)

    # Minimum allowable element size to avoid degenerate elements
    h_floor_small = max(0.6 * uniform_size, min(0.20 * tk_min, 0.01 * r_max), 5e-4)

    # ---------------- Core band refinement ----------------
    # Finer mesh near the innermost sphere
    r0 = r_min
    h_core_target = min(max(tk_min, r0) / max(n_thick_loc, 1), h_circ_near_fun(max(r0, tk_min)))
    h_core = max(h_floor_small, min(h_core_target, h_max_abs))

    d_core = fld.add("MathEval"); last_field_ids.append(d_core)
    fld.setString(d_core, "F", f"abs({r_expr} - {r0})")

    # Threshold field constrains mesh size transition outward from the core
    t_core = fld.add("Threshold"); last_field_ids.append(t_core)
    try:
        fld.setNumber(t_core, "IField", d_core)
    except Exception:
        fld.setNumber(t_core, "InField", d_core)
    fld.setNumber(t_core, "LcMin", h_core)
    fld.setNumber(t_core, "LcMax", h_max_abs)
    fld.setNumber(t_core, "DistMin", 0.0)
    fld.setNumber(t_core, "DistMax", r0 + eps)
    fld.setNumber(t_core, "StopAtDistMax", 1)
    to_min.append(t_core)

    # ---------------- Interface band refinement ----------------
    # For every pair of spheres, generate mesh refinement near their boundaries
    ULTRA_THIN_RATIO = float(os.getenv("GMESH_ULTRA_THIN_RATIO", "0.01"))
    ABS_THIN = float(os.getenv("GMESH_ABS_THIN", "0.6"))

    for k in range(1, len(radii)):
        rin, rout = radii[k - 1], radii[k]
        tk = max(rout - rin, eps)
        rmid = 0.5 * (rin + rout)
        is_ultra_thin = (tk < ABS_THIN) or (tk / max(rout, eps) < ULTRA_THIN_RATIO)

        if is_ultra_thin:
            # Use finer mesh for ultra-thin shells
            hk_target = min(tk, h_circ_near_fun(rmid, ncirc=8))
            hk = max(h_floor_small, min(hk_target, h_max_abs))
            band = max(0.50 * tk, 1.50 * hk)
        else:
            # Allow coarser meshing for thicker shells
            is_thick = (tk > max(1.0, 0.25 * rin))
            hk_target = min(tk / (3.0 if is_thick else max(n_thick_loc, 1)), h_circ_near_fun(rmid))
            hk = max(h_floor_small, min(hk_target, h_max_abs))
            band = max(band_coeff * tk, band_extra_mult * hk)

        # Define threshold fields for inner and outer shell interfaces
        din = fld.add("MathEval"); last_field_ids.append(din)
        fld.setString(din, "F", f"abs({r_expr} - {rin})")
        dout = fld.add("MathEval"); last_field_ids.append(dout)
        fld.setString(dout, "F", f"abs({r_expr} - {rout})")

        tin = fld.add("Threshold"); last_field_ids.append(tin)
        try:
            fld.setNumber(tin, "IField", din)
        except Exception:
            fld.setNumber(tin, "InField", din)
        fld.setNumber(tin, "LcMin", hk)
        fld.setNumber(tin, "LcMax", h_max_abs)
        fld.setNumber(tin, "DistMin", 0.0)
        fld.setNumber(tin, "DistMax", band)
        fld.setNumber(tin, "StopAtDistMax", 1)

        tout = fld.add("Threshold"); last_field_ids.append(tout)
        try:
            fld.setNumber(tout, "IField", dout)
        except Exception:
            fld.setNumber(tout, "InField", dout)
        fld.setNumber(tout, "LcMin", hk)
        fld.setNumber(tout, "LcMax", h_max_abs)
        fld.setNumber(tout, "DistMin", 0.0)
        fld.setNumber(tout, "DistMax", band)
        fld.setNumber(tout, "StopAtDistMax", 1)

        # Combine inner and outer interface fields
        tpair = fld.add("Min"); last_field_ids.append(tpair)
        fld.setNumbers(tpair, "FieldsList", [tin, tout])
        to_min.append(tpair)

    # Combine all refinement regions into one background field
    f_bg = fld.add("Min"); last_field_ids.append(f_bg)
    fld.setNumbers(f_bg, "FieldsList", to_min)
    fld.setAsBackgroundMesh(f_bg)

    # Enforce overall mesh size bounds
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h_floor_small)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h_max_abs)
    return f_bg  # Return the ID of the background field

# Create mesh field configuration with slightly coarser settings by default
build_fields(
    n_thick_loc=int(os.getenv("GMESH_N_THICK", "3")),
    n_circ_near_loc=int(os.getenv("GMESH_N_CIRC_NEAR", "12")),
    n_circ_far_loc=int(os.getenv("GMESH_N_CIRC_FAR", "6")),
    band_coeff=0.55,
    band_extra_mult=1.4,
)

# --------------------- Preflight scaling ------------------
# Determines mesh scaling factor to control total node count and maintain performance.

BUDGET_2D_NODES = int(os.getenv("GMESH_BUDGET_2D_NODES", "12000"))
SCALE_INIT = float(os.getenv("GMESH_SCALE_INIT", "1.0"))
SCALE_MAX = float(os.getenv("GMESH_SCALE_MAX", "80.0"))
MAX_ITERS = int(os.getenv("GMESH_SCALE_ITERS", "6"))
PANIC_SCALE_TRIG = float(os.getenv("GMESH_PANIC_TRIG", "3.0"))
SKIP_PREFLIGHT = bool(int(os.getenv("GMESH_SKIP_PREFLIGHT", "0")))
PREFLIGHT_TIME_BUDGET = float(os.getenv("GMESH_PREFLIGHT_SEC", "8"))

# Automatically tune the scale to stay within a node budget
def preflight_scale_by_2d_budget():
    scale = SCALE_INIT
    t0 = time.monotonic()
    n2d = 0
    for _ in range(MAX_ITERS):
        if time.monotonic() - t0 > PREFLIGHT_TIME_BUDGET:
            return min(scale * 2.0, SCALE_MAX), n2d
        gmsh.option.setNumber("Mesh.MeshSizeFactor", scale)
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(1)
        gmsh.model.mesh.generate(2)
        nodeTags, _, _ = gmsh.model.mesh.getNodes()
        n2d = len(nodeTags)
        if n2d <= BUDGET_2D_NODES or scale >= SCALE_MAX:
            return scale, n2d
        over = n2d / max(BUDGET_2D_NODES, 1)
        scale *= max(1.10, (over ** 0.6))
        scale = min(scale, SCALE_MAX)
    return scale, n2d

# Apply preflight logic or skip if instructed
thin_ratio = (min_tk / max(r_max, eps)) if math.isfinite(min_tk) else 1.0
if SKIP_PREFLIGHT:
    scale_used, n2d_final = (max(SCALE_INIT, 3.0 if thin_ratio < 0.01 else 1.8), 0)
else:
    scale_used, n2d_final = preflight_scale_by_2d_budget()
    if scale_used >= PANIC_SCALE_TRIG:
        # Retry with looser mesh constraints if mesh is too dense
        build_fields(
            n_thick_loc=max(2, int(os.getenv("GMESH_N_THICK", "3")) - 1),
            n_circ_near_loc=max(10, int(os.getenv("GMESH_N_CIRC_NEAR", "12")) // 2),
            n_circ_far_loc=int(os.getenv("GMESH_N_CIRC_FAR", "6")),
            band_coeff=0.45,
            band_extra_mult=1.25,
        )
        scale_used, n2d_final = preflight_scale_by_2d_budget()

# ----------------------- Mesh generation -------------------------
# Generate final mesh at the determined element size scale

gmsh.model.mesh.clear()
gmsh.option.setNumber("Mesh.MeshSizeFactor", scale_used)
try:
    gmsh.model.mesh.generate(1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.generate(3)
except Exception:
    # If fast mesher fails, fallback to classical Delaunay for stability
    gmsh.model.mesh.clear()
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.generate(3)

# -------------------------- Save mesh ----------------------------
# Save generated mesh and optionally open in GUI

os.makedirs(out_dir, exist_ok=True)
outfile = os.path.join(out_dir, f"{model_name}_{N+1}_shells.msh")
gmsh.write(outfile)
print(f"Mesh written to: {outfile}")

if SHOW_POPUP:
    gmsh.fltk.run()

gmsh.finalize()  # Cleanly close Gmsh session
