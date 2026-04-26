from ast import literal_eval
import gmsh
import sys
import os
import math


def parse_radii(text):
    text = text.strip()
    if text.startswith("["):
        values = literal_eval(text)
    else:
        values = [float(x) for x in text.split(",") if x.strip()]
    radii = [float(x) for x in values]

    if len(radii) < 1:
        raise ValueError("Need at least one radius to form a sphere.")
    if any(r <= 0.0 for r in radii):
        raise ValueError("All radii must be positive.")
    if any(r2 <= r1 for r1, r2 in zip(radii, radii[1:])):
        raise ValueError("Radii must be strictly increasing.")

    return radii


def safe_out_dir(data_path):
    if len(data_path) > 9:
        return data_path[:-9]
    parent = os.path.dirname(data_path)
    return parent if parent else "."


def set_options():
    num_threads = int(os.getenv("GMESH_THREADS", "1"))

    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.NumThreads", max(1, num_threads))
    gmsh.option.setNumber(
        "Mesh.MaxNumThreads2D",
        max(1, int(os.getenv("GMESH_MAX_THREADS_2D", str(num_threads))))
    )
    gmsh.option.setNumber(
        "Mesh.MaxNumThreads3D",
        max(1, int(os.getenv("GMESH_MAX_THREADS_3D", str(num_threads))))
    )

    gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
    gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)
    gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
    gmsh.option.setNumber("Geometry.OCCSewFaces", 1)
    gmsh.option.setNumber("Geometry.Tolerance", float(os.getenv("GEOMETRY_TOL", "1e-8")))

    gmsh.option.setNumber("Mesh.MinimumCirclePoints", int(os.getenv("GMESH_MIN_CIRCLE_POINTS", "12")))
    gmsh.option.setNumber("Mesh.RandomFactor", float(os.getenv("GMESH_RANDOM_FACTOR", "1e-6")))

    gmsh.option.setNumber("Mesh.Algorithm", 6)      # Frontal-Delaunay for surfaces
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)    # Delaunay first
    gmsh.option.setNumber("Mesh.ElementOrder", 1)

    gmsh.option.setNumber("Mesh.Smoothing", int(os.getenv("GMESH_SMOOTHING", "6")))
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)

    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)


def build_geometry_and_groups(radii):
    n_shells = len(radii) - 1
    sphere_tags = [gmsh.model.occ.addSphere(0.0, 0.0, 0.0, r) for r in radii]
    gmsh.model.occ.synchronize()

    inner_volume = sphere_tags[0]
    shell_volumes = []

    for k in range(1, len(radii)):
        copied_outer = gmsh.model.occ.copy([(3, sphere_tags[k])])
        gmsh.model.occ.synchronize()

        out_dimtags, _ = gmsh.model.occ.cut(
            copied_outer,
            [(3, sphere_tags[k - 1])],
            removeObject=True,
            removeTool=False,
        )
        gmsh.model.occ.synchronize()

        new_vols = [tag for dim, tag in out_dimtags if dim == 3]
        if not new_vols:
            raise RuntimeError(f"Cut failed for shell {k}: r={radii[k - 1]} -> {radii[k]}")

        shell_volumes.append(new_vols[0])

    try:
        gmsh.model.occ.removeAllDuplicates()
    except Exception:
        pass
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(3, [inner_volume], tag=1)
    gmsh.model.setPhysicalName(3, 1, "Inner")

    for k, vol in enumerate(shell_volumes, start=1):
        phys_tag = k + 1
        gmsh.model.addPhysicalGroup(3, [vol], tag=phys_tag)
        gmsh.model.setPhysicalName(3, phys_tag, f"Shell{k}")

    return n_shells


def build_smooth_background_field(radii):
    r_min = radii[0]
    r_max = radii[-1]
    eps = 1e-12

    if len(radii) > 1:
        min_tk = min(r2 - r1 for r1, r2 in zip(radii, radii[1:]))
    else:
        min_tk = float("inf")

    base_lc_env = os.getenv("GMESH_BASE_LC", "").strip()
    if base_lc_env:
        base_lc = float(base_lc_env)
    else:
        if math.isfinite(min_tk):
            base_lc = min(0.85 * min_tk, 0.18 * r_max)
            base_lc = max(base_lc, 0.04 * r_max)
        else:
            base_lc = 0.12 * r_max

    inner_scale = float(os.getenv("GMESH_INNER_SCALE", "0.90"))
    outer_scale = float(os.getenv("GMESH_OUTER_SCALE", "1.12"))

    lc_inner = max(1e-4, inner_scale * base_lc)
    lc_outer = max(lc_inner, outer_scale * base_lc)

    span = max(r_max - r_min, eps)
    r_expr = "sqrt(x*x + y*y + z*z)"

    expr = f"{lc_inner} + ({lc_outer - lc_inner}) * (({r_expr} - {r_min}) / {span})"

    field = gmsh.model.mesh.field
    bg = field.add("MathEval")
    field.setString(bg, "F", expr)
    field.setAsBackgroundMesh(bg)

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc_inner)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_outer)

    return base_lc, lc_inner, lc_outer


def clear_mesh():
    try:
        gmsh.model.mesh.clear()
    except Exception:
        pass


def generate_mesh_and_count(scale, optimize_after=True):
    clear_mesh()
    gmsh.option.setNumber("Mesh.MeshSizeFactor", float(scale))

    last_err = None
    for alg3d in (1, 4, 10):  # Delaunay, Frontal, HXT fallback
        try:
            gmsh.option.setNumber("Mesh.Algorithm3D", alg3d)
            gmsh.model.mesh.generate(3)

            if optimize_after:
                try:
                    gmsh.model.mesh.optimize("")
                except Exception:
                    pass
                try:
                    gmsh.model.mesh.optimize("Netgen")
                except Exception:
                    pass

            node_tags, _, _ = gmsh.model.mesh.getNodes()
            elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)

            n_nodes = len(node_tags)
            n_elems_3d = sum(len(tags) for tags in elem_tags)

            return n_nodes, n_elems_3d, alg3d
        except Exception as err:
            last_err = err
            clear_mesh()

    raise last_err


def tune_scale_to_target(target_nodes, node_tol_frac):
    enable_tune = bool(int(os.getenv("GMESH_ENABLE_3D_TUNE", "1")))
    scale_init = float(os.getenv("GMESH_SCALE_INIT", "1.0"))
    scale_min = float(os.getenv("GMESH_SCALE_MIN", "0.20"))
    scale_max = float(os.getenv("GMESH_SCALE_MAX", "8.0"))
    max_iters = int(os.getenv("GMESH_SCALE_ITERS", "10"))

    target_lo = max(50, int((1.0 - node_tol_frac) * target_nodes))
    target_hi = max(target_lo + 1, int((1.0 + node_tol_frac) * target_nodes))

    scale = min(max(scale_init, scale_min), scale_max)
    trials = []

    n_nodes, n_elems_3d, alg3d = generate_mesh_and_count(scale)
    trials.append((scale, n_nodes, n_elems_3d, alg3d))
    print(f"[tune] scale={scale:.6g} -> nodes={n_nodes}, elems3D={n_elems_3d}, alg3D={alg3d}")

    if (not enable_tune) or (target_lo <= n_nodes <= target_hi):
        return min(trials, key=lambda t: abs(t[1] - target_nodes))

    too_many = None
    too_few = None

    if n_nodes > target_hi:
        too_many = (scale, n_nodes, n_elems_3d, alg3d)

        for _ in range(max_iters):
            step = min(2.0, max(1.20, (n_nodes / max(target_nodes, 1)) ** (1.0 / 3.0)))
            new_scale = min(scale * step, scale_max)
            if abs(new_scale - scale) < 1e-12:
                break

            scale = new_scale
            n_nodes, n_elems_3d, alg3d = generate_mesh_and_count(scale)
            trials.append((scale, n_nodes, n_elems_3d, alg3d))
            print(f"[tune] scale={scale:.6g} -> nodes={n_nodes}, elems3D={n_elems_3d}, alg3D={alg3d}")

            if target_lo <= n_nodes <= target_hi:
                return min(trials, key=lambda t: abs(t[1] - target_nodes))

            if n_nodes > target_hi:
                too_many = (scale, n_nodes, n_elems_3d, alg3d)
            else:
                too_few = (scale, n_nodes, n_elems_3d, alg3d)
                break
    else:
        too_few = (scale, n_nodes, n_elems_3d, alg3d)

        for _ in range(max_iters):
            step = min(2.0, max(1.20, (target_nodes / max(n_nodes, 1)) ** (1.0 / 3.0)))
            new_scale = max(scale / step, scale_min)
            if abs(new_scale - scale) < 1e-12:
                break

            scale = new_scale
            n_nodes, n_elems_3d, alg3d = generate_mesh_and_count(scale)
            trials.append((scale, n_nodes, n_elems_3d, alg3d))
            print(f"[tune] scale={scale:.6g} -> nodes={n_nodes}, elems3D={n_elems_3d}, alg3D={alg3d}")

            if target_lo <= n_nodes <= target_hi:
                return min(trials, key=lambda t: abs(t[1] - target_nodes))

            if n_nodes < target_lo:
                too_few = (scale, n_nodes, n_elems_3d, alg3d)
            else:
                too_many = (scale, n_nodes, n_elems_3d, alg3d)
                break

    if too_many is not None and too_few is not None:
        s_many = too_many[0]
        s_few = too_few[0]

        for _ in range(max_iters):
            mid = math.sqrt(s_many * s_few)
            n_nodes, n_elems_3d, alg3d = generate_mesh_and_count(mid)
            trials.append((mid, n_nodes, n_elems_3d, alg3d))
            print(f"[tune] scale={mid:.6g} -> nodes={n_nodes}, elems3D={n_elems_3d}, alg3D={alg3d}")

            if target_lo <= n_nodes <= target_hi:
                return min(trials, key=lambda t: abs(t[1] - target_nodes))

            if n_nodes > target_hi:
                s_many = mid
            else:
                s_few = mid

    return min(trials, key=lambda t: abs(t[1] - target_nodes))


def main():
    radlist = sys.argv[1]
    data_path = sys.argv[2]

    radii = parse_radii(radlist)
    out_dir = safe_out_dir(data_path)

    model_name = "n_shells_sphere"
    show_popup = bool(int(os.getenv("GMESH_SHOW_POPUP", "0")))
    target_nodes = int(os.getenv("GMESH_TARGET_NODES", "1000"))
    node_tol_frac = float(os.getenv("GMESH_NODE_TOL_FRAC", "0.08"))

    print(out_dir)
    print(radii)

    gmsh.initialize()
    try:
        gmsh.model.add(model_name)
        set_options()

        n_shells = build_geometry_and_groups(radii)
        base_lc, lc_inner, lc_outer = build_smooth_background_field(radii)

        print(f"[size] base_lc={base_lc:.6g}, lc_inner={lc_inner:.6g}, lc_outer={lc_outer:.6g}")

        best_scale, best_nodes, best_elems_3d, best_alg3d = tune_scale_to_target(
            target_nodes=target_nodes,
            node_tol_frac=node_tol_frac,
        )

        print(f"[final-choice] target_nodes={target_nodes}")
        print(
            f"[final-choice] scale={best_scale:.6g}, "
            f"nodes={best_nodes}, elems3D={best_elems_3d}, alg3D={best_alg3d}"
        )

        final_nodes, final_elems_3d, final_alg3d = generate_mesh_and_count(best_scale)
        print(
            f"[final-mesh] scale={best_scale:.6g}, "
            f"nodes={final_nodes}, elems3D={final_elems_3d}, alg3D={final_alg3d}"
        )

        os.makedirs(out_dir, exist_ok=True)
        outfile = os.path.join(out_dir, f"{model_name}_{n_shells + 1}_shells.msh")
        gmsh.write(outfile)
        print(f"Mesh written to: {outfile}")

        if show_popup:
            gmsh.fltk.run()

    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()