import os
import re
import math
import xml.etree.ElementTree as ET

# ---- Filepaths ----
ROOT = r"../../icsbep_original"
OUT  = r"../../spherical_cases"
# -------------------


# Return a sorted list of floats with near-duplicates removed within a tolerance.
def unique_sorted(vals, tol=1e-9):
    vals = sorted(float(v) for v in vals)
    out = []
    for v in vals:
        if not out or abs(v - out[-1]) > tol:
            out.append(v)
    return out


# Convert an openmc-relative directory path into "<prefix-before-openmc>/case-<n>".
def derive_dest_relpath(openmc_dir_rel):
    parts = [p for p in os.path.normpath(openmc_dir_rel).split(os.sep) if p and p != "."]

    # Find the last "openmc" segment so we handle nested paths safely.
    j = -1
    for i in range(len(parts) - 1, -1, -1):
        if parts[i].lower() == "openmc":
            j = i
            break
    if j < 0:
        return None

    prefix = parts[:j]
    after = parts[j + 1:]

    # Pull the first digits out of a "case-like" segment (case-2, c2, case_2, etc.).
    case_num = None
    if after:
        seg = after[0]
        if re.match(r"^(case|c)[-_]?\d+", seg, flags=re.IGNORECASE):
            m = re.search(r"(\d+)", seg)
            if m:
                try:
                    case_num = int(m.group(1))
                except Exception:
                    case_num = None
    if case_num is None:
        case_num = 1

    if prefix:
        return os.path.join(*prefix, f"case-{case_num}")
    return f"case-{case_num}"


# Parse a geometry.xml and return (True, radii) only if every surface is a sphere with a valid radius.
def analyze_geometry_xml(geom_path):
    try:
        root = ET.parse(geom_path).getroot()
    except Exception:
        return (False, [])

    radii = []
    any_surface = False

    for elem in root.iter():
        # Only look at <surface> tags (including namespaced ones).
        if elem.tag.split("}")[-1] != "surface":
            continue

        any_surface = True
        stype = (elem.attrib.get("type") or "").strip().lower()
        if stype != "sphere":
            return (False, [])

        # Prefer coeffs="x0 y0 z0 r" if present; otherwise use r= or radius=.
        r = None
        coeffs = elem.attrib.get("coeffs")
        if coeffs:
            parts = coeffs.replace(",", " ").split()
            if len(parts) >= 4:
                try:
                    r = float(parts[3])
                except Exception:
                    r = None
        else:
            for key in ("r", "radius"):
                if key in elem.attrib:
                    try:
                        r = float(elem.attrib[key])
                        break
                    except Exception:
                        r = None

        if r is None or (not math.isfinite(r)) or r <= 0.0:
            return (False, [])

        radii.append(r)

    if not any_surface:
        return (False, [])

    return (True, radii)


def main():
    if not os.path.isdir(ROOT):
        raise FileNotFoundError(f"ROOT not found: {ROOT}")
    os.makedirs(OUT, exist_ok=True)

    seen_geom = 0
    spherical = 0
    copied_geom = 0
    wrote_radii = 0
    copied_mat = 0
    missing_mat = 0

    root_abs = os.path.abspath(ROOT)

    for dirpath, _, filenames in os.walk(root_abs):
        # Only consider folders at/below an "openmc" segment.
        dir_parts = [p for p in os.path.normpath(dirpath).split(os.sep) if p]
        if not any(p.lower() == "openmc" for p in dir_parts):
            continue

        if "geometry.xml" not in filenames:
            continue

        seen_geom += 1
        geom_path = os.path.join(dirpath, "geometry.xml")

        ok, radii = analyze_geometry_xml(geom_path)
        if not ok:
            continue

        radii = unique_sorted(radii)
        if not radii:
            continue

        openmc_dir_rel = os.path.relpath(dirpath, root_abs)
        dest_rel = derive_dest_relpath(openmc_dir_rel)
        if not dest_rel:
            continue

        dest_case_dir = os.path.join(OUT, dest_rel)

        # ---- mesh outputs ----
        mesh_dir = os.path.join(dest_case_dir, "mesh")
        os.makedirs(mesh_dir, exist_ok=True)

        # Copy geometry.xml -> mesh/geometry.xml (simple byte copy).
        dest_geom = os.path.join(mesh_dir, "geometry.xml")
        with open(geom_path, "rb") as fsrc, open(dest_geom, "wb") as fdst:
            fdst.write(fsrc.read())
        copied_geom += 1

        # Write mesh/radii.txt
        radii_path = os.path.join(mesh_dir, "radii.txt")
        with open(radii_path, "w", encoding="utf-8") as f:
            for r in radii:
                f.write(f"{r:.17g}\n")
        wrote_radii += 1

        # ---- materials copy (if present next to geometry.xml) ----
        src_mat = os.path.join(dirpath, "materials.xml")
        if os.path.isfile(src_mat):
            mat_dir = os.path.join(dest_case_dir, "materials")
            os.makedirs(mat_dir, exist_ok=True)

            dest_mat = os.path.join(mat_dir, "materials.xml")
            with open(src_mat, "rb") as fsrc, open(dest_mat, "wb") as fdst:
                fdst.write(fsrc.read())
            copied_mat += 1
        else:
            missing_mat += 1

        spherical += 1

    print("ROOT:", os.path.abspath(ROOT))
    print("OUT :", os.path.abspath(OUT))
    print("geometry.xml found under openmc:", seen_geom)
    print("spherical cases copied        :", spherical)
    print("geometry.xml copied to mesh/  :", copied_geom)
    print("radii.txt written             :", wrote_radii)
    print("materials.xml copied          :", copied_mat)
    print("spherical cases missing materials.xml:", missing_mat)


if __name__ == "__main__":
    main()
