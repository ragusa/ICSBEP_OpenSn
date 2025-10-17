#!/usr/bin/env python3
"""
Scan a directory tree for OpenMC geometries that are exclusively spherical,
then create a mirrored mesh tree with normalized case folders.

What it does
------------
- Walks the filesystem starting at ROOT.
- Processes any directory at or below an 'openmc' segment (case-insensitive).
- In each such directory, looks for 'geometry.xml', parses <surface> elements,
  and confirms every surface is type="sphere".
- If exclusively spherical, records radii and writes them into a mirrored tree
  under MESH_ROOT at: <prefix-before-openmc>/case-<n>,
  where <n> is taken from the first case-like segment after 'openmc' (case-2, c2 → case-2),
  or defaults to case-1 if missing/non-standard.

No third-party deps.
"""

import os
import sys
import math
import pickle
import re
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple, Optional


# ----------------- XML helpers -----------------

def _is_surface_tag(elem: ET.Element) -> bool:
    return elem.tag.split("}")[-1] == "surface"


def _to_float_or_none(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        return None


def _parse_coeffs(coeffs: str) -> Optional[Tuple[float, float, float, float]]:
    """
    Parse coeffs="x0 y0 z0 r" into floats.
    Returns (x0, y0, z0, r) or None if invalid.
    """
    try:
        parts = coeffs.replace(",", " ").split()
        if len(parts) < 4:
            return None
        x0, y0, z0, r = map(float, parts[:4])
        return (x0, y0, z0, r)
    except Exception:
        return None


def _sphere_from_surface_elem(s: ET.Element) -> Optional[Dict[str, float]]:
    """
    Build a sphere dict from a <surface type='sphere'> element.
    Prefer 'coeffs'. Fallback to attributes r/radius and x0,y0,z0.
    Returns {'r','x0','y0','z0'} or None if invalid/missing data.
    """
    coeffs = s.attrib.get("coeffs")
    if coeffs is not None:
        parsed = _parse_coeffs(coeffs)
        if parsed is None:
            return None
        x0, y0, z0, r = parsed
        if not (math.isfinite(r) and r > 0.0):
            return None
        return {"r": r, "x0": x0, "y0": y0, "z0": z0}

    r = _to_float_or_none(s.attrib.get("r"))
    if r is None:
        r = _to_float_or_none(s.attrib.get("radius"))
    if r is None or not (math.isfinite(r) and r > 0.0):
        return None
    x0 = _to_float_or_none(s.attrib.get("x0"))
    y0 = _to_float_or_none(s.attrib.get("y0"))
    z0 = _to_float_or_none(s.attrib.get("z0"))
    x0 = 0.0 if x0 is None else x0
    y0 = 0.0 if y0 is None else y0
    z0 = 0.0 if z0 is None else z0

    return {"r": float(r), "x0": float(x0), "y0": float(y0), "z0": float(z0)}


def analyze_geometry_xml(path: str) -> Tuple[bool, List[Dict[str, float]]]:
    """
    Analyze geometry.xml at `path`.
    Returns (is_exclusively_spherical, spheres).
    """
    try:
        tree = ET.parse(path)
        root = tree.getroot()
    except Exception:
        return (False, [])
    surfaces = [e for e in root.iter() if _is_surface_tag(e)]
    if not surfaces:
        return (False, [])
    spheres: List[Dict[str, float]] = []
    for s in surfaces:
        stype = (s.attrib.get("type") or "").strip().lower()
        if stype != "sphere":
            return (False, [])
        sp = _sphere_from_surface_elem(s)
        if sp is None:
            return (False, [])
        spheres.append(sp)
    return (True, spheres)


# ----------------- Geometry summarizers -----------------

def unique_sorted_radii(
    spheres: List[Dict[str, float]], tol: float = 1e-9
) -> List[float]:
    """Return unique radii sorted ascending, deduplicated within tolerance."""
    vals = sorted(s["r"] for s in spheres)
    uniq: List[float] = []
    for r in vals:
        if not uniq or abs(r - uniq[-1]) > tol:
            uniq.append(r)
    return uniq


# ----------------- Scanner -----------------

def scan_for_openmc_spheres(root_dir: str) -> List[dict]:
    """
    Walk `root_dir`, find any directories that are at or below an 'openmc' segment,
    parse geometry.xml there, and collect spherical-only results.
    """
    results: List[dict] = []
    root_dir = os.path.abspath(root_dir)

    for dirpath, dirnames, filenames in os.walk(root_dir):
        parts = os.path.normpath(dirpath).split(os.sep)
        if not any(seg.lower() == "openmc" for seg in parts):
            continue

        geom_path = os.path.join(dirpath, "geometry.xml")
        if not os.path.isfile(geom_path):
            continue

        is_spherical, spheres = analyze_geometry_xml(geom_path)
        rec = {
            "openmc_dir": os.path.relpath(dirpath, root_dir),  # may be .../openmc or .../openmc/case-*
            "geometry_xml": geom_path,                         # absolute path
            "exclusively_spherical": bool(is_spherical),
        }
        if is_spherical:
            rec["radii"] = unique_sorted_radii(spheres)
        results.append(rec)

    return results


def print_report(results: List[dict]) -> None:
    """Print only the folders whose geometry is exclusively spherical, with radii."""
    print("\nFolders with exclusively spherical geometries:\n" + "-" * 60)
    any_printed = False
    for rec in results:
        if not rec["exclusively_spherical"]:
            continue
        any_printed = True
        print(rec["openmc_dir"])
        rads = rec.get("radii", [])
        print(
            "  Shell radii (sorted): "
            + (", ".join(f"{r:g}" for r in rads) if rads else "(none)")
        )
        print()
    if not any_printed:
        print("(none)\n")


# ----------------- Mesh tree builder (normalized case-#) -----------------

# case-like leader and digits
CASE_LEAD_RE = re.compile(r'^(case|c)[-_]?\d+', re.IGNORECASE)
DIGITS_RE = re.compile(r'(\d+)')


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_radii_file(path: str, radii, fname: str = "radii.txt") -> None:
    # ensure the "mesh" subfolder exists and write radii.txt inside it
    target_dir = os.path.join(path, "mesh")
    ensure_dir(target_dir)
    vals = sorted(set(float(r) for r in radii))
    out_path = os.path.join(target_dir, fname)
    with open(out_path, "w", encoding="utf-8") as f:
        for r in vals:
            f.write(f"{r:.17g}\n")



def _split_any(path: str) -> List[str]:
    # Platform-agnostic splitter: normalize slashes to '/'
    path = (path or "").replace("\\", "/")
    return [p for p in path.split("/") if p and p != "."]


def _strip_filename(parts: List[str]) -> List[str]:
    # Drop trailing filename (e.g., geometry.xml) if present
    return parts[:-1] if parts and "." in parts[-1] else parts


def _find_openmc_index(parts: List[str]) -> int:
    for i in range(len(parts) - 1, -1, -1):
        if parts[i].lower() == "openmc":
            return i
    return -1


def _extract_case_number(seg: str) -> Optional[int]:
    if not CASE_LEAD_RE.match(seg or ""):
        return None
    m = DIGITS_RE.search(seg)
    if not m:
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None


def _derive_target_rel_from_record(rec: dict) -> Optional[str]:
    """
    Build target relative path under mesh_root as:
        <prefix-before-openmc> / case-<n>
    where <n> is extracted from the first case-like segment immediately after 'openmc',
    or defaults to 1 if no case segment exists or is non-standard.
    """
    path_hint = (
        rec.get("openmc_dir")
        or rec.get("openmc_dir_rel")
        or rec.get("geometry_xml")
        or rec.get("geometry_xml_rel")
        or ""
    )
    parts = _strip_filename(_split_any(path_hint))
    j = _find_openmc_index(parts)
    if j == -1:
        return None

    prefix = parts[:j]
    after = parts[j + 1 :]

    case_num = None
    if after:
        case_num = _extract_case_number(after[0])
    if case_num is None:
        case_num = 1

    case_dir = f"case-{case_num}"
    rel = os.path.join(*prefix, case_dir) if prefix else case_dir
    return os.path.normpath(rel)


def create_mesh_tree_from_records(
    records: List[dict], mesh_root: str, radii_filename: str = "radii.txt"
) -> None:
    mesh_root = os.path.abspath(mesh_root)
    ensure_dir(mesh_root)
    created = 0
    for rec in records:
        if not isinstance(rec, dict):
            continue
        if not rec.get("exclusively_spherical", False):
            continue
        target_rel = _derive_target_rel_from_record(rec)
        if not target_rel:
            continue
        target_dir = os.path.join(mesh_root, target_rel)
        radii = rec.get("radii") or []
        if not radii:
            continue
        write_radii_file(target_dir, radii, fname=radii_filename)
        created += 1
    print(f"Created {created} folder(s) with {radii_filename} under: {mesh_root}")


# ----------------- Main -----------------

if __name__ == "__main__":
    # Set paths here
    ROOT = r"../../icsbep_original/"   # starting folder to scan
    MESH_ROOT = r"../../spherical_cases"          # output mirrored root

    if not os.path.isdir(ROOT):
        print(f"Error: '{ROOT}' is not a directory.", file=sys.stderr)
        sys.exit(2)

    # 1) Scan and report
    results = scan_for_openmc_spheres(ROOT)
    print_report(results)

    # 2) Save pickle for archival/debug
    out_pkl = "spherical_results.pkl"
    with open(out_pkl, "wb") as f:
        pickle.dump(results, f)
    print(f"Saved pickle results to {out_pkl}")

    # 3) Create mesh tree with normalized case folders
    create_mesh_tree_from_records(results, MESH_ROOT, radii_filename="radii.txt")
