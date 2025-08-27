#!/usr/bin/env python3
"""
Scan a directory tree for OpenMC geometries that are exclusively spherical.

What it does
------------
- Walks the filesystem starting at ROOT (set in __main__).
- Only inspects directories whose name equals 'openmc' case-insensitively.
- In each such directory, looks for 'geometry.xml'.
- Parses all <surface> elements (namespace-agnostic).
- Confirms every surface has type="sphere".
- Extracts sphere definitions from either:
    (a) coeffs="x0 y0 z0 r"   (preferred if present), or
    (b) attributes r[/radius] and optional x0,y0,z0 (defaults to 0).
- If exclusively spherical, prints the folder path and the sorted shell radii.

No third-party deps.
"""

import os
import sys
import math
import pickle
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple, Optional

# ----------------- XML helpers -----------------


def _is_surface_tag(elem: ET.Element) -> bool:
    """True if element's tag ends with 'surface' (namespace-agnostic)."""
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
    # Prefer coeffs
    coeffs = s.attrib.get("coeffs")
    if coeffs is not None:
        parsed = _parse_coeffs(coeffs)
        if parsed is None:
            return None
        x0, y0, z0, r = parsed
        if not (math.isfinite(r) and r > 0.0):
            return None
        return {"r": r, "x0": x0, "y0": y0, "z0": z0}
    # Fallback to r / radius and optional x0,y0,z0
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

    Returns
    -------
    (is_exclusively_spherical, spheres)
      - True if all <surface> elements have type='sphere' and are valid.
      - spheres: list of {'r','x0','y0','z0'} for each spherical surface.
    """
    try:
        tree = ET.parse(path)
        root = tree.getroot()
    except Exception:
        return (False, [])
    surfaces = [e for e in root.iter() if _is_surface_tag(e)]
    if not surfaces:
        # No surfaces found → cannot assert “exclusively spherical”
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
    Walk `root_dir`, find directories named 'openmc' (case-insensitive),
    parse geometry.xml, and collect spherical-only results.

    Returns list of dicts:
      {
        'openmc_dir': <path>,
        'geometry_xml': <full path>,
        'exclusively_spherical': <bool>,
        'radii': [float, ...]    # only if exclusively_spherical
      }
    """
    results: List[dict] = []
    root_dir = os.path.abspath(root_dir)

    for dirpath, dirnames, filenames in os.walk(root_dir):
        if os.path.basename(dirpath).lower() != "openmc":
            continue
        geom_path = os.path.join(dirpath, "geometry.xml")
        if not os.path.isfile(geom_path):
            continue
        is_spherical, spheres = analyze_geometry_xml(geom_path)
        rec = {
            "openmc_dir": os.path.relpath(dirpath, root_dir),
            "geometry_xml": geom_path,
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


# ----------------- Main (set ROOT here) -----------------

if __name__ == "__main__":
    # Set the root directory you want to scan:
    ROOT = r"../icsbep_original/"  # <-- change this to your starting folder

    if not os.path.isdir(ROOT):
        print(f"Error: '{ROOT}' is not a directory.", file=sys.stderr)
        sys.exit(2)
    results = scan_for_openmc_spheres(ROOT)
    print_report(results)

    out_pkl = "spherical_results.pkl"
    with open(out_pkl, "wb") as f:
        pickle.dump(results, f)
    print(f"Saved pickle results to {out_pkl}")
