#!/usr/bin/env python3
import os
import sys
import re
import pickle
from typing import Dict, Any, Iterable, List, Optional


CASE_SUBDIR_RE = re.compile(r"^(case-\d+|c-\d+)$", re.IGNORECASE)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_radii_file(path: str, radii: Iterable[float], fname: str = "radii.txt") -> None:
    ensure_dir(path)
    # numeric only, one per line, sorted & deduped
    vals = sorted({float(r) for r in radii})
    with open(os.path.join(path, fname), "w", encoding="utf-8") as f:
        for r in vals:
            f.write(f"{r:.17g}\n")


def load_results_from_pickle(pkl_path: str) -> List[Dict[str, Any]]:
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    # Normalize to a list of record dicts
    if isinstance(data, list):
        records = data
    elif isinstance(data, dict):
        if "results" in data and isinstance(data["results"], list):
            records = data["results"]
        else:
            # Assume dict-of-records -> take values
            records = list(data.values())
    else:
        raise TypeError("Pickled object must be a list or dict containing result records.")
    return records


def _normalize(path: str) -> str:
    """Normalize slashes and collapse .., . etc."""
    return os.path.normpath(path.replace("\\", "/"))


def _pick_case_subdir(parts_after_openmc: List[str]) -> Optional[str]:
    """
    From the list of path parts after 'openmc', return the first subfolder that matches
    'case-N' or 'c-N' (N integer), case-insensitive. If none match, return None.
    """
    for p in parts_after_openmc:
        if CASE_SUBDIR_RE.match(p):
            return p
    return None


def _derive_rel_path_from_geometry_xml(geom_path: str) -> Optional[str]:
    """
    Given a geometry.xml path like:
        ".../<case>/openmc/case-1/geometry.xml"  or  ".../<case>/openmc/geometry.xml"
    return:
        "<case>/openmc/case-1"                   or  "<case>/openmc"
    Only keeps a single child subfolder under 'openmc' if it matches (case-insensitive):
        'case-N' or 'c-N'.
    """
    if not geom_path:
        return None
    parts = geom_path.replace("\\", "/").split("/")
    if not parts:
        return None

    # Remove the filename (e.g., geometry.xml) if present
    if parts[-1].lower().endswith(".xml"):
        parts = parts[:-1]

    # Find the last occurrence of 'openmc' in the path
    openmc_idx = None
    for i, p in enumerate(parts):
        if p.lower() == "openmc":
            openmc_idx = i
    if openmc_idx is None:
        return None

    # The directory right before 'openmc' is the "<case>" folder (if present)
    start_idx = max(0, openmc_idx - 1)

    # Look for a matching subfolder right under 'openmc'
    after = parts[openmc_idx + 1 :]  # parts after 'openmc'
    chosen = _pick_case_subdir(after)

    rel_parts = parts[start_idx : openmc_idx + 1]  # "<case>/openmc"
    if chosen:
        rel_parts.append(chosen)  # keep only the matching 'case-N' or 'c-N' layer

    rel_path = "/".join(rel_parts)
    return rel_path or None


def _sanitize_openmc_dir_path(odir: str) -> str:
    """
    For openmc_dir[_rel] fallbacks, keep:
        "<case>/openmc" or "<case>/openmc/<case-N|c-N>" if present.
    If odir already points to deeper paths, strip any extra layers that don't match the rule.
    """
    parts = odir.replace("\\", "/").strip("/").split("/")
    # find 'openmc'
    try:
        openmc_idx = max(i for i, p in enumerate(parts) if p.lower() == "openmc")
    except ValueError:
        # no 'openmc' at all; return original normalized
        return odir.replace("\\", "/").rstrip("/")

    start_idx = max(0, openmc_idx - 1)  # include "<case>"
    after = parts[openmc_idx + 1 :]
    chosen = _pick_case_subdir(after)

    rel_parts = parts[start_idx : openmc_idx + 1]
    if chosen:
        rel_parts.append(chosen)
    return "/".join(rel_parts)


def _derive_rel_openmc_dir(rec: Dict[str, Any]) -> Optional[str]:
    """
    Best-effort derivation of a relative path of the form:
      "<case>/openmc[/case-N|c-N]"
    Preference order:
      1) geometry_xml_rel
      2) geometry_xml
      3) openmc_dir_
