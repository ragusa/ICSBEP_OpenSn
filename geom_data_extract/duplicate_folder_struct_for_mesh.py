#!/usr/bin/env python3
import os
import sys
import pickle


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_radii_file(path: str, radii, fname: str = "radii.txt") -> None:
    ensure_dir(path)
    # numeric only, one per line, sorted & deduped
    vals = sorted(set(float(r) for r in radii))
    with open(os.path.join(path, fname), "w", encoding="utf-8") as f:
        for r in vals:
            f.write(f"{r:.17g}\n")


def load_results_from_pickle(pkl_path: str):
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
        raise TypeError(
            "Pickled object must be a list or dict containing result records."
        )
    return records


def create_mesh_tree_from_pickle(
    pkl_path: str, mesh_root: str, radii_filename: str = "radii.txt"
) -> None:
    records = load_results_from_pickle(pkl_path)
    mesh_root = os.path.abspath(mesh_root)
    ensure_dir(mesh_root)

    created = 0
    for rec in records:
        if not isinstance(rec, dict):
            continue
        # Require spherical-only entries
        if not rec.get("exclusively_spherical", False):
            continue
        # Prefer the relative path field captured during scanning
        rel_path = rec.get("openmc_dir") or rec.get("openmc_dir_rel")

        # If missing, attempt to infer from geometry_xml (…/<case>/openmc/geometry.xml)
        if not rel_path:
            geom = rec.get("geometry_xml") or rec.get("geometry_xml_rel") or ""
            if geom:
                parts = geom.replace("\\", "/").split("/")
                # try to capture "<case>/openmc"
                # geometry.xml path ends with ".../<case>/openmc/geometry.xml"
                if len(parts) >= 3:
                    rel_path = "/".join(parts[-3:-1])  # "<case>/openmc"
            if not rel_path:
                continue
        rel_path_norm = os.path.normpath(rel_path)
        target_dir = os.path.join(mesh_root, rel_path_norm)

        radii = rec.get("radii") or []
        if not radii:
            continue
        write_radii_file(target_dir, radii, fname=radii_filename)
        created += 1
    print(f"Created {created} folder(s) with {radii_filename} under: {mesh_root}")


if __name__ == "__main__":
    # === Set these two values ===
    PKL_PATH = r"./spherical_results.pkl"  # your saved results (pickle)
    MESH_ROOT = r"../mesh"  # new mirrored root name
    # ============================

    if not os.path.isfile(PKL_PATH):
        print(f"Error: pickle not found: {PKL_PATH}", file=sys.stderr)
        sys.exit(2)
    create_mesh_tree_from_pickle(PKL_PATH, MESH_ROOT, radii_filename="radii.txt")
