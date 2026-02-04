from pathlib import Path
import shutil
import pickle
import xml.etree.ElementTree as ET  # stdlib XML parser

def build_material_nuclide_map(xml_path: Path) -> dict:
    """
    Returns a dict mapping material 'id' -> list of (name, ao) tuples.
    Only collects <nuclide> children that define an 'ao' attribute.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    material_map = {}
    for mat in root.findall("./material"):
        mat_id = mat.get("id")
        if not mat_id:
            continue

        pairs = []
        for nuc in mat.findall("./nuclide"):
            name = nuc.get("name")
            ao = nuc.get("ao")  # collect only atom fraction
            if name is None or ao is None:
                continue
            try:
                ao_val = float(ao)
            except ValueError:
                ao_val = ao  # keep raw string if it can't be parsed as float
            pairs.append((name, ao_val))

        material_map[mat_id] = pairs

    return material_map

# Assume this script lives in mat_extract; repo root is one level up via "../"
script_dir = Path(__file__).resolve().parent
repo_root = script_dir / "../.."

spherical_cases_dir = repo_root / "spherical_cases"
icsbep_original_dir = repo_root / "icsbep_original"

if not spherical_cases_dir.is_dir():
    raise FileNotFoundError(f"Missing folder: {spherical_cases_dir}")
if not icsbep_original_dir.is_dir():
    raise FileNotFoundError(f"Missing folder: {icsbep_original_dir}")

# Collect immediate child folders in spherical_cases as the destination case names
case_names = sorted([p.name for p in spherical_cases_dir.iterdir() if p.is_dir()])

print(f"Found {len(case_names)} case folders in spherical_cases")

copied = 0
skipped_missing_source = []
skipped_missing_materials = []

for name in case_names:
    # Source layout expected:
    # icsbep_original/<name>/openmc/(case-#?)/materials.xml
    src_case_dir = icsbep_original_dir / name
    if not src_case_dir.is_dir():
        skipped_missing_source.append(name)
        continue

    openmc_dir = src_case_dir / "openmc"
    if not openmc_dir.is_dir():
        skipped_missing_source.append(name)
        continue

    # Prefer explicit case-* subfolders; if none, treat openmc itself as "case-1"
    case_dirs = sorted([d for d in openmc_dir.iterdir() if d.is_dir() and d.name.lower().startswith("case-")])
    if not case_dirs:
        case_dirs = [openmc_dir]  # will map to "case-1" on destination

    any_copied_for_name = False

    for cd in case_dirs:
        # Resolve source materials.xml path
        src_materials = cd / "materials.xml"

        # Determine destination case label
        if cd == openmc_dir:
            dest_case_label = "case-1"
        else:
            dest_case_label = cd.name

        if not src_materials.is_file():
            skipped_missing_materials.append(f"{name}/{dest_case_label}")
            continue

        # Destination: spherical_cases/<name>/<case-#>/materials/materials.xml
        dest_case_dir = spherical_cases_dir / name / dest_case_label / "materials"
        dest_case_dir.mkdir(parents=True, exist_ok=True)

        dest_xml = dest_case_dir / "materials.xml"
        shutil.copy2(src_materials, dest_xml)
        copied += 1
        any_copied_for_name = True
        print(f"Copied: {name}/openmc/{dest_case_label}/materials.xml -> spherical_cases/{name}/{dest_case_label}/materials/materials.xml")

        # Build and save the pickle alongside the copied XML
        try:
            nuclide_map = build_material_nuclide_map(dest_xml)
            dest_pkl = dest_case_dir / "materials.pkl"
            with open(dest_pkl, "wb") as f:
                pickle.dump(nuclide_map, f, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Wrote pickle: spherical_cases/{name}/{dest_case_label}/materials/materials.pkl")
        except Exception as e:
            print(f"Warning: failed to parse/pickle {dest_xml}: {e}")

    # If there were explicit case-* folders but none had materials.xml, track at top level too
    if not any_copied_for_name and case_dirs and case_dirs != [openmc_dir]:
        skipped_missing_materials.append(f"{name}/(all cases)")

print(f"\nCompleted. Files copied: {copied}")
if skipped_missing_source:
    print("No matching folder in icsbep_original for:", ", ".join(skipped_missing_source))
if skipped_missing_materials:
    print("Missing materials.xml under openmc for:", ", ".join(skipped_missing_materials))
