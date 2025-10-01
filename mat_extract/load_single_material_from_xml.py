from pathlib import Path
from typing import Optional, Union
import openmc


def load_single_material_from_xml(xml_path: Union[str, Path],
                                  *,
                                  name: Optional[str] = None,
                                  id: Optional[int] = None,
                                  index: Optional[int] = None) -> openmc.Material:
    """Load exactly one openmc.Material from a materials XML file.

    The caller must specify exactly one of (name, id, index).

    Parameters
    ----------
    xml_path : str or pathlib.Path
        Path to an OpenMC materials XML (e.g., 'materials.xml').
    name : str, optional
        Material name to select (exact match).
    id : int, optional
        Material ID to select.
    index : int, optional
        Zero-based index in the file order.

    Returns
    -------
    openmc.Material

    Raises
    ------
    ValueError
        If zero or more than one selector is given, or if no match is found.
    """
    xml_path = Path(xml_path)
    if sum(x is not None for x in (name, id, index)) != 1:
        raise ValueError("Specify exactly one of: name=..., id=..., index=...")

    # Parse XML -> Materials collection
    mats = openmc.Materials.from_xml(str(xml_path))

    # Selection logic
    if name is not None:
        matches = [m for m in mats if (m.name == name)]
        key_desc = f'name="{name}"'
    elif id is not None:
        matches = [m for m in mats if (getattr(m, "id", None) == id)]
        key_desc = f'id={id}'
    else:  # index is not None
        if index < 0 or index >= len(mats):
            raise ValueError(f"index={index} is out of range (0..{len(mats)-1})")
        matches = [mats[index]]
        key_desc = f'index={index}'

    if not matches:
        # Build a helpful message listing options
        listing = "\n".join([f"- idx {i:2d}: id={getattr(m, 'id', None)}, name={repr(m.name)}"
                             for i, m in enumerate(mats)])
        raise ValueError(
            f"No material found for {key_desc} in {xml_path}.\n"
            f"Available materials:\n{listing}"
        )

    # Return the unique selection (do not clone unless you need to mutate)
    return matches[0]
