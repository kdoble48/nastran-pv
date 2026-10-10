"""Fixed code tables for categorical cell arrays.

``ElementType``, ``RigidType`` and ``ConnectorType`` codes are global (not per block), so one
indexed lookup table colours every block consistently. Unknown card types are
appended at run time.
"""
from __future__ import annotations

ELEMENT_TYPES = [
    "CTRIA3", "CTRIAR", "CTRIA6", "CQUAD4", "CQUADR", "CQUAD8", "CQUAD", "CSHEAR",
    "CHEXA", "CPENTA", "CTETRA", "CPYRAM",
    "CBAR", "CBEAM", "CROD", "CONROD", "CTUBE", "CVISC", "CBEND",
]
RIGID_TYPES = ["RBE2", "RBE3", "RBAR", "RBAR1", "RROD", "RBE1", "RTRPLT", "RTRPLT1", "RSPLINE", "RSSCON"]
CONNECTOR_TYPES = ["MPC", "CELAS1", "CELAS2", "CDAMP1", "CDAMP2", "CGAP", "CFAST"]
BUSH_AXES = ["ellipsoid", "K1", "K2", "K3"]  # Axis cell array of the bush glyphs


def code_of(table: list[str], name: str) -> int:
    try:
        return table.index(name)
    except ValueError:
        table.append(name)
        return len(table) - 1
