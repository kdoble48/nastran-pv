"""docs/api.md: every name it documents imports, has a docstring and is in its module's __all__."""
import importlib
import inspect
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
API = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")
DOCUMENTED = re.findall(r"^\.\. auto(?:function|class|data):: (nastran_pv\.[\w.]+)", API, re.M)
# modules whose __all__ is exactly what the page shows (topology, post, fields and results_h5 list
# a few internal helpers too: the page is the contract, docs/versioning.md#public-api)
EXACT = ("nastran_pv.deck", "nastran_pv.assemble", "nastran_pv.results", "nastran_pv.vtkio", "nastran_pv.vtkhdf",
         "nastran_pv.result_types", "nastran_pv.fields.derived", "nastran_pv.fields.tables",
         "nastran_pv.fields.calc", "nastran_pv.sections.props")


def _split(dotted):
    module, _, name = dotted.rpartition(".")
    return module, name


def test_documented_names_are_public():
    assert len(DOCUMENTED) > 40
    for dotted in DOCUMENTED:
        module, name = _split(dotted)
        mod = importlib.import_module(module)
        obj = getattr(mod, name)
        assert name in getattr(mod, "__all__", ()), f"{dotted} is on api.md but not in {module}.__all__"
        if not isinstance(obj, tuple):  # REGISTRY's doc is a #: comment
            assert inspect.getdoc(obj), f"{dotted} has no docstring"


def test_exact_modules_list_only_the_page():
    shown = {}
    for dotted in DOCUMENTED:
        module, name = _split(dotted)
        shown.setdefault(module, set()).add(name)
    for module in EXACT:
        assert set(importlib.import_module(module).__all__) == shown[module], module
