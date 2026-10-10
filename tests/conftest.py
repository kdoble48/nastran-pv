from pathlib import Path

import numpy as np
import pytest
from vtkmodules.util.numpy_support import vtk_to_numpy

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def run_deck():
    return load_deck(DATA / "tiny_run.bdf", echo=False)


@pytest.fixture(scope="session")
def bulk_deck():
    return load_deck(DATA / "tiny_bulk.bdf", echo=False)


@pytest.fixture(scope="session")
def run_mb():
    return build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))


@pytest.fixture(scope="session")
def topo_base():
    return build_multiblock(load_deck(DATA / "topology.bdf", echo=False))


@pytest.fixture()
def topo_mb(topo_base):
    """Fresh attribute containers per test: topology operations add arrays in place."""
    from nastran_pv.fields.common import shallow_clone

    return shallow_clone(topo_base)


def arr(dataset, name, where="point"):
    data = getattr(dataset, {"point": "GetPointData", "cell": "GetCellData", "field": "GetFieldData"}[where])()
    a = data.GetArray(name)
    assert a is not None, f"missing {where} array {name}"
    return vtk_to_numpy(a)


def strings(dataset, name):
    a = dataset.GetFieldData().GetAbstractArray(name)
    return [a.GetValue(i) for i in range(a.GetNumberOfValues())]


def row_of(dataset, key_array, key, where="point"):
    keys = arr(dataset, key_array, where)
    hits = np.nonzero(keys == key)[0]
    assert len(hits) == 1, f"{key_array}={key} found {len(hits)} times"
    return hits[0]


__all__ = ["arr", "strings", "row_of", "get_block", "DATA"]
