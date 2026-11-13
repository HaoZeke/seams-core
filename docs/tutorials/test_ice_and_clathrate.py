"""Ice and clathrate labels through pydseamslib.

Install the package the paper uses::

    pip install pydseamslib

``import pydseams`` is the module. ``import pydseamslib`` is the same
package. ``seams chill-plus`` and ``Frame.chill_plus`` both use the
cutoff neighbour list. ``--graph mutual`` is the other list. On the
cubic mW frame both lists label every site cubic. The 2.11 wheel of
pydseamslib still prints 12 interfacial clathrate on the mixed frame.
"""

import pytest

pydseams = pytest.importorskip("pydseams")


def test_mixed_tip4p_cutoff_chill_plus(tip4p):
    frame = pydseams.read(tip4p, atom_type=2)
    ice = frame.chill_plus()
    assert ice.interClathrate == 9
    assert ice.water == 241
    assert ice.clathrate == 0
    cages = frame.cages()
    assert cages.n_ih == 0
    assert cages.n_ic == 0
    assert cages.n_water == 250


def test_cubic_mw_is_all_cubic(cubic):
    frame = pydseams.read(cubic, frame=1, atom_type=1)
    assert frame.n_atoms == 4096
    assert frame.chill_plus().cubic == 4096
    assert frame.cages().n_ic == 4096
    assert frame.cages().n_water == 0


def test_sI_512_cages(si):
    frame = pydseams.read(si, atom_type=1)
    found = frame.cages_by_signature("512")
    atoms = set()
    for cage in found:
        atoms.update(cage["vertices"])
    assert len(found) == 2
    assert len(atoms) == 40
