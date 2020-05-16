from __future__ import absolute_import, print_function
from sage.all import factorial
from moment_polytopes import *
import pytest


def test_hrepr_irred():
    # unit square as an irredundant set of inequalities
    unit_square = HRepr(ieqs=[((1, 0), 0), ((0, 1), 0), ((-1, 0), -1), ((0, -1), -1),])
    assert unit_square == unit_square.irred()
    assert (0, 0) in unit_square
    assert (-1, 0) not in unit_square

    # add a redundant inequality
    unit_square_redund = unit_square & HRepr(ieqs=[((1, 1), 0)])
    assert unit_square != unit_square_redund
    assert unit_square == unit_square_redund.irred()


@pytest.mark.xfail
def test_irred_linearities():
    line = HRepr(ieqs=[((1, -1), 0), ((-1, 1), 0),])

    # irred does *NOT* convert the two inequalities into a linearity :-(
    line_irred = line.irred()
    assert len(line_irred.eqns) == 1
    assert len(line_irred.ieqs) == 0

    # line_irred = line.vrepr().hrepr().irred() <-- this works...


def test_hrepr_vrepr():
    # triangle
    triangle_hrepr = HRepr(ieqs=[((1, 0), 0), ((0, 1), 0), ((-1, -1), -1),])
    triangle_vrepr = VRepr(vertices=[(0, 0), (1, 0), (0, 1),])
    assert triangle_hrepr.vrepr() == triangle_vrepr
    assert triangle_vrepr.hrepr() == triangle_hrepr


def test_hrepr_sage():
    triangle_hrepr = HRepr(ieqs=[((1, 0), 0), ((0, 1), 0), ((-1, -1), -1),])

    assert HRepr.from_sage(triangle_hrepr.to_sage()) == triangle_hrepr


def test_vrepr_sage():
    triangle_vrepr = VRepr(vertices=[(0, 0), (1, 0), (0, 1),])

    assert VRepr.from_sage(triangle_vrepr.to_sage()) == triangle_vrepr


def test_hrepr_vertices():
    # fail if there are extremal rays
    quadrant = HRepr(ieqs=[((1, 0), 0), ((0, 1), 0),])
    with pytest.raises(ValueError):
        quadrant.vertices()


def test_overflow_expectations():
    # check my expectations with regards to numpy and long integers
    import numpy as np

    N = factorial(171)
    # assert isinstance(N, long)

    with pytest.raises(OverflowError):
        np.array([N], dtype=np.int)

    assert np.array([N]).dtype == object


def test_ambient_dim():
    # usually, ambient dimensions are detected automatically
    half_space = HRepr(ieqs=[((1, 0), 0)])
    assert half_space.ambient_dim == 2

    origin = VRepr(vertices=[(0, 0, 0)])
    assert origin.ambient_dim == 3

    # if no data is given, ambient dimension cannot be detected
    with pytest.raises(ValueError):
        HRepr()
    with pytest.raises(ValueError):
        VRepr()

    # it works if we specify the ambient dimension manually
    assert HRepr(ambient_dim=2).ambient_dim == 2
    assert VRepr(ambient_dim=2).ambient_dim == 2


def test_default_hrepr():
    assert not HRepr(ambient_dim=1).to_sage().is_empty()


def test_default_vrepr():
    assert VRepr(ambient_dim=1).to_sage().is_empty()
