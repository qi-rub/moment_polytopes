from __future__ import absolute_import, print_function
import pytest
from sage.all import vector, QQ, ZZ, matrix, RootSystem
from moment_polytopes import *


def test_dual_root_primitive():
    root_system = 'A1xA2'

    # evaluated on the simple roots, this is the vector (2, 1, 1)
    H = vector([1, -1, 1, 0, -1])
    assert is_dual_root_primitive(root_system, H)

    # evaluated on the simple roots, this is the vector (6, 0, 15)
    H = vector([3, -3, 5, 5, -10])
    H_prim = vector([1, -1, QQ('5/3'), QQ('5/3'), QQ('-10/3')])
    assert not is_dual_root_primitive(root_system, H)
    assert dual_root_primitive(root_system, H) == H_prim


def test_GL4_positive_weyl_chamber():
    # positive Weyl chamber
    pwc_hrepr_irred = positive_weyl_chamber_hrepr('A3').irred()
    pwc_hrepr_expected = HRepr(ieqs=[
        ((1, -1, 0, 0), 0),
        ((0, 1, -1, 0), 0),
        ((0, 0, 1, -1), 0),
    ])
    assert pwc_hrepr_irred == pwc_hrepr_expected


def test_GL3_fundamental():
    # too many parts
    with pytest.raises(AssertionError):
        R = weyl_module(3, [4, 3, 2, 1])

    # fundamental representation
    R = weyl_module(3, [1])

    # weights
    assert R.weights == map(vector, [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ])

    # affine hull of weights
    assert R.dimension_affine_hull_weights == 2
    assert R.reduced_positive_weyl_chamber_hrepr.eqns == [
        (vector([1, 1, 1]), 1),
    ]

    # negative root action
    V = vector
    assert R.negative_roots == map(vector, [
        [-1, 1, 0],
        [-1, 0, 1],
        [0, -1, 1],
    ])
    assert R.negative_root_action(0) == matrix([
        [0, 0, 0],
        [1, 0, 0],
        [0, 0, 0],
    ])
    assert R.negative_root_action(1) == matrix([
        [0, 0, 0],
        [0, 0, 0],
        [QQ('1/2'), 0, 0],
    ])
    assert R.negative_root_action(2) == matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, QQ('1/2'), 0],
    ])


def test_GL12_weyl_module_21():
    # representation with highest weight [2,1]
    R = weyl_module(12, [2, 1])
    V = vector

    # affine hull of weights
    assert R.dimension_affine_hull_weights == 11
    assert R.reduced_positive_weyl_chamber_hrepr.eqns == [
        (V([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]), 3),
    ]

    # check eqns. (3.21) in my thesis -- NB: we are using zero-based indexing here!
    alpha = R.negative_roots.index(V([0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]))
    v = R.tableau_vector([[1, 1], [2]])
    w = R.tableau_vector([[1, 1], [8]]) / QQ(720)
    assert R.negative_root_action(alpha) * v == w

    alpha = R.negative_roots.index(V([-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    v = R.tableau_vector([[1, 1], [2]])
    w = R.tableau_vector([[1, 2], [2]])
    assert R.negative_root_action(alpha) * v == w


def test_GL6_third_antisymmetric():
    # representation with highest weight [1,1,1]
    R = weyl_module(6, [1, 1, 1])
    V = vector

    # check dimension
    assert R.dimension == 20

    # affine hull of weights
    assert R.dimension_affine_hull_weights == 5
    assert R.reduced_positive_weyl_chamber_hrepr.eqns == [
        (V([1, 1, 1, 1, 1, 1]), 3),
    ]

    # negative root action
    alpha = R.negative_roots.index(V([-1, 0, 0, 1, 0, 0]))
    v = R.tableau_vector([[1], [2], [3]])
    w = R.tableau_vector([[2], [3], [4]])
    assert R.negative_root_action(alpha) * v == w

    alpha = R.negative_roots.index(V([0, -1, 0, 1, 0, 0]))
    v = R.tableau_vector([[1], [2], [3]])
    w = R.tableau_vector([[1], [3], [4]])
    assert R.negative_root_action(alpha) * v == -w  # note the negative sign


def test_GL4_two_det():
    # representation with highest weight [-2,-2,-2,-2]
    R = weyl_module(4, [-2, -2, -2, -2])
    V = vector

    # check dimension
    assert R.dimension == 1

    # check weights
    assert R.weights == [
        V([-2, -2, -2, -2]),
    ]

    # affine hull of weights
    assert R.dimension_affine_hull_weights == 0
    assert R.reduced_positive_weyl_chamber_hrepr.eqns == [
        (V([1, 1, 1, 1]), -8),
    ]

    # negative root action
    for idx, alpha in enumerate(R.negative_roots):
        assert R.negative_root_action(idx) == matrix([[0]])


def test_GL5_adjoint():
    R = weyl_module(5, [1, 0, 0, 0, -1])
    V = vector

    # check highest weight and partition
    assert R.highest_weight == V([1, 0, 0, 0, -1])

    # check dimension
    assert R.dimension == (5 - 1) * (5 + 1)

    # check weights
    got = map(tuple, R.weights)
    expected = [(0, 0, 0, 0, 0)] * (5 - 1) + map(
        tuple, R.positive_roots + R.negative_roots)
    assert sorted(got) == sorted(expected)

    # affine hull of weights
    assert R.dimension_affine_hull_weights == 4
    assert R.reduced_positive_weyl_chamber_hrepr.eqns == [
        (V([1, 1, 1, 1, 1]), 0),
    ]

    # negative root action
    alpha = R.negative_roots.index(V([-1, 0, 0, 1, 0]))  # alpha_{4,1}
    v = R.tableau_vector([[1, 1], [2], [4], [5]])  # |1><3|
    w = R.tableau_vector([[1, 4], [2], [4], [5]]) / QQ(8)  # |4><3|
    assert R.negative_root_action(alpha) * v == w

    v = R.tableau_vector([[1, 1], [2], [3], [5]])  # |1><4|
    w = R.tableau_vector(
        [[1, 2], [3], [4], [5]]) * QQ('1/2') + R.tableau_vector(
            [[1, 3], [2], [4], [5]]) * QQ('-1/6') + R.tableau_vector(
                [[1, 4], [2], [3], [5]]) * QQ('1/6')  # |4><4| - |1><1|
    assert R.negative_root_action(alpha) * v == w


def test_external_tensor_product():
    R = external_tensor_product([2, 3])
    V = vector

    # root system
    assert R.root_system == RootSystem('A1xA2')

    # negative roots
    assert R.negative_roots == [
        V([-1, 1, 0, 0, 0]),
        V([0, 0, -1, 1, 0]),
        V([0, 0, -1, 0, 1]),
        V([0, 0, 0, -1, 1]),
    ]

    # weights
    assert R.weights == [
        V([1, 0, 1, 0, 0]),
        V([1, 0, 0, 1, 0]),
        V([1, 0, 0, 0, 1]),
        V([0, 1, 1, 0, 0]),
        V([0, 1, 0, 1, 0]),
        V([0, 1, 0, 0, 1]),
    ]

    # affine hull of weights
    assert R.dimension_affine_hull_weights == 3
    assert R.reduced_positive_weyl_chamber_hrepr.eqns == [
        (V([0, 0, 1, 1, 1]), 1),
        (V([1, 1, 0, 0, 0]), 1),
    ]

    # negative root action]
    alpha = R.negative_roots.index(V([0, 0, -1, 0, 1]))
    assert R.negative_root_action(alpha) == matrix([
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [QQ('1/2'), 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, QQ('1/2'), 0, 0],
    ])
