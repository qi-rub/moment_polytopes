from __future__ import absolute_import, print_function
import pytest
from moment_polytopes import *


def test_two_three_six(algorithm):
    R = external_tensor_product([2, 3, 6])
    T = ressayre_tester(R, algorithm=algorithm)

    # one of the many inequalites for 2 x 3 x 6 (cf. Wernli/Klyachko)
    assert T.is_ressayre(((-1, 0, -1, 1, 0, 3, 2, 1, 2, 1, 0), 1))


def test_three_three_nine(algorithm):
    R = external_tensor_product([3, 3, 9])
    T = ressayre_tester(R, algorithm=algorithm)

    # the last inequality of Klyachko
    assert T.is_ressayre(((3, 0, -3, -5, 4, 1, 8, 5, 2, 2, -1, -4, -7, -4, -1), 0))

    # another inequality by Klyachko (the one that M.V. checked)
    assert T.is_ressayre(((-1, 0, 1, 0, -1, 1, 2, 1, 0, 1, 0, 0, -1, -1, -2), 0))


def test_fermi_four_eight(algorithm):
    R = weyl_module(8, [1, 1, 1, 1])
    T = ressayre_tester(R, algorithm=algorithm)

    assert T.is_ressayre(((-1, 0, 0, 0, 0, 0, 0, 0), -1))
    assert T.is_ressayre(((0, 0, 0, 0, -1, 1, 1, 1), 0))
    assert T.is_ressayre(((-1, 1, 0, 0, 0, 0, 1, 1), 0))
    assert T.is_ressayre(((-1, 0, 1, 0, 0, 1, 0, 1), 0))
    assert T.is_ressayre(((-1, 0, 0, 1, 0, 1, 1, 0), 0))
    assert T.is_ressayre(((-1, 0, 0, 1, 1, 0, 0, 1), 0))
    assert T.is_ressayre(((0, 0, -1, 1, 0, 0, 1, 1), 0))
    assert T.is_ressayre(((0, -1, 0, 1, 0, 1, 0, 1), 0))
    assert T.is_ressayre(((0, -1, -1, 0, -1, 0, 0, 1), -2))
    assert T.is_ressayre(((-1, 0, -1, 0, 0, -1, 0, 1), -2))
    assert T.is_ressayre(((-1, -1, 0, 0, 0, 0, -1, 1), -2))
    assert T.is_ressayre(((-1, -1, -1, 1, 0, 0, 0, 0), -2))
    assert T.is_ressayre(((-1, 0, 0, -1, -1, 0, 0, 1), -2))
    assert T.is_ressayre(((-1, -1, 0, 0, -1, 1, 0, 0), -2))
    assert T.is_ressayre(((-1, 0, -1, 0, -1, 0, 1, 0), -2))


def test_fermi_three_eight(algorithm):
    R = weyl_module(8, [1, 1, 1])
    T = ressayre_tester(R, algorithm=algorithm)

    # two of out many
    assert T.is_ressayre(((-1, -2, 3, 1, 2, 1, 0, -1), 0))
    assert T.is_ressayre(((1, 0, -1, -2, 3, 2, 1, -1), 0))


@pytest.mark.parametrize("d", [6, 7])
def test_spin_orbit(algorithm, d):
    R = external_tensor_product([weyl_module(d, [2, 1]), 2])
    T = ressayre_tester(R, algorithm=algorithm)

    # special case of Eqn. (3.19) in my thesis
    assert T.is_ressayre(((-2, 2, 0, 0) + (0,) * (d - 4) + (-1, 1), -3))

    # some other inequalities found by Klyachko
    assert T.is_ressayre(((0, -2, 2, 0) + (0,) * (d - 4) + (-1, 1), -3))
    assert T.is_ressayre(((-2, 0, 2, 0) + (0,) * (d - 4) + (1, -1), -3))
    assert T.is_ressayre(((-1, 1, 1, 0) + (0,) * (d - 4) + (0, 0), -1))
    assert T.is_ressayre(((-2, 1, 0, -1) + (0,) * (d - 4) + (0, -1), -4))
