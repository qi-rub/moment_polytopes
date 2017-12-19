from __future__ import absolute_import, print_function
import pytest
from moment_polytopes import *


# prior work
def test_borland_dennis(algorithm):
    R = weyl_module(6, [1, 1, 1])
    C = ressayre_tester(R, algorithm=algorithm)
    assert C.is_ressayre(((0, 0, 0, -1, 1, 1), 0))


def test_8_qubits(algorithm):
    R = external_tensor_product([2] * 8)
    C = ressayre_tester(R, algorithm=algorithm)

    for i in range(8):
        assert C.is_ressayre(([-1, 0] * i + [1, 0] + [-1, 0] * (8 - 1 - i),
                              2 - 8))


def test_bravyi(algorithm):
    R = external_tensor_product([2, 2, 4])
    C = ressayre_tester(R, algorithm=algorithm)

    assert C.is_ressayre(((-1, 0, 0, 0, 1, 1, 0, 0), 0))
    assert C.is_ressayre(((0, 0, -1, 0, 1, 1, 0, 0), 0))
    assert C.is_ressayre(((-1, 0, -1, 0, 1, 0, 0, -1), -1))
    assert C.is_ressayre(((1, 0, -1, 0, 1, 0, -1, 0), 0))
    assert C.is_ressayre(((1, 0, -1, 0, 0, 1, 0, -1), 0))
    assert C.is_ressayre(((-1, 0, 1, 0, 1, 0, -1, 0), 0))
    assert C.is_ressayre(((-1, 0, 1, 0, 0, 1, 0, -1), 0))


# examples in my thesis
def test_qubits(algorithm):
    # 7 qubits
    R = external_tensor_product([2] * 7)
    C = ressayre_tester(R, algorithm=algorithm)
    assert C.is_ressayre(([-1, 1] * 6 + [1, -1], 2 - 7))

    # 3 qubits -- redundant inequality!
    R = external_tensor_product([2, 2, 2])
    C = ressayre_tester(R, algorithm=algorithm)
    assert C.is_ressayre(((0, 0, -1, 1, 0, 0), -1))


def test_mixed_state_of_two_qubits(algorithm):
    R = external_tensor_product([2, 2, 4])
    C = ressayre_tester(R, algorithm=algorithm)

    assert C.is_ressayre(((-1, 1, 1, -1, 2, 0, -2, 0), 0))
    assert C.is_ressayre(((-1, 1, 1, -1, 0, 2, 0, -2), 0))


def test_three_qutrits(algorithm):
    R = external_tensor_product([3, 3, 3])
    C = ressayre_tester(R, algorithm=algorithm)

    assert C.is_ressayre(((0, -1, 1, -1, 0, 1, 1, 0, -1), -1))


@pytest.mark.parametrize("d", [6, 7, 8])
def test_three_fermions_with_total_spin_one_half(algorithm, d):
    R = external_tensor_product([weyl_module(d, [2, 1]), 2])
    C = ressayre_tester(R, algorithm=algorithm)

    assert C.is_ressayre(([-2, 2] + [0] * (d - 2) + [-1, 1], -3))
