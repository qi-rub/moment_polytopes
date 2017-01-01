from __future__ import absolute_import, print_function
from sage.all import vector, QQ
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
