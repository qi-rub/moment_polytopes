from __future__ import absolute_import, print_function
from sage.all import Integer, vector, gcd, ZZ, QQ, RootSystem

__all__ = ["is_dual_root_primitive", "dual_root_primitive"]


def _dual_root_primitive_gcd(root_system, H):
    """Return greatest common divisor of coefficients of H with respect to dual roots."""
    # get ambient space
    ambient_space = RootSystem(root_system).ambient_space()
    assert len(H) == ambient_space.dimension(), 'Dimension mismatch'

    # compute coefficients
    H = vector(H)
    cs = [
        vector(alpha).dot_product(H) for alpha in ambient_space.simple_roots()
    ]
    assert all(c in ZZ for c in cs), 'Vector is not in the dual root lattice.'

    # divide by gcd
    return gcd(cs)


def is_dual_root_primitive(root_system, H):
    """Determine if the vector ``H`` is primitive in the dual root lattice of the given root system.

    :param root_system: root system.
    :param H: vector to be rescaled. Should be an element of the dual root lattice (in ambient coordinates).
    :type root_system: :class:`sage.RootSystem`
    :type H: :class:`sage.vector`
    :rtype: :class:`sage.vector`
    """
    return _dual_root_primitive_gcd(root_system, H) in [1, -1]


def dual_root_primitive(root_system, H):
    """Rescale the vector ``H`` so that it is primitive in the dual root lattice of the given root system.

    :param root_system: root system.
    :param H: vector to be rescaled. Should be an element of the dual root lattice (in ambient coordinates).
    :type root_system: :class:`sage.RootSystem`
    :type H: :class:`sage.vector`
    :rtype: :class:`sage.vector`
    """
    c = _dual_root_primitive_gcd(root_system, H)
    return vector(H) / c
