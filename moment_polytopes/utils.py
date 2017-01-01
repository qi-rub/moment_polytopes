from __future__ import absolute_import, print_function
from sage.all import matrix

__all__ = ['dim_affine_hull']


def dim_affine_hull(points):
    """Return dimension of affine hull of given collection of points."""
    return matrix([p - points[0] for p in points[1:]]).rank()
