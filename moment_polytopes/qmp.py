from __future__ import absolute_import, print_function
from sage.all import QQ, lcm, gcd, vector

__all__ = ['facet_normal_form']


def facet_normal_form(dims, ieq):
    r"""Given a facet :math:`H \cdot \lambda \geq c` for the quantum marginal problem with given dimensions, where :math:`H = (H_A, H_B, \dots)`, make each component :math:`H_A, H_B, \dots` traceless, integral, and primitive.

    This generically makes the inequality unique.

    :param dims: the dimensions :math:`d_1,\dots,d_n`.
    :param ieq: the inequality :math:`(H,c)` defining the facet.
    :rtype: the """
    H, c = ieq
    assert sum(dims) == len(H)
    hs = [H[sum(dims[:i]):sum(dims[:i + 1])] for i in range(len(dims))]

    # (1,...,1) with sum "a" corresponds to a shift of -1 on z (since z is on the right-hand side of the equations)
    subs = [QQ((-sum(h), d)) for (h, d) in zip(hs, dims)]
    H = []
    for (h, s, d) in zip(hs, subs, dims):
        H += list(vector(h) + vector([s] * d))
    H = vector(H)
    c += sum(subs)

    # make all components integral with gcd = 1
    f = lcm([x.denominator() for x in H] + [c.denominator()])
    H, c = f * H, f * c
    f = gcd([x for x in H] + [c])
    if f:
        H, c = H / f, c / f
    return H, c
