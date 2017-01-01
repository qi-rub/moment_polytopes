from __future__ import absolute_import, print_function
import logging
from collections import defaultdict
from sage.all import gcd, vector, matrix, QQ, Polyhedron, cartesian_product, Permutation, Permutations
from . import HRepr, times
from .disk_cache import disk_cache

__all__ = [
    "rect_tableaux",
    "cubicle",
    "cubicle_tableaux",
    "is_dominant",
    "is_extremal_edge",
    "is_extremal_edge_ieq",
    "extremal_edges",
    "perms_of_length",
    "length_tuples",
    "is_shuffle",
    "is_antishuffle",
    "shuffles",
    "antishuffles",
    "perm_action",
    "StabilizerGroup",
]

logger = logging.getLogger(__name__)


@disk_cache
def rect_tableaux(a, b):
    """Return all rectangular standard Young tableaux of shape :math:`a \\times b`.

    :param a: number of rows
    :param b: number of columns
    :rtype: list of :class:`sage.StandardTableau`
    """
    from sage.all import StandardTableaux
    return map(list, StandardTableaux([b] * a))


def cubicle(T):
    """Return Sage :class:`sage.Polyhedron` representing the cubicle that corresponds to the given rectangular tableaux T (if any).
    This is the maximal-dimensional polytope defined by

    .. math::

        \{ (a,b) : a_i + b_j \leq a_k + b_l \Leftrightarrow T_{i,j} \geq T_{k,l} \}

    together with the equations :math:`\sum_i a_i = \sum_j b_j = 0`.

    :param T: a rectangular standard Young tableaux
    :type T: :class:`sage.StandardTableau`
    :rtype: :class:`sage.Polyhedron` or None, if the tableaux is not additive (i.e., does not correspond to a cubicle).
    """
    # assert rectangular shape
    a, b = len(T), len(T[0])
    assert all(len(row) == b for row in T), 'T should be rectangular'

    # set up constraints for a+b-dimensional polyhedron
    # - traces should be zero
    eqns = [
        [0] + [1] * a + [0] * b,
        [0] + [0] * a + [1] * b,
    ]

    # - ordering constraints imposed by Weyl chamber
    ieqs = [[0] + [0] * i + [1, -1] + [0] * (a - 2 - i) + [0] * b
            for i in range(a - 1)]
    ieqs += [[0] + [0] * a + [0] * j + [1, -1] + [0] * (b - 2 - j)
             for j in range(b - 1)]

    # - ordering constraints imposed by cubicle
    for i in range(a):
        for j in range(b):
            for k in range(a):
                for l in range(b):
                    if (i, j) != (k, l):
                        if T[i][j] >= T[k][l]:
                            H = [0] * (a + b)
                            H[i] -= 1
                            H[a + j] -= 1
                            H[k] += 1
                            H[a + l] += 1
                            ieqs.append([0] + H)

    # build sage polyhedron
    P = Polyhedron(ieqs=ieqs, eqns=eqns)

    # compute codimension
    codim = a + b - P.dim()
    assert codim >= 2, 'Expect codimension to be at least two (because of the two trace constraints).'
    return P if codim == 2 else None


@disk_cache
def cubicle_tableaux(a, b):
    """Return list of tableaux corresponding to cubicles for :math:`a \\times b`.

    :param a: number of rows
    :param b: number of columns
    :rtype: list of :class:`sage.StandardTableau`
    """
    return [T for T in rect_tableaux(a, b) if cubicle(T)]


def is_dominant(v):
    r"""Determine if vector :math:`v` is dominant, i.e., :math:`v_0 \geq v_1 \geq \dots`.

    :param v: vector to test.
    :rtype: bool.
    """
    v = list(v)
    return v == sorted(v, reverse=True)


def is_extremal_edge(dims,
                     V,
                     assert_dominant=True,
                     assert_primitive=True,
                     assert_traceless=True):
    """Determine whether given vector is an extremal edge. That is, verify whether :math:`V=(H_1,\dots,H_n)` is

    - dominant
    - primitive
    - traceless
    - admissible for the system :math:`C(d_1,\dots,d_n)` of restricted roots for :math:`SU(d_1) \\times \dots \\times SU(d_n) \\to SU(\prod_i d_i)`

    :param dims: dimensions :math:`d_1,\dots,d_n` of local unitary group.
    :param V: vector of length :math:`\sum_i d_i` to test.
    :param assert_dominant: verify that vector is dominant.
    :param assert_primitive: verify that vector is primitive.
    :param assert_traceless: verify that vector is traceless.
    """
    assert len(V) == sum(dims)

    # extract components and ensure that they are sorted and traceless
    vs = [list(V[sum(dims[:i]):sum(dims[:i + 1])]) for i in range(len(dims))]
    if assert_dominant:
        assert all(is_dominant(v) for v in vs), 'Expected dominant.'
    else:
        vs = [sorted(v, reverse=True) for v in vs]
    if assert_traceless:
        assert all(sum(v) == 0 for v in vs), 'Expect trace to be zero.'

    # check that V is primitive
    if assert_primitive:
        diffs = [x - y for v in vs for (x, y) in zip(v, v[1:])]
        c = gcd(diffs)
        assert c == 1, 'Expect vector that is primitive in dual of root lattice, but gcd is %s' % c

    # build dictionary storing all indices with same sum of components
    indices_by_sum = defaultdict(list)
    for indices in cartesian_product(map(range, dims)):
        s = sum(v[i] for (v, i) in zip(vs, indices))
        indices_by_sum[s].append(indices)

    # build corresponding equations
    # - traces equal to zero
    eqns = []
    for i in range(len(dims)):
        before = sum(dims[:i])
        after = sum(dims[i + 1:])
        eqn = vector([0] * before + [1] * dims[i] + [0] * after)
        eqns.append(eqn)

    # - equality constraints (a_i + b_j + ... == a_k + b_l + ...)
    for indices in indices_by_sum.values():
        for idx1, idx2 in zip(indices, indices[1:]):
            eqn = vector(QQ, sum(dims))
            for i in range(len(dims)):
                offset1 = sum(dims[:i]) + idx1[i]
                offset2 = sum(dims[:i]) + idx2[i]
                eqn[offset1] += 1
                eqn[offset2] -= 1
            eqns.append(eqn)

    # compute dimension of kernel
    codim = sum(dims) - matrix(eqns).rank()
    assert codim >= 1, "Kernel should have dimension at least one (since the given vector is contained in it by construction)."
    return codim == 1


def is_extremal_edge_ieq(dims,
                         ieq,
                         assert_dominant=True,
                         assert_primitive=True,
                         assert_traceless=True):
    """Check whether given inequality corresponds to an extremal edge (see :func:`is_extremal_edge`).

    :param dims: dimensions :math:`d_1,\dots,d_n` of local unitary group.
    :param ieq: inequality :math:`(H,c)` to test. We require that :math:`c=0`.
    :param assert_dominant: verify that vector is dominant.
    :param assert_primitive: verify that vector is primitive.
    :param assert_traceless: verify that vector is traceless.
    """
    H, c = ieq
    assert len(H) == sum(dims)
    assert c == 0

    # extract parts
    hs = [list(H[sum(dims[:i]):sum(dims[:i + 1])]) for i in range(len(dims))]

    # check that last component consists of negative partial sums
    h_last_expected = map(lambda xs: -sum(xs), cartesian_product(hs[:-1]))
    h_last = H[sum(dims[:-1]):]
    assert sorted(h_last) == sorted(h_last_expected)  # return False

    # check if the first components form an extremal edge
    V = H[:sum(dims[:-1])]
    return is_extremal_edge(
        dims[:-1],
        V,
        assert_dominant=assert_dominant,
        assert_primitive=assert_primitive,
        assert_traceless=assert_traceless)


@disk_cache
def _extremal_edges_bipartite(a, b, include_perms=True):
    """Returns extremal edges for :math:`a \\times b`, i.e., the finite list of all vectors :math:`(H_A, H_B)` that are

    - dominant
    - primitive
    - admissible for the system :math:`C(a, b)` of restricted roots for :math:`SU(a) \\times SU(b) \\to SU(ab)`
    """
    # collect extremal edges
    edges = set()
    Ts = rect_tableaux(a, b)
    for i, T in enumerate(Ts):
        if i % 100 == 0:
            logger.debug('Computing extremal edges (%s/%s)', i + 1, len(Ts))
        P = cubicle(T)
        if P:
            edges |= {tuple(ray.vector()) for ray in P.rays()}

    # exclude permutations of the parties?
    if a == b and not include_perms:

        def sort(H):
            H_A, H_B = H[:a], H[a:]
            H_A, H_B = sorted((H_A, H_B))
            return H_A + H_B

        edges = {sort(H) for H in edges}

    # make edges primitive in dual root lattice
    G = times([a, b])
    return map(G.make_dual_root_primitive, edges)


@disk_cache
def _extremal_edges_generic(dims, include_perms=True):
    """Generic implementation that works for any dimension tuple."""
    from sage.all import Subsets, Polyhedron, vector

    # NOTE: An alternative implementation could use solid standard Young tableaux a la http://www.math.rutgers.edu/~zeilberg/mamarim/mamarimhtml/ssyt.html

    # basis vector
    def e(i, d):
        v = [0] * d
        v[i] = 1
        return v

    # weight for given indices
    def omega(indices):
        return vector(sum(map(e, indices, dims), []))

    # compute all restricted roots (up to sign)
    basis = map(tuple, cartesian_product(map(range, dims)))
    restricted_roots = []
    for i, j in Subsets(basis, 2):
        alpha = tuple(omega(i) - omega(j))
        restricted_roots.append(alpha)

    # compute inequalities for Weyl chamber
    ieqs = []
    for i in range(len(dims)):
        for j in range(dims[i] - 1):
            v = [0] * dims[i]
            v[j] = 1
            v[j + 1] = -1
            ieqs.append([0] + [0] * sum(dims[:i]) + v + [0] * sum(dims[i +
                                                                       1:]))
    # choose all possible subsets of equations
    num_eqns = sum(d - 1 for d in dims) - 1
    subsets = Subsets(restricted_roots, num_eqns)
    logger.debug('%s restricted roots => %s subsets each containing %s',
                 len(restricted_roots), subsets.cardinality(), num_eqns)
    rays = set()
    for roots in subsets:
        # trace equations
        eqns = []
        for i in range(len(dims)):
            eqns.append([-1] + [0] * sum(dims[:i]) + [1] * dims[i] + [0] * sum(
                dims[i + 1:]))

        # add orthogonality constraints
        for root in roots:
            eqns.append((0, ) + root)

        # if the space of solutions is one-dimensional, we have found an extremal edge
        P = Polyhedron(ieqs=ieqs, eqns=eqns)
        if P.dim() == 1:
            assert len(P.rays()) == 1
            rays.add(tuple(P.rays()[0]))
    logger.debug("%s distinct extreme rays found", len(rays))

    # remove permutations?
    if not include_perms:
        stab = StabilizerGroup(dims)
        ray_nfs = set()
        for ray in rays:
            H = [
                ray[sum(dims[:i]):sum(dims[:i + 1])] for i in range(len(dims))
            ]
            ray_nf = sum(stab.normal_form(H), ())
            ray_nfs.add(ray_nf)
        rays = ray_nfs

    # post-process edges
    G = times(dims)
    return map(G.make_dual_root_primitive, rays)


def extremal_edges(dims, include_perms=True, algorithm=None):
    """Returns extremal edges for ``dims`` (see :func:`is_extremal_edge`).

    :param dims: dimensions :math:`d_1,\dots,d_n` of local unitary group.
    :param algorithm: ``None``, ``'bipartite'``, or ``'generic'``.
    :rtype: list of :class:`sage.vector`
    """
    # use optimized bipartite implementation (if possible)
    if algorithm is None:
        algorithm = 'bipartite' if len(dims) == 2 else 'generic'
    elif algorithm == 'bipartite':
        assert len(dims) == 2

    if algorithm == 'bipartite':
        return _extremal_edges_bipartite(
            dims[0], dims[1], include_perms=include_perms)
    elif algorithm == 'generic':
        return _extremal_edges_generic(dims, include_perms=include_perms)
    raise Exception('Unknown algorithm "%s"' % algorithm)


def perms_of_length(n, length):
    """Return all permutations in :math:`S_n` of the given length (i.e., with the specified number of inversion).

    This uses the algorithm in `<http://webhome.cs.uvic.ca/~ruskey/Publications/Inversion/InversionCAT.pdf>`_.

    :param n: specifies the permutation group :math:`S_n`.
    :param length: number of inversions.
    :rtype: list of :class:`sage.Permutation`
    """
    result = []

    def gen(S, l, suffix=[]):
        if l == 0:
            result.append(Permutation(S + suffix))
            return

        n = len(S)
        bin = (n - 1) * (n - 2) / 2
        for i in range(n):
            if n - (i + 1) <= l <= bin + n - (i + 1):
                x = S[i]
                gen(S[0:i] + S[i + 1:], l - n + (i + 1), [x] + suffix)

    gen(S=range(1, n + 1), l=length)
    return result


def length_tuples(dims, total):
    r"""Return integer tuples :math:`(\ell_1, ..., \ell_n)` such that

    - each component :math:`\ell_i` is in :math:`\{0,\dots,{d_i \choose 2}\}`
    - their sum is equal to ``total``

    :param dims: dimensions :math:`d_1,\dots,d_n`.
    :param total: total length
    :rtype: generator of tuples of integers
    """
    # compute maximal lengths
    max_lengths = [d * (d - 1) / 2 for d in dims]
    ranges = [range(0, min(m + 1, total + 1)) for m in max_lengths[:-1]]
    for most in cartesian_product(ranges):
        last = total - sum(most)
        if last >= 0 and last <= max_lengths[-1]:
            yield tuple(most) + (last, )


def is_shuffle(pi, v):
    """Check if the permutation ``pi`` is a shuffle with respect to the dominant element ``v``, i.e.,

    .. math::

        v_i = v_{i+1} \Rightarrow \pi_i < pi_{i+1}.

    :param pi: the permutation :math:`pi`.
    :param v: the dominant vector :math:`v`.
    :rtype: bool
    """
    assert is_dominant(v)
    return all(pi[n] < pi[n + 1] for n in range(len(v) - 1)
               if v[n] == v[n + 1])


def is_antishuffle(pi, v):
    """Check if the permutation ``pi`` is an antishuffle with respect to the dominant element ``v``, i.e.,

    .. math::

        v_i = v_{i+1} \Rightarrow \pi_i > pi_{i+1}.

    :param pi: the permutation :math:`pi`.
    :param v: the dominant vector :math:`v`.
    :rtype: bool
    """
    assert is_dominant(v)
    return all(pi[n] > pi[n + 1] for n in range(len(v) - 1)
               if v[n] == v[n + 1])


def shuffles(v, length):
    r"""Return all permutations in :math:`S_{\lvert v \rvert}` that are shuffles with respect to the dominant element ``v`` (see :func:`is_shuffle`) and have the desired length.

    :param v: the dominant vector :math:`v`.
    :param length: the desired length.
    :rtype: :class:`sage.Permutation`
    """
    result = []
    assert is_dominant(v)

    def gen(S, l, suffix=[]):
        if not S and l == 0:
            result.append(Permutation(S + suffix))
            return

        n = len(S)
        bin = (n - 1) * (n - 2) / 2
        for i in range(n):
            if n - (i + 1) <= l <= bin + n - (i + 1):
                x = S[i]

                # filter non-shuffles
                if suffix and v[n - 1] == v[n]:
                    # if pi(n) > pi(n + 1)
                    if x > suffix[0]:
                        continue

                gen(S[0:i] + S[i + 1:], l - n + (i + 1), [x] + suffix)

    gen(S=range(1, len(v) + 1), l=length)
    return result


def antishuffles(v, antilength):
    r"""Return all permutations in :math:`S_{\lvert v \rvert}` that are antishuffles with respect to the dominant element ``v`` (see :func:`is_antishuffle`) and have the desired antilength.

    :param v: the dominant vector :math:`v`.
    :param antilength: the desired antilength.
    :rtype: :class:`sage.Permutation`
    """
    d = len(v)
    pis = shuffles(v, length=antilength)
    return [Permutation([d + 1 - pi[i] for i in range(d)]) for pi in pis]


def perm_action(pi, v, zero_based=False):
    r"""Left action of a permutation :math:`\pi \in S_n` on :math:`v \in \mathbb R^n`:

    .. math::

        (\pi \cdot v)_i = v_{\pi^{-1}(i)}

    i.e.

    .. math::

        (\pi \cdot v)_{\pi(i)} = v_i.

    :param pi: the permutation :math:`\pi`.
    :param v: the vector :math:`v`.
    :para, zero_based: if ``True`` then permutation is in :math:`S_{\{0,\dots,n-1\}}` instead.
    :rtype: :class:`sage.vector`
    """
    # make zero-based
    pi = list(pi) if zero_based else [x - 1 for x in pi]
    assert len(pi) == len(v) and sorted(pi) == range(len(v))

    # permute components
    d = len(pi)
    v_out = [0] * d
    for i in range(d):
        v_out[pi[i]] = v[i]
    return tuple(v_out)


class StabilizerGroup(object):
    r"""Stabilizer group :math:`S_v \subseteq S_{\lvert v \rvert}` of a vector :math:`v`.

    :param v: the vector.
    """

    def __init__(self, v):
        #: The vector defining the stabilizer group.
        self.v = v

        #: The blocks of indices of components that are equal.
        self.blocks = defaultdict(list)
        for k, d in enumerate(v):
            self.blocks[d].append(k)
        self.blocks = self.blocks.values()

    def normal_form(self, hs):
        """Returns the unique normal form of a vector ``hs`` in its orbit under the stabilizer group.

        :param hs: a vector of the same length as :attr:`v`.
        :rtype: tuple
        """
        assert len(hs) == len(self.v)
        hs_nf = [None] * len(self.v)
        for indices in self.blocks:
            # hs[indices] = sorted(hs[indices])
            group = sorted([hs[k] for k in indices])
            for j, k in enumerate(indices):
                hs_nf[k] = group[j]
        return tuple(hs_nf)

    def orbit(self, hs_iterable):
        """Returns the orbit of a vectors ``hs_iterable`` under the stabilizer group.

        :param hs_iterable: a collection of vectors
        :rtype: set of tuples
        """
        # generate all permutations of indices
        blocks_perms = cartesian_product(map(Permutations, self.blocks))

        # apply all permutations
        hs_perm = [None] * len(self.v)
        orbit = set()
        for blocks_perm in blocks_perms:
            for hs in hs_iterable:
                # permute according to (idx1 -> idx2) in each block
                for (idx1, idx2) in zip(self.blocks, blocks_perm):
                    for (i, j) in zip(idx1, idx2):
                        hs_perm[j] = hs[i]
                orbit.add(tuple(hs_perm))
        return orbit
