from __future__ import absolute_import, print_function
import logging, string
from sage.all import QQ, lcm, gcd, vector, prod, cartesian_product
from . import extremal_edges, StabilizerGroup, external_tensor_product, c_candidates, length_tuples, antishuffles, perm_action, ressayre_tester, HRepr
from .disk_cache import disk_cache

__all__ = [
    'DEFAULT_SUBSYSTEM_LABELS',
    'H_AB_dominant',
    'H_ABC_dominant',
    'H_dominant_admissible',
    'H_candidates',
    'H_ressayre',
    'hrepr',
    'vrepr',
    'facet_normal_form',
    'ieqs_wo_perms',
    'vertices_wo_perms',
    'pretty',
]

logger = logging.getLogger(__name__)

#: Default subsystem labels used by :func:`pretty`.
DEFAULT_SUBSYSTEM_LABELS = string.ascii_uppercase


def H_AB_dominant(a, b, include_perms=True):
    """Candidates for dominant and primitive :math:`(H_A, H_B)`, i.e., extremal edges.

    :param a: dimension of first tensor factor.
    :param b: dimension of second tensor factor.
    :param include_perms: if ``True``, include permutations of the two subsystems.
    :rtype: set of tuples :math:`(H_A,H_B)`.
    """
    edges = extremal_edges((a, b), include_perms=include_perms)
    return {(tuple(H[:a]), tuple(H[a:])) for H in edges}


@disk_cache
def H_ABC_dominant(a, b, c, include_perms=True):
    """Candidates for dominant and primitive :math:`(H_A, H_B, H_C)`.

    :param a: dimension of first tensor factor.
    :param b: dimension of second tensor factor.
    :param c: dimension of third tensor factor.
    :param include_perms: if ``True``, include permutations of the three subsystems.
    :rtype: set of tuples :math:`(H_A,H_B,H_C)`.
    """
    dims = [a, b, c]
    stab = StabilizerGroup(dims)

    # zero vectors
    zero_A, zero_B, zero_C = [(0, ) * d for d in dims]

    # fetch extremal edges
    H_ABs = H_AB_dominant(a, b, include_perms=False)
    H_ACs = H_AB_dominant(a, c, include_perms=False) | {(zero_A, zero_C)}
    H_BCs = H_AB_dominant(b, c, include_perms=False) | {(zero_B, zero_C)}
    H_Cs = {H_C for (H_A, H_C) in H_ACs} & {H_C for (H_B, H_C) in H_BCs}

    # add candidates with (H_A, H_B) != 0
    candidates = set()
    for (H_A, H_B) in H_ABs:
        for H_C in H_Cs:
            if (H_A, H_C) in H_ACs and (H_B, H_C) in H_BCs:
                candidates.add(stab.normal_form((H_A, H_B, H_C)))

    # add candidates with (H_A, H_B) == 0
    for H_C in H_Cs:
        if H_C != zero_C and (zero_A, H_C) in H_ACs and (zero_B, H_C) in H_BCs:
            candidates.add(stab.normal_form((zero_A, zero_B, H_C)))

    # include permutations?
    if include_perms:
        candidates = stab.orbit(candidates)

    return candidates


@disk_cache
def H_dominant_admissible(dims, include_perms=True):
    r"""Candidates for dominant and admissible :math:`((H_A, H_B, H_C,\dots),c)`.

    :param dims: dimensions :math:`d_1,\dots,d_n` of the tensor factors.
    :param include_perms: if ``True``, include permutations of the :math:`n` subsystems.
    :rtype: set of tuples :math:`((H_A,H_B,H_C,\dots),c)`.
    """
    assert not include_perms
    assert list(dims) == sorted(
        dims), 'Dimensions should be sorted increasingly: %s' % str(dims)

    # HEURISTICS ONE: extremal edges, partial sums, z = 0
    if prod(dims[:-1]) == dims[-1]:
        candidates = set()
        for H in extremal_edges(dims[:-1], include_perms=False):
            # build chunks and extend by the vector of partial sums
            most = [
                tuple(H[sum(dims[:i]):sum(dims[:i + 1])])
                for i in range(len(dims) - 1)
            ]
            last = tuple(
                sorted(
                    [-sum(entries) for entries in cartesian_product(most)],
                    reverse=True))
            candidates.add((tuple(most + [last]), 0))

        # do not forget to add the *dominant*, *traceless* version of lambda_{AB...,ab...} >= 0 as another constraint
        #hs_positivity_dominant = [(0,) * dim for dim in dims[:-1]] + [(1,) + (0,) * (dims[-1] - 1)]
        hs_positivity_dominant = [(0, ) * dim
                                  for dim in dims[:-1]] + [(dims[-1] - 1, ) +
                                                           (-1, ) *
                                                           (dims[-1] - 1)]
        candidates.add((tuple(hs_positivity_dominant), -1))
        return candidates

    # HEURISTICS TWO: tripartite via bipartite
    if len(dims) == 3:
        R = external_tensor_product(dims)
        candidates = set()
        for hs in H_ABC_dominant(*dims, include_perms=False):
            for z in c_candidates(R, sum(hs, ())):
                candidates.add((hs, z))
        return candidates

    raise Exception("No suitable heuristics for %s." % (dims, ))


@disk_cache
def H_candidates(dims, include_perms=True):
    """Return candidates for Ressayre elements (all conditions except determinant condition).

    :param dims: dimensions :math:`d_1,\dots,d_n` of the tensor factors.
    :param include_perms: if ``True``, include permutations of the :math:`n` subsystems.
    :rtype: set of tuples :math:`((H_A,H_B,H_C,\dots),c)`.
    """
    if include_perms:
        raise NotImplementedError
    R = external_tensor_product(dims)
    stab = StabilizerGroup(dims)
    antilength_max = len(R.negative_roots)

    # for all dominant admissible H
    candidates = set()
    for (hs_dominant, z) in H_dominant_admissible(dims, include_perms=False):
        # compute desired antilength
        H_dominant = vector(sum(hs_dominant, ()))
        antilength_desired = sum(1 for omega in R.weights
                                 if omega.dot_product(H_dominant) < z)
        if antilength_desired > antilength_max:
            continue

        # determine possible tuples of antilengths (up to permutations of the hs)
        hs_stab = StabilizerGroup(hs_dominant)
        antilength_tuples = {
            hs_stab.normal_form(tuple)
            for tuple in length_tuples(dims, total=antilength_desired)
        }

        # determine possible tuples of antishuffles (up to permutations of the hs)
        antishuffle_tuples = set()
        for antilength_tuple in antilength_tuples:
            for antishuffle_tuple in cartesian_product(
                    map(antishuffles, hs_dominant, antilength_tuple)):
                antishuffle_tuples.add(hs_stab.normal_form(antishuffle_tuple))

        # act by each tuple of antishuffles, and store result (up to permutations of the subsystems)
        for antishuffle_tuple in antishuffle_tuples:
            hs = tuple(map(perm_action, antishuffle_tuple, hs_dominant))
            hs = stab.normal_form(hs)
            candidates.add((hs, z))

    return candidates


@disk_cache
def H_ressayre(dims, include_perms=True, **kwargs):
    """Return all Ressayre elements for representation of :math:`\times_i GL(d_i)` on :math:`\bigotimes_i \mathbb C^{d_i}`.

    :param dims: dimensions :math:`d_1,\dots,d_n` of the tensor factors.
    :param include_perms: if ``True``, include permutations of the :math:`n` subsystems.

    All other arguments are forwarded to :func:`moment_polytopes.ressayre_tester`.
    """
    assert not include_perms
    R = external_tensor_product(dims)
    T = ressayre_tester(R, **kwargs)

    # test all candidates inequalities...
    ieqs = set()
    candidates = H_candidates(dims, include_perms=False)

    for i, (hs, z) in enumerate(candidates):
        logger.debug('checking candidate inequality (%d/%d)', i + 1,
                     len(candidates))
        H = sum(hs, ())
        if T.is_ressayre((H, z)):
            ieqs.add((hs, z))

    return ieqs


@disk_cache
def hrepr(dims, irred=True, **kwargs):
    r"""Return H-representation of moment polytope for representation of :math:`\times_i GL(d_i)` on :math:`\bigotimes_i \mathbb C^{d_i}`.

    :param dims: dimensions :math:`d_1,\dots,d_n` of the tensor factors.
    :param irred: if ``True`` then an irredunant H-representation is returned.

    All other arguments are forwarded to :func:`moment_polytopes.ressayre_tester`.

    :rtype: :class:`moment_polytopes.HRepr`
    """
    stab = StabilizerGroup(dims)
    R = external_tensor_product(dims)
    ieqs_flat = []
    for (hs, z) in H_ressayre(dims, include_perms=False, **kwargs):
        for hs_permuted in stab.orbit([hs]):
            H_permuted = sum(hs_permuted, ())
            ieqs_flat.append((H_permuted, z))

    # build H-representation
    hrepr = HRepr(ieqs=ieqs_flat)
    hrepr = hrepr & R.reduced_positive_weyl_chamber_hrepr
    return hrepr.irred() if irred else hrepr


@disk_cache
def vrepr(dims, **kwargs):
    r"""Return V-representation of moment polytope for representation of :math:`\times_i GL(d_i)` on :math:`\bigotimes_i \mathbb C^{d_i}`.

    :param dims: dimensions :math:`d_1,\dots,d_n` of the tensor factors.
    :param irred: if ``True`` then an irredunant H-representation is returned.

    All other arguments are forwarded to :func:`moment_polytopes.ressayre_tester`.

    :rtype: :class:`moment_polytopes.VRepr`
    """
    return hrepr(dims, irred=True, **kwargs).vrepr()


def facet_normal_form(dims, ieq):
    r"""Given a facet :math:`H \cdot \lambda \geq c` for the quantum marginal problem with given dimensions, where :math:`H = (H_A, H_B, \dots)`, make each component :math:`H_A, H_B, \dots` traceless, integral, and primitive.

    This generically makes the inequality unique.

    :param dims: the dimensions :math:`d_1,\dots,d_n`.
    :param ieq: the inequality :math:`(H,c)` defining the facet.
    """
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


def ieqs_wo_perms(dims, ieqs):
    """Returns list of inequalities up to permutations of subsystems.

    :param dims: the dimensions :math:`d_1,\dots,d_n`.
    :param ieqs: list of inequalities :math:`(H,c)`.
    """
    assert sum(dims) == len(ieqs[0][0]), 'Total dimension mismatch.'
    stab = StabilizerGroup(dims)
    result = set()
    for (H, z) in ieqs:
        # extract and verify that traceless
        hs = [
            tuple(H[sum(dims[:i]):sum(dims[:i + 1])])
            for i in range(len(dims))
        ]
        assert all(sum(h) == 0 for h in hs), 'Not traceless: %s' % (hs, )

        # bring into normal form, compress, and check that we have integers with gcd = 1
        hs_nf = stab.normal_form(hs)
        H_nf = sum(hs_nf, ())
        assert all(QQ(int(x)) == x for x in H_nf), 'not integers'
        f = gcd([x for x in H_nf] + [z])
        assert f == 1, 'not simplified (%s)' % f
        result.add((H_nf, z))
    return result


def vertices_wo_perms(dims, vertices):
    """Returns list of vertices up to permutations of subsystems.

    :param dims: the dimensions :math:`d_1,\dots,d_n`.
    :param vertices: list of vertices.
    """
    assert sum(dims) == len(vertices[0]), 'Total dimension mismatch.'
    stab = StabilizerGroup(dims)
    result = set()
    for V in vertices:
        # extract
        vs = [
            tuple(V[sum(dims[:i]):sum(dims[:i + 1])])
            for i in range(len(dims))
        ]

        # bring into normal form and compress
        vs_nf = stab.normal_form(vs)
        V_nf = sum(vs_nf, ())
        result.add(V_nf)
    return result


class PrettyPrinter(object):
    """Pretty-print moment polytope for pure-state quantum marginal problem.

    :param dims: dimensions :math:`d_1,\dots,d_n` of the tensor factors.
    :param show_hrepr: show H-representation.
    :param show_vrepr: show V-representation.
    :param include_perms: if ``True``, include permutations of the :math:`n` subsystems.
    :param subsystem_labels: custom subsystem labels.

    All other arguments are forwarded to :func:`moment_polytopes.ressayre_tester`.
    """

    def __init__(self,
                 dims,
                 show_hrepr=True,
                 show_vrepr=True,
                 include_perms=False,
                 subsystem_labels=None,
                 **kwargs):
        self.dims = dims
        self.subsystem_labels = subsystem_labels if subsystem_labels else DEFAULT_SUBSYSTEM_LABELS

        # compute H-representation
        if show_hrepr:
            ieqs = hrepr(dims, **kwargs).ieqs
            ieqs = ieqs if include_perms else ieqs_wo_perms(dims, ieqs)
            self.ieqs = sorted((tuple(H), c) for (H, c) in ieqs)
        else:
            self.ieqs = None

        # compute V-representation
        if show_vrepr:
            vertices = vrepr(dims, **kwargs).vertices
            vertices = vertices if include_perms else vertices_wo_perms(
                dims, vertices)
            self.vertices = sorted(tuple(V) for V in vertices)
        else:
            self.vertices = None

    def __repr__(self):
        """Pretty-print quantum marginal problem polytope."""
        from tabulate import tabulate as _tabulate

        # title
        title = 'C(%s)' % ','.join(str(d) for d in self.dims)
        lines = [title, '=' * len(title)]

        # facets
        if self.ieqs is not None:
            headers = ['#'] + [
                'H_%s' % s for (_, s) in zip(self.dims, self.subsystem_labels)
            ] + ['z', 'Remarks']
            lines += [
                '', 'Facets', '------', '', _tabulate(
                    self._facets_table(), headers=headers)
            ]
            lines += [
                '',
                'Facet format is (H_A,lambda_A) + ... + z >= 0. The last column states',
                'whether the facet includes the origin (o) or the highest weight (*).'
            ]

        # facets
        if self.vertices is not None:
            headers = ['#'] + [
                'V_%s' % s for (_, s) in zip(self.dims, self.subsystem_labels)
            ]
            lines += [
                '', 'Vertices', '--------', '', _tabulate(
                    self._vertices_table(), headers=headers)
            ]

        lines += ['', 'All data is up to permutations of subsystems.']
        return '\n'.join(lines)

    def _facets_table(self):
        facets = []
        for idx, (H, c) in enumerate(self.ieqs):
            hs = [
                H[sum(self.dims[:i]):sum(self.dims[:i + 1])]
                for i in range(len(self.dims))
            ]

            # collect remarks
            assert all(sum(h) == 0 for h in hs)
            remarks = []
            if c == 0: remarks.append('o')
            if sum(h[0] for h in hs) == c: remarks.append('*')

            # format
            hs = ['(%s)' % ', '.join(map(str, h)) for h in hs]
            facets.append([idx + 1] + hs + [-c, ', '.join(remarks)])
        return facets

    def _vertices_table(self):
        vertices = []
        for idx, V in enumerate(self.vertices):
            vs = [
                V[sum(self.dims[:i]):sum(self.dims[:i + 1])]
                for i in range(len(self.dims))
            ]
            vs = ['(%s)' % ', '.join(map(str, v)) for v in vs]
            vertices.append([idx + 1] + vs)
        return vertices


pretty = PrettyPrinter
