from __future__ import absolute_import, print_function
import logging
from collections import defaultdict
from sage.all import vector, matrix, ceil, log, random_vector, ZZ, QQ, PolynomialRing, mathematica
from .utils import dim_affine_hull

__all__ = [
    'DEFAULT_FAILURE_PROBABILITY',
    'RessayreTester',
    'ressayre_tester',
    'is_ressayre',
    'is_admissible',
    'c_candidates',
]

logger = logging.getLogger(__name__)

#: Default failure probability of ``probabilistic`` algorithm for testing Ressayre elements (see :func:`ressayre_tester`).
DEFAULT_FAILURE_PROBABILITY = 1e-10


class RessayreTester(object):
    """Base class for testing Ressayre elements. Use :func:`is_ressayre` and :func:`ressayre_tester`.

    :param R: the representation.
    :type R: :class:`Representation`
    """

    def __init__(self, R):
        #: The representation.
        self.representation = R

    def is_ressayre(self, ieq):
        """Check whether given :math:`(H,c)` is a Ressayre element (see :func:`is_ressayre`).

        :param ieq: the inequality :math:`(H,c)` to be tested.
        """
        """Check hypothetical inequality (-,H) >= c."""
        R = self.representation
        H, c = vector(ieq[0]), ieq[1]

        # determine n_-(H < 0) and group according to H-weight
        n_neg = defaultdict(list)
        for i, alpha in enumerate(R.negative_roots):
            weight = alpha.dot_product(H)
            if weight < 0:
                n_neg[weight].append(i)

        # determine M(H < c) and group according to H-weight
        M_neg = defaultdict(list)
        for j, omega in enumerate(R.weights):
            weight = omega.dot_product(H)
            if weight < c:
                M_neg[weight - c].append(j)

        # weight blocks mismatch?
        if sorted(n_neg.keys()) != sorted(M_neg.keys()):
            logger.debug("H-weights do not match: %r != %r", n_neg.keys(),
                         M_neg.keys())
            return False
        for weight in n_neg:
            if len(n_neg[weight]) != len(M_neg[weight]):
                logger.debug("H-weight %r multiplicity mismatch: %d != %d",
                             weight, len(n_neg[weight]), len(M_neg[weight]))
                return False

        # determine M(H = c)
        M_null = [
            j for j, omega in enumerate(R.weights) if omega.dot_product(H) == c
        ]
        if not M_null:
            logger.debug("H-weight space %d is empty", c)
            return False

        # proceed block by block...
        logger.debug('need to compute %d determinants...', len(n_neg))
        for weight, n_neg_block in n_neg.iteritems():
            # only nonzero blocks without mismatches should appear
            M_neg_block = M_neg[weight]
            assert len(M_neg_block) == len(n_neg_block) > 0

            # compute tangent map for each basis vector of M(H = c)
            P_neg = matrix([R.weight_vector(j) for j in M_neg_block])
            Ts = {
                j: matrix(
                    [
                        P_neg * R.negative_root_action(i, j)
                        for i in n_neg_block
                    ],
                    sparse=True)
                for j in M_null
            }

            # compute determinant
            logger.debug(
                'processing H-weight block %d -> %d: size %dx%d with %d variables',
                weight, weight + c, len(n_neg_block), len(n_neg_block),
                len(M_null))

            # compute determinant of polynomial matrix
            det = self.det(Ts, len(n_neg_block))
            logger.debug('... det() returned %s', det)
            if not det:
                return False
        logger.debug('all %d determinants nonzero', len(n_neg))
        return True

    def det(self, Ts, d):
        """Return nonzero value if and only if determinant of polynomial matrix :math:`sum_j z_j T_j` is nonzero.

        :param Ts: quadratic matrices :math:`T_j`, indexed by weight vector index :math:`j`.
        :param d: size of the matrices :math:`T_j`
        """
        raise NotImplementedError


class ProbabilisticRessayreTester(RessayreTester):
    """Probabilistic algorithm for testing Ressayre elements, based on the Schwartz-Zippel lemma.

    If :meth:`det` returns nonzero then the symbolic determinant is guaranteed to be nonzero.
    However, :meth:`det` may fail and zero even though the symbolic determinant is nonzero.

    :param failure_probability: ``None``, or desired probability of failure.
    """

    def __init__(self, R, failure_probability=None):
        super(ProbabilisticRessayreTester, self).__init__(R)

        #: The failure probability.
        self.failure_probability = failure_probability if failure_probability is not None else DEFAULT_FAILURE_PROBABILITY

    def det(self, Ts, d):
        # the determinant polynomial is of degree d. therefore, the Schwartz-Zippel lemma asserts that by choosing a random
        # number in {0, ..., N*d - 1} we obtain a failure probability of <= 1/N
        # => we repeat this process a sufficient number of times to achieve the desired probability of failure
        # XXX: we might want to keep N * d bounded to avoid expensive arithmetic
        N = 16
        num_repetitions = ceil(-log(self.failure_probability) / log(N))
        logger.debug(
            '%d repetitions to achieve desired failure probability %r',
            num_repetitions, self.failure_probability)

        for _ in range(num_repetitions):
            # draw random base point zs in M(H - c = 0)
            zs = random_vector(ZZ, self.representation.dimension, N * d)

            # build tangent map n_-(H < 0 | weight) ----> M(H - c < 0 | weight) at that point
            T = sum(T * zs[j] for (j, T) in Ts.iteritems())

            # return if non-zero (otherwise keep trying)
            det = T.determinant()
            assert det in QQ, 'Assuming that integer determinants are computed exactly.'
            if det:
                return det
        return 0


class SageRessayreTester(RessayreTester):
    def __init__(self, R):
        super(SageRessayreTester, self).__init__(R)

        PR = PolynomialRing(ZZ, len(R.weights), b"z")
        self.zs = PR.gens()
        """Polynomial variables for each weight vector."""

    def det(self, Ts, d):
        T = sum(T * self.zs[j] for (j, T) in Ts.iteritems())
        return T.determinant()


class MathematicaRessayreTester(RessayreTester):
    def __init__(self, R):
        super(MathematicaRessayreTester, self).__init__(R)

        PR = PolynomialRing(ZZ, len(R.weights), b"z")
        self.zs = PR.gens()
        """Polynomial variables for each weight vector."""

    def det(self, Ts, d):
        T = sum(T * self.zs[j] for (j, T) in Ts.iteritems())
        return mathematica(T).Det().sage()


class CompositeRessayreTester(RessayreTester):
    """Test Ressayre elements by deferring tries multiple :class:`RessayreTester` until one reports a non-zero determinant.

    Useful for combining fast probabilistic testers with a deterministic tester.
    """

    def __init__(self, testers):
        assert len(testers) > 0, 'Need to specify at least one tester.'
        R = testers[0].representation
        assert all(
            T.representation is R for T in
            testers), 'All testers should work on the same representation.'
        super(CompositeRessayreTester, self).__init__(R)

        #: The delegate testers.
        self.testers = testers

    def det(self, Ts, d):
        for T in self.testers:
            det = T.det(Ts, d)
            if det:
                return det
        return 0


def ressayre_tester(R, algorithm=None, failure_probability=None):
    r"""Create object for batch-testing Ressayre elements (see :func:`is_ressayre`).

    A **Ressayre element** :math:`(H,c)` for a :math:`G`-representation :math:`R` satisfies the following three properties:

    1. the hyperplane :math:`H \cdot \omega = c` is spanned by weights (admissibility, see :func:`is_admissible`),
    2. the number of negative roots :math:`\alpha` with :math:`H \cdot \alpha < c` is equal to the number of weights :math:`\varphi` of the representation :math:`R` such that :math:`H \cdot \varphi < c`,
    3. for some weight vector with weight on the hyperplane :math:`H \cdot \omega = c`, the correspondingly restricted tangent map is an isomorphism.

    By `Vergne and Walter (2014) <https://arxiv.org/abs/1410.8144>`_, the Ressayre elements form a complete set of inequalities moment polytope for the :math:`G`-action on the projective space :math:`\mathbb P(R)` (except for the Weyl chamber constraints).

    :param R: the representation.
    :param algorithm: ``None``, ``'sage'``, ``'mathematica'``, or ``'probabilistic'``.
    :param failure_probability: ``None``, or desired probability of failure. This only affects the correctness when using the ``'probabilistic'`` algorithm.
    :type R: :class:`Representation`
    :rtype: :class:`RessayreTester`
    """
    # default algorithm
    if algorithm is None:
        algorithm = 'sage'

    # first use probabilistic algorithm
    T = ProbabilisticRessayreTester(R, failure_probability)
    if algorithm == 'probabilistic':
        return T

    # create RessayreTester instance
    if algorithm == 'sage':
        return CompositeRessayreTester([T, SageRessayreTester(R)])
    elif algorithm == 'mathematica':
        return CompositeRessayreTester([T, MathematicaRessayreTester(R)])
    raise NotImplementedError('Unknown algorithm "%s".' % algorithm)


def is_ressayre(R, ieq, **kwargs):
    """Determine whether given inequality :math:`(H,c)` is a Ressayre element for the given representation.

    Convenience function that accepts the same optional argument as :func:`ressayre_tester`.

    :param R: the representation.
    :param ieq: the inequality :math:`(H,c)` to be tested.
    :type R: :class:`Representation`
    """
    return ressayre_tester(R, **kwargs).is_ressayre(ieq)


def is_admissible(R, ieq):
    r"""Determine whether the given element :math:`(H,c)` is admissible (i.e., that the hyperplane :math:`H \cdot \omega = c` is spanned by weights of the representation).

    :param R: the representation.
    :param ieq: the element :math:`(H,c)`.
    :type R: :class:`Representation`
    """
    # extract element
    H, c = ieq
    H = vector(H)

    # determine weights on the hyperplane
    weights = [omega for omega in R.weights if omega.dot_product(H) == c]

    # impossible?
    if len(weights) < R.dimension_affine_hull_weights:
        logger.debug('not enough weights on the hyperplane')
        return False

    # check if affine hull has the correct dimension
    dim = dim_affine_hull(weights)
    assert dim < R.dimension_affine_hull_weights, 'Absurd'
    logger.debug('dimension of affine hull is %d, supposed to be %d', dim,
                 R.dimension_affine_hull_weights - 1)
    return dim == R.dimension_affine_hull_weights - 1


def c_candidates(R, H):
    """Return possible :math:`c` such that :math:`(H,c)` is admissible for the given representation (see :func:`is_admissible`).

    :param R: the representation.
    :param H: the normal vector.
    :type R: :class:`Representation`
    :type H: :class:`sage.vector`
    :rtype: set of integers
    """
    # sort weights by H-weight (= possible values of c)
    H = vector(H)
    weights_by_c = defaultdict(list)
    for omega in R.weights:
        weights_by_c[omega.dot_product(H)].append(omega)

    # test all possible values of c
    cs = set()
    for c, weights in weights_by_c.iteritems():
        # impossible?
        if len(weights) < R.dimension_affine_hull_weights:
            continue

        # check if affine hull has correct dimension
        dim = dim_affine_hull(weights)
        assert dim < R.dimension_affine_hull_weights
        if dim == R.dimension_affine_hull_weights - 1:
            cs.add(c)

    return cs
