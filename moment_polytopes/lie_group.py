from __future__ import absolute_import, print_function
from sage.all import Integer, vector, gcd, ZZ, QQ, RootSystem, Partition, SemistandardTableaux, matrix, copy, Tableau, cartesian_product, GelfandTsetlinPatterns, GelfandTsetlinPattern, prod, factorial, diagonal_matrix
from . import HRepr
from .utils import dim_affine_hull

__all__ = [
    'Representation',
    'weyl_module',
    'external_tensor_product',
    'is_dual_root_primitive',
    'dual_root_primitive',
    'positive_weyl_chamber_hrepr',
]


class Representation(object):
    """Base class for representations of complex reductive Lie groups.

    See :func:`weyl_module` and :func:`external_tensor_product` for concrete representations.

    .. attribute:: root_system

       Root system of the corresponding Lie group.

    .. attribute:: ambient_dim

       Dimension of the ambient space containing all negative roots and weights.

    .. attribute:: negative_roots

       List of negative roots of the corresponding Lie group.

    .. attribute:: weights

       List of weights of the corresponding Lie group (repeated according to multiplicity).

    .. attribute:: reduced_eqns

       Equations :math:`(H,c)` defining the affine subspace spanned by the weights.
    """

    def __init__(self):
        self._dimension_affine_hull_weights = None

    @property
    def dimension(self):
        """The dimension of the representation."""
        return len(self.weights)

    @property
    def positive_roots(self):
        """List of positive roots of the corresponding Lie group."""
        return [-alpha for alpha in self.negative_roots]

    @property
    def dimension_affine_hull_weights(self):
        """The dimension of the affine subspace spanned by the weights."""
        # compute once and cache for future invocations
        if self._dimension_affine_hull_weights is None:
            self._dimension_affine_hull_weights = dim_affine_hull(self.weights)
        return self._dimension_affine_hull_weights

    @property
    def reduced_positive_weyl_chamber_hrepr(self):
        """Return intersection of positive Weyl chamber and the affine subspace spanned by the weights.

        :type: :class:`HRepr`
        """
        pwc = positive_weyl_chamber_hrepr(self.root_system)
        eqns = HRepr(eqns=self.reduced_eqns, ambient_dim=self.ambient_dim)
        return pwc & eqns

    def negative_root_action(self, idx_negative_root, idx_weight_vector=None):
        r"""Apply lowering operator corresponding to a negative root.

        The indices refer to :attr:`negative_roots` and :attr:`weights`, respectively.

        :param idx_negative_root: index of the negative root.
        :param idx_weight_vector: index of the weight vector, or ``None``. If ``None`` then the matrix representation of the lowering operator is returned.
        :rtype: :class:`sage.vector` or :class:`sage.matrix`
        """
        raise NotImplementedError

    def weight_vector(self, idx_weight_vector):
        """Return weight vector for given index (this is just a standard basis vector with a single ``1`` at the given index).

        :param idx_weight_vector: index of the weight vector (indexes the :attr:`weights` list).
        :rtype: :class:`sage.vector`
        """
        return vector(ZZ, self.dimension, {idx_weight_vector: 1})


class WeylModule(Representation):
    """Return rational irreducible representation of :math:`GL(d)` with given highest weight.

    >>> weyl_module(4, [2, 1])
    WeylModule(4, (2, 1, 0, 0))

    We work with an orthogonal basis labeled by Gelfand-Tsetlin pattern that has the properties that
    1. The highest weight is a unit vector.
    2. Actions of positive and negative root operators are adjoint of each other.
    3. The matrix elements of the positive and negative root operators (as well as of the torus generators) are integral.

    :param d: the rank of the Lie group :math:`GL(d)`.
    :param highest_weight: the highest weight. Should have no more than ``d`` parts. Will be padded by 0 if has fewer than ``d`` parts.
    :type highest_weight: `list`
    :rtype: :class:`Representation`
    """

    def __init__(self, d, highest_weight):
        super(WeylModule, self).__init__()

        #: The rank of GL(d).
        self.d = d

        # normalize highest weight
        highest_weight = list(highest_weight)
        assert len(
            highest_weight) <= d, 'Highest weight has more than %s parts.' % d
        assert len(
            highest_weight
        ) == d or highest_weight[-1] >= 0, 'Cannot pad highest weight by zero.'
        highest_weight = highest_weight + [0] * (d - len(highest_weight))
        highest_weight = vector(highest_weight)

        #: The highest weight
        self.highest_weight = highest_weight

        #: The partition corresponding to the highest weight (tensored with a minimal number of determinants)
        shift = max(-self.highest_weight[d - 1], 0)
        self.partition = Partition(self.highest_weight + vector([shift] * d))

        #: The Gelfand-Tsetlin pattern labeling the basis vectors.
        self.patterns = GelfandTsetlinPatterns(top_row=highest_weight).list()
        self.patterns = sorted(self.patterns, reverse=True)

        #: The tableaux labeling the basis vectors (see `partition`)
        self.tableaux = []
        for p in self.patterns:
            q = [[x + shift for x in row] for row in p]
            t = GelfandTsetlinPattern(q).to_tableau()
            self.tableaux.append(t)

        # implement properties of base class
        self.ambient_dim = d
        self.root_system = RootSystem(['A', d - 1])
        ambient_space = self.root_system.ambient_space()
        self._simple_roots = [
            vector(ZZ, d, {
                i: 1,
                i + 1: -1
            }) for i in range(d - 1)
        ]
        self.negative_roots = map(vector, ambient_space.negative_roots())
        self.weights = [vector(p.weight()) for p in self.patterns]
        self.reduced_eqns = [(vector([1] * d), sum(self.highest_weight))]

        # make each pattern a list of lists (since GelfandTsetlinPattern do not compare well)
        self.patterns = [map(list, p) for p in self.patterns]

        # compute norm square of (unnormalized) Gelfand-Tsetlin basis vectors
        self._norm_squares = []
        for pattern in self.patterns:

            def l(k, j):
                return pattern[self.d - k][j - 1] - j + 1

            ns = prod(
                prod(
                    QQ((factorial(l(k, i) - l(k - 1, j)),
                        factorial(l(k - 1, i) - l(k - 1, j))))
                    for i in range(1, k) for j in range(i, k)) * prod(
                        QQ((factorial(l(k, i) - l(k, j) - 1),
                            factorial(l(k - 1, i) - l(k, j) - 1)))
                        for i in range(1, k + 1) for j in range(i + 1, k + 1))
                for k in range(2, self.d + 1))
            self._norm_squares.append(ns)

        # precompute action of negative simple roots in the (unnormalized) Gelfand-Tsetlin basis
        # (see https://arxiv.org/pdf/math/0211289.pdf, Theorem 2.3, (2.7))
        self._negative_simple_root_actions = []
        for k in range(1, len(self._simple_roots) + 1):
            # E_{k+1,k}
            d = {}
            for col, pattern in enumerate(self.patterns):

                def l(k, j):
                    return pattern[self.d - k][j - 1] - j + 1

                for i in range(1, k + 1):
                    # build new pattern (clone the list of lists)
                    new_pattern = list(map(list, pattern))
                    new_pattern[self.d - k][i - 1] -= 1
                    try:
                        GelfandTsetlinPattern(new_pattern)
                    except AssertionError:
                        continue
                    row = self.patterns.index(new_pattern)

                    # compute prefactor
                    numer = prod([l(k, i) - l(k - 1, j) for j in range(1, k)])
                    denom = prod(
                        [l(k, i) - l(k, j) for j in range(1, k + 1) if j != i])
                    assert (row, col) not in d
                    d[(row, col)] = QQ(numer) / QQ(denom)

            m = matrix(self.dimension, self.dimension, d, sparse=True)
            self._negative_simple_root_actions.append(m)

    def __repr__(self):
        return 'WeylModule(%d, %r)' % (self.d, self.highest_weight)

    def negative_root_action(self, idx_negative_root, idx_weight_vector=None):
        # write alpha as sum of negative simple roots such that all partial sums are again negative roots
        decomposition = []
        alpha = self.negative_roots[idx_negative_root]
        while alpha:
            for i, beta in enumerate(self._simple_roots):
                if beta.dot_product(alpha) < 0:
                    decomposition = [i] + decomposition
                    alpha = alpha + beta
        assert len(decomposition) > 0

        # as a consequence, X_alpha is equal to the iterated commutator [[X_1, X_2], ...] where X_i correspond to the simple roots alpha_i in the decomposition
        def X(alpha):
            return matrix(
                self.d,
                self.d, {(list(alpha).index(1), list(alpha).index(-1)): 1},
                sparse=True)

        m = X(-self._simple_roots[decomposition[0]])
        for i in decomposition[1:]:
            m = m.commutator(X(-self._simple_roots[i]))
        assert m == X(self.negative_roots[idx_negative_root])

        # compute iterated commutator
        m = self._negative_simple_root_actions[decomposition[0]]
        for i in decomposition[1:]:
            m = m.commutator(self._negative_simple_root_actions[i])

        # return result
        if idx_weight_vector is not None:
            return m.column(idx_weight_vector)
        else:
            return m

    def positive_root_action(self, idx_positive_root, idx_weight_vector=None):
        # https://arxiv.org/pdf/math/0211289.pdf, p. 10
        m = self.negative_root_action(idx_positive_root)

        # not an orthonormal basis
        D = diagonal_matrix(self._norm_squares)
        D_inv = diagonal_matrix([1 / ns for ns in self._norm_squares])
        m = D_inv * m.transpose() * D

        if idx_weight_vector is not None:
            return m.column(idx_weight_vector)
        else:
            return m

    def torus_action(self, idx_positive_root):
        alpha = self.positive_roots[idx_positive_root]
        return diagonal_matrix(map(alpha.dot_product, self.weights))

    def tableau_norm_squared(self, tableau):
        """Return norm squared of weight vector for given tableau.

        :param tableau: the semistandard tableau labeling the basis vector.
        :type tableau: :class:`sage.Tableau`
        :rtype: `int`
        """
        idx = self.tableaux.index(Tableau(tableau))
        return self._norm_squares[idx]

    def tableau_vector(self, tableau):
        """Return basis vector for given tableau.

        :param tableau: the semistandard tableau labeling the basis vector.
        :type tableau: :class:`sage.Tableau`
        :rtype: :class:`sage.vector`
        """
        i = self.tableaux.index(Tableau(tableau))
        return vector(ZZ, self.dimension, {i: 1})


weyl_module = WeylModule


def _embed_vector(v, k, dims):
    """Inject vector into k-th summand of direct sum."""
    before = sum(dims[:k])
    after = sum(dims[k + 1:])
    return vector([0] * before + list(v) + [0] * after)


class ExternalTensorProduct(Representation):
    """Construct external tensor product of given representations.

    >>> external_tensor_product([2,3])
    ExternalTensorProduct([WeylModule(2, (1, 0)), WeylModule(3, (1, 0, 0))])

    :param Rs: the representations. Integers :math:`d` are interpreted as fundamental representations of :math:`GL(d)`.
    :type Rs: list of :class:`Representation` or integers.
    :rtype: :class:`Representation`
    """

    def __init__(self, Rs):
        super(ExternalTensorProduct, self).__init__()

        #: The tensor factor representations.
        self.factors = Rs = [weyl_module(R, [1]) if R in ZZ else R for R in Rs]

        # implement properties of base class
        self.root_system = RootSystem([R.root_system for R in Rs])

        ambient_dims = [R.ambient_dim for R in Rs]
        self.ambient_dim = sum(ambient_dims)

        # collect negative roots and build reverse look-up table
        self.negative_roots = []
        self._negative_root_table = []
        for k, R in enumerate(Rs):
            for i, alpha in enumerate(R.negative_roots):
                self.negative_roots.append(
                    _embed_vector(alpha, k, ambient_dims))
                self._negative_root_table.append((k, i))

        # collect weights and build reverse look-up table
        self.weights = [
            vector(sum(map(list, ws), []))
            for ws in cartesian_product([R.weights for R in Rs])
        ]
        self._weight_table = [
            tuple(js)
            for js in cartesian_product([range(len(R.weights)) for R in Rs])
        ]

        # equations satisfied by weights
        self.reduced_eqns = []
        for k, R in enumerate(Rs):
            self.reduced_eqns += [(_embed_vector(H, k, ambient_dims), c)
                                  for (H, c) in R.reduced_eqns]

    def __repr__(self):
        return 'ExternalTensorProduct(%r)' % self.factors

    def negative_root_action(self, idx_negative_root, idx_weight_vector=None):
        # entire matrix requested? (XXX: slow - use tensor product instead of this hack)
        if idx_weight_vector is None:
            columns = [
                self.negative_root_action(idx_negative_root, i)
                for i in range(self.dimension)
            ]
            return matrix(columns).transpose()

        # look-up indices
        k, i = self._negative_root_table[idx_negative_root]
        js = self._weight_table[idx_weight_vector]

        # apply action to tensor factor
        v = self.factors[k].negative_root_action(i, js[k])

        # build result
        w = vector(QQ, self.dimension)
        for jj in range(self.factors[k].dimension):
            idx = self._weight_table.index(js[:k] + (jj, ) + js[k + 1:])
            w[idx] += v[jj]
        return w


external_tensor_product = ExternalTensorProduct


def positive_weyl_chamber_hrepr(root_system):
    """Return H-representation of positive Weyl chamber for given root system.

    :param root_system: the root system.
    :type root_system: :class:`sage.RootSystem`
    :rtype: :class:`HRepr`
    """
    ambient_space = RootSystem(root_system).ambient_space()
    simple_roots = map(vector, ambient_space.simple_roots())
    return HRepr(
        ieqs=[(alpha, 0) for alpha in simple_roots],
        ambient_dim=ambient_space.dimension())


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
