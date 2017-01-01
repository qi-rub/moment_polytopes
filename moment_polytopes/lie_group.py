from __future__ import absolute_import, print_function
from sage.all import Integer, vector, gcd, QQ

__all__ = ["Group", "GL", "times", "CartesianProduct"]


class Group(object):
    """Compact, connected Lie group"""

    @property
    def pwc_hrepr(self):
        """Return H-representation of positive Weyl chamber."""
        if hasattr(self, 'simple_roots'):
            roots = self.simple_roots
        else:
            roots = [-alpha for alpha in self.negative_roots]
        return HRepr(
            ieqs=[(alpha, 0) for alpha in roots], ambient_dim=self.rank)

    def make_dual_root_primitive(self, H):
        """Make H primitive in the dual root lattice."""
        H = vector(H)
        c = gcd(alpha.dot_product(H) for alpha in self.simple_roots)
        assert c, 'Expected non-zero vector.'
        return vector([QQ((x, c)) for x in H])

    def is_dual_root_primitive(self, H):
        """Determine if H is primitive in the dual lattice of the root lattice."""
        H = vector(H)
        c = gcd(alpha.dot_product(H) for alpha in self.simple_roots)
        return c == 1 or c == -1


class GL(Group):
    """The general linear group GL(d)."""

    def __init__(self, d):
        super(GL, self).__init__()
        self.rank = d
        self.d = d

        # compute simple roots
        self.simple_roots = []
        for i in range(d - 1):
            root = vector(QQ, d, {i: 1, i + 1: -1}, sparse=False)
            self.simple_roots.append(root)

        # compute *lists* of positive and negative roots (order matters since we will later address them by index)
        self.negative_roots = []
        self.negative_root_ij = []
        for j in range(d):
            for i in range(j + 1, d):
                root = vector(QQ, d, {i: 1, j: -1}, sparse=False)
                self.negative_roots.append(root)
                self.negative_root_ij.append((i, j))

    def fundamental(self, n=1):
        """Return n-th fundamental (i.e., antisymmetric) representation."""
        assert 1 <= n <= self.rank
        return self.antisymmetric(n)

    def antisymmetric(self, n):
        """Return n-th antisymmetric power."""
        return AntisymmetricRepr(self, n)

    def weyl_module_21(self):
        """Return Wey module with Young diagram (2,1)."""
        return WeylModule21(self)


def _embed_vector(v, k, dims):
    before = sum(dims[:k])
    after = sum(dims[k + 1:])
    return vector([0] * before + list(v) + [0] * after)


def _embed_ieq((H, c), k, dims):
    return (_embed_vector(H, k, dims), c)


class CartesianProduct(Group):
    """Cartesian product of Lie groups."""

    def __init__(self, Gs):
        super(CartesianProduct, self).__init__()
        ranks = [G.rank for G in Gs]
        self.Gs = Gs
        self.rank = sum(ranks)

        # collect simple roots
        self.simple_roots = []
        for k, G in enumerate(Gs):
            self.simple_roots += [
                _embed_vector(alpha, k, ranks) for alpha in G.simple_roots
            ]

        # collect negative roots and build reverse look-up table
        self.negative_roots = []
        self.negative_root_table = []
        for k, G in enumerate(Gs):
            for i, alpha in enumerate(G.negative_roots):
                self.negative_roots.append(_embed_vector(alpha, k, ranks))
                self.negative_root_table.append((k, i))


def times(Gs):
    """Return Cartesian product of given Lie groups."""
    Gs = list(GL(G) if isinstance(G, (int, Integer)) else G for G in Gs)
    return CartesianProduct(Gs)
