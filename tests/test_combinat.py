from __future__ import absolute_import, print_function
from collections import defaultdict
from sage.all import prod, factorial, QQ, vector, Permutation, Permutations
from moment_polytopes import *
import pytest


def test_rect_tableaux_22():
    tableaux = list(rect_tableaux(2, 2))
    assert len(tableaux) == 2
    assert [(1, 2), (3, 4)] in tableaux
    assert [(1, 3), (2, 4)] in tableaux


@pytest.mark.parametrize("a, b, count", [
    (2, 2, 2),
    (2, 3, 5),
    (3, 2, 5),
    (2, 4, 14),
    (4, 2, 14),
    (3, 3, 42),
    (4, 3, 462),
    (3, 4, 462),
    (4, 4, 24024),
])
def test_rect_tableaux_count(a, b, count):
    assert len(rect_tableaux(a, b)) == count


@pytest.mark.parametrize("a, b, count", [
    (2, 2, 2),
    (2, 3, 5),
    (3, 2, 5),
    (2, 4, 14),
    (4, 2, 14),
    (3, 3, 36),
    (3, 4, 295),
    (4, 3, 295),
    (4, 4, 6660),
])
def test_cubicle_tableaux_counts(a, b, count):
    assert len(cubicle_tableaux(a, b)) == count


@pytest.mark.parametrize("v, dominant", [
    ((3, 2, 4), False),
    ((4, 3, 2), True),
])
def test_is_dominant(v, dominant):
    assert is_dominant(v) == dominant


@pytest.mark.parametrize("dims, count", [
    ((2, 2), 3),
    ((2, 3), 6),
    ((3, 2), 6),
    ((2, 4), 11),
    ((4, 2), 11),
    ((3, 3), 17),
    ((3, 4), 56),
    ((4, 3), 56),
    ((4, 4), 457),
    ((2, 2, 3), 39),
])
def test_extremal_edges(dims, count):
    # count
    edges = extremal_edges(dims, include_perms=True)
    assert len(edges) == count

    # check that they are indeed extremal edges
    for edge in edges:
        assert is_extremal_edge(dims, edge)


def test_extremal_edges_implementations():
    for dims in [(2, 2), (2, 3), (3, 2), (3, 3)]:
        assert set(
            map(tuple, extremal_edges(dims, algorithm='generic'))) == set(
                map(tuple, extremal_edges(dims, algorithm='bipartite')))


@pytest.mark.parametrize(
    "dims, count",
    [
        ((2, 2), 2),
        ((2, 3), 6),
        ((3, 2), 6),  # same as above
        ((2, 4), 11),
        ((4, 2), 11),  # same as above
        ((3, 3), 10),
        ((3, 4), 56),
        ((4, 3), 56),  # same as above
        ((4, 4), 233),
        ((2, 2, 3), 25),
        ((2, 2, 2), 4),
        ((2, 2, 2, 2), 12),
    ])
def test_extremal_edges_up_to_perms(dims, count):
    # count
    edges = extremal_edges(dims, include_perms=False)
    assert len(edges) == count

    # check that they are indeed extremal edges
    for edge in edges:
        assert is_extremal_edge(dims, edge)


def test_is_extremal_edge():
    # non-example of Klyachko
    dims = (2, 2, 3)
    H = (2, -2, 1, -1, QQ('5/3'), QQ('2/3'), QQ('-7/3'))
    assert not is_extremal_edge(dims, H)


def test_is_extremal_edge_ieq():
    dims = (2, 2, 3, 12)
    ieq = vector(
        [0, 0, 0, 0, -2, 1, 1, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1]), 0

    # the first two asserts will fail since ieq is neither dominant nor primitive
    with pytest.raises(AssertionError):
        assert not is_extremal_edge_ieq(dims, ieq)
    with pytest.raises(AssertionError):
        assert not is_extremal_edge_ieq(dims, ieq, assert_dominant=False)
    assert is_extremal_edge_ieq(
        dims, ieq, assert_dominant=False, assert_primitive=False)


@pytest.mark.parametrize("a, b", [
    (2, 2),
    (2, 3),
    (3, 2),
    (2, 4),
    (4, 2),
    (3, 3),
    (3, 4),
    (4, 3),
    (4, 4),
])
def test_primitivity_lemma(a, b):
    """Test Lemma 6.3 in Vergne and Walter (2014)."""
    edges = extremal_edges((a, b), include_perms=False)

    root_system_A = ["A", a - 1]
    root_system_B = ["A", b - 1]
    root_system_AB = [root_system_A, root_system_B]

    for H in edges:
        # edge (H_A, H_B) should be primitive...
        assert is_dual_root_primitive(root_system_AB, H)

        # ...and by our lemma this should imply primitivity of H_A and H_B (if non-zero)
        H_A, H_B = H[:a], H[a:]
        assert H_A.is_zero() or is_dual_root_primitive(root_system_A, H_A)
        assert H_B.is_zero() or is_dual_root_primitive(root_system_B, H_B)


def P(*pi):
    return Permutation(list(pi))


@pytest.mark.parametrize("n, length, perms", [
    (3, 0, {P(1, 2, 3)}),
    (3, 1, {P(2, 1, 3), P(1, 3, 2)}),
    (3, 2, {P(2, 3, 1), P(3, 1, 2)}),
    (3, 3, {P(3, 2, 1)}),
    (4, 0, {P(1, 2, 3, 4)}),
    (4, 1, {P(2, 1, 3, 4), P(1, 3, 2, 4), P(1, 2, 4, 3)}),
    (4, 2, {
        P(2, 3, 1, 4), P(2, 1, 4, 3), P(1, 3, 4, 2), P(3, 1, 2, 4),
        P(1, 4, 2, 3)
    }),
    (4, 3, {
        P(4, 1, 2, 3), P(3, 1, 4, 2), P(1, 4, 3, 2), P(2, 3, 4, 1),
        P(3, 2, 1, 4), P(2, 4, 1, 3)
    }),
    (4, 4, {
        P(3, 2, 4, 1), P(4, 1, 3, 2), P(4, 2, 1, 3), P(3, 4, 1, 2),
        P(2, 4, 3, 1)
    }),
    (4, 5, {P(3, 4, 2, 1), P(4, 2, 3, 1), P(4, 3, 1, 2)}),
    (4, 6, {P(4, 3, 2, 1)}),
])
def test_perms_of_length(n, length, perms):
    assert set(perms_of_length(n, length=length)) == perms


@pytest.mark.parametrize("dims, total, tuples", [
    ((2, 2, 2), -1, set()),
    ((2, 2, 2), 0, {(0, 0, 0)}),
    ((2, 2, 2), 1, {(1, 0, 0), (0, 1, 0), (0, 0, 1)}),
    ((2, 2, 2), 2, {(1, 1, 0), (1, 0, 1), (0, 1, 1)}),
    ((2, 2, 2), 3, {(1, 1, 1)}),
    ((2, 2, 2), 4, set()),
    ((2, 2, 2, 2), -1, set()),
    ((2, 2, 2, 2), 0, {(0, 0, 0, 0)}),
    ((2, 2, 2, 2), 1,
     {(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)}),
    ((2, 2, 2, 2), 2, {(1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1), (0, 1, 1, 0),
                       (0, 1, 0, 1), (0, 0, 1, 1)}),
    ((2, 2, 2, 2), 3,
     {(1, 1, 1, 0), (1, 1, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1)}),
    ((2, 2, 2, 2), 4, {(1, 1, 1, 1)}),
    ((2, 2, 2, 2), 5, set()),
])
def test_length_tuples(dims, total, tuples):
    assert set(length_tuples(dims, total=total)) == tuples


@pytest.mark.parametrize(
    "v, n, expected_shuffles, expected_antishuffles",
    [
        # regular element
        ([3, 2, 1], 0, {P(1, 2, 3)}, {P(1, 2, 3)}),
        ([3, 2, 1], 1, {P(2, 1, 3), P(1, 3, 2)}, {P(2, 1, 3), P(1, 3, 2)}),
        ([3, 2, 1], 2, {P(2, 3, 1), P(3, 1, 2)}, {P(2, 3, 1), P(3, 1, 2)}),
        ([3, 2, 1], 3, {P(3, 2, 1)}, {P(3, 2, 1)}),
        # one degeneracy
        ([2, 2, 1], 0, {P(1, 2, 3)}, set()),
        ([2, 2, 1], 1, {P(1, 3, 2)}, {P(2, 1, 3)}),
        ([2, 2, 1], 2, {P(2, 3, 1)}, {P(3, 1, 2)}),
        ([2, 2, 1], 3, set(), {P(3, 2, 1)}),
        ([2, 1, 1], 0, {P(1, 2, 3)}, set()),
        ([2, 1, 1], 1, {P(2, 1, 3)}, {P(1, 3, 2)}),
        ([2, 1, 1], 2, {P(3, 1, 2)}, {P(2, 3, 1)}),
        ([2, 1, 1], 3, set(), {P(3, 2, 1)}),
        # completely degenerate
        ([1, 1, 1], 0, {P(1, 2, 3)}, set()),
        ([1, 1, 1], 1, set(), set()),
        ([1, 1, 1], 2, set(), set()),
        ([1, 1, 1], 3, set(), {P(3, 2, 1)}),
    ])
def test_shuffles_S3(v, n, expected_shuffles, expected_antishuffles):
    assert set(shuffles(v, length=n)) == expected_shuffles
    assert set(antishuffles(v, antilength=3 - n)) == expected_antishuffles


def _card_coset(v):
    blocks = defaultdict(int)
    for x in v:
        blocks[x] += 1

    return factorial(len(v)) / prod([factorial(l) for l in blocks.values()])


@pytest.mark.parametrize("v", [
    (4, 3, 2, 1),
    (4, 4, 2, 1),
    (4, 4, 2, 2),
    (4, 3, 3, 1),
    (4, 3, 2, 2),
    (4, 4, 4, 2),
    (4, 3, 3, 3),
    (4, 4, 4, 4),
])
def test_shuffles_S4(v):
    d = len(v)
    len_max = d * (d - 1) // 2

    # test shuffles
    num_shuffles = 0
    for l in range(len_max + 1):
        pis = shuffles(v, length=l)
        num_shuffles += len(pis)

        # check that permutation is shuffle of the correct length
        assert all(is_shuffle(pi, v) for pi in pis)
        assert all(pi.number_of_inversions() == l for pi in pis)

    # test shuffles
    num_antishuffles = 0
    for l in range(len_max + 1):
        pis = antishuffles(v, antilength=l)
        num_antishuffles += len(pis)

        # check that permutation is shuffle of the correct length
        assert all(is_antishuffle(pi, v) for pi in pis)
        assert all(pi.number_of_noninversions(2) == l for pi in pis)

    # check that the number of shuffles is the expected one
    card_permutations = factorial(d)
    assert num_shuffles == num_antishuffles == _card_coset(v)


def test_perm_action():
    def perm0(pi):
        return [i - 1 for i in pi]

    H = ['X', 'Y', 'Z', 'W']
    for pi in Permutations(4):
        for tau in Permutations(4):
            # 1-based action
            assert perm_action(pi, perm_action(tau, H)) == perm_action(
                pi.left_action_product(tau), H)

            # 0-based action
            pi0 = perm0(pi)
            tau0 = perm0(tau)
            pitau0 = perm0(pi.left_action_product(tau))
            assert perm_action(
                pi0, perm_action(tau0, H, zero_based=True),
                zero_based=True) == perm_action(
                    pitau0, H, zero_based=True)


def test_stabilizer_group():
    # all dimensions distinct
    stab = StabilizerGroup([2, 3, 4])
    assert stab.blocks == [[0], [1], [2]]
    assert stab.normal_form("BCA") == tuple("BCA")
    assert stab.orbit({"BCA"}) == {tuple("BCA")}
    assert stab.orbit({"XYX"}) == {tuple("XYX")}
    assert stab.orbit({"ZZZ"}) == {tuple("ZZZ")}

    # first and last dimension equal
    stab = StabilizerGroup([2, 3, 2])
    assert stab.blocks == [[0, 2], [1]]
    assert stab.normal_form("BCA") == tuple("ACB")
    assert stab.orbit({"BCA"}) == {tuple("ACB"), tuple("BCA")}
    assert stab.orbit({"XYX"}) == {tuple("XYX")}
    assert stab.orbit({"ZZZ"}) == {tuple("ZZZ")}

    # all dimension equal
    stab = StabilizerGroup([3, 3, 3])
    assert stab.blocks == [[0, 1, 2]]
    assert stab.normal_form("BCA") == tuple("ABC")
    assert stab.orbit({"BCA"}) == {
        tuple("ABC"), tuple("ACB"), tuple("BAC"), tuple("BCA"), tuple("CAB"),
        tuple("CBA")
    }
    assert stab.orbit({"XYX"}) == {tuple("YXX"), tuple("XYX"), tuple("XXY")}
    assert stab.orbit({"ZZZ"}) == {tuple("ZZZ")}
