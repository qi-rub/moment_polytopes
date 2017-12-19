from __future__ import absolute_import, print_function
import logging, pytest
from sage.all import QQ, vector
from moment_polytopes import *

logger = logging.getLogger(__name__)


def test_fermi_3_6():
    # check bare inequalities
    bare_hrepr = third_party.klyachko_fermi_hrepr(3, 6, bare=True)
    bare_hrepr_expected = HRepr(
        eqns=[
            ((-1, 0, 0, 0, 0, -1), -1),
            ((0, -1, 0, 0, -1, 0), -1),
            ((0, 0, -1, -1, 0, 0), -1),
        ],
        ieqs=[
            ((0, 0, 0, -1, 1, 1), 0),
        ])
    assert bare_hrepr == bare_hrepr_expected

    # check complete inequalities
    hrepr = third_party.klyachko_fermi_hrepr(3, 6).irred()
    assert len(
        hrepr.eqns
    ) == 3  # three Borland--Dennis equations (that together imply the trace equation!)
    assert len(
        hrepr.ieqs
    ) == 1 + 3  # Borland-Dennis facet plus three Weyl chamber inequalities

    # ...and vertices
    vertices = {tuple(v[:3]) for v in hrepr.vertices()}
    vertices_expected = {
        ("1", "1", "1"),
        ("3/4", "3/4", "1/2"),
        ("1/2", "1/2", "1/2"),
        ("1", "1/2", "1/2"),
    }
    vertices_expected = {tuple(QQ(x) for x in v) for v in vertices_expected}
    assert vertices == vertices_expected


@pytest.mark.parametrize("n, d", third_party.KLYACHKO_FERMI_SCENARIOS)
def test_fermi_ieqs(algorithm, n, d):
    # retrieve inequalities according to Klyachko (w/o Weyl chamber)
    hrepr = third_party.klyachko_fermi_hrepr(n, d, bare=True)

    # verify all inequalities
    R = weyl_module(d, [1] * n)
    T = ressayre_tester(R, algorithm=algorithm)
    for i, ieq in enumerate(hrepr.ieqs):
        logger.debug("testing inequality %d/%d...", i + 1, len(hrepr.ieqs))
        assert T.is_ressayre(ieq)


def test_qmp_333():
    hrepr = third_party.klyachko_qmp_hrepr([3, 3, 3])
    assert len(hrepr.ieqs) == 45
    assert len(hrepr.vertices()) == 33


@pytest.mark.parametrize(
    "dims",
    [
        (3, 2, 6),
        pytest.param((4, 2, 8), marks=pytest.mark.slow),
        pytest.param((3, 3, 9), marks=pytest.mark.slow),
        # there is a mistake in Klyachko's table:
        pytest.param(
            (2, 2, 3, 12), marks=[pytest.mark.slow, pytest.mark.xfail]),
        pytest.param((2, 2, 2, 2, 16), marks=pytest.mark.slow),
    ])
def test_klyachko_qmp_hrepr_bare(algorithm, dims):
    # retrieve inequalities according to Klyachko (w/o Weyl chamber and permutations)
    bare_hrepr = third_party.klyachko_qmp_hrepr(dims, bare=True, irred=False)

    # verify inequalities
    R = external_tensor_product(dims)
    T = ressayre_tester(R, algorithm=algorithm)
    for i, ieq in enumerate(bare_hrepr.ieqs):
        logger.debug("testing inequality %d/%d...", i + 1, len(
            bare_hrepr.ieqs))
        assert T.is_ressayre(ieq)


@pytest.mark.parametrize("dims", [
    (3, 2, 6),
    pytest.param((4, 2, 8), marks=pytest.mark.slow),
    pytest.param((3, 3, 9), marks=pytest.mark.slow),
    pytest.param((2, 2, 3, 12), marks=pytest.mark.slow),
    pytest.param((2, 2, 2, 2, 16), marks=pytest.mark.slow),
])
def test_klyachko_qmp_hrepr_positivity_bug(dims):
    # retrieve QMP polytope according to Klyachko
    hrepr = third_party.klyachko_qmp_hrepr(dims, irred=False)

    # check that lambda_{AB,ab} >= 0 holds for all vertices (and hence for the polytope) -- there was a bug where we forgot to add this (which also implied that the polytope was unbounded!)
    vertices = hrepr.vertices()
    assert all(v[-1] >= 0 for v in vertices)


def test_klyachko_vs_higuchi():
    klyachko3 = third_party.klyachko_qmp_hrepr([2, 2, 2]).vrepr()
    higuchi3 = third_party.higuchi_hrepr(num_qubits=3).vrepr()
    assert klyachko3 == higuchi3

    klyachko5 = third_party.klyachko_qmp_hrepr([2, 2, 2, 2, 2]).vrepr()
    higuchi5 = third_party.higuchi_hrepr(num_qubits=5).vrepr()
    assert klyachko5 == higuchi5


def test_klyachko_vs_franz():
    klyachko = third_party.klyachko_qmp_hrepr([3, 3, 3]).vrepr()
    franz = third_party.franz_hrepr().vrepr()
    assert klyachko == franz


def test_klyachko_vs_bravyi():
    klyachko = third_party.klyachko_qmp_hrepr([2, 2, 4]).vrepr()
    bravyi = third_party.bravyi_hrepr().vrepr()
    assert klyachko == bravyi


def test_klyachko_qmp_22312_wrong():
    # load the inequalities for 3x2x6 and slice with 2x2x3 to get the Bravyi inequalities for mixed states of rank <= 3
    def proj((H, c)):
        return list(H[:2]) + list(H[3:8]), c

    bravyi_correct = third_party.klyachko_qmp_hrepr([3, 2,
                                                     6]).map(proj).irred()

    # attempt to get the same by loading the inequalities for 2x2x3x12 and slicing with ... x {1}
    dims = (2, 2, 3, 12)
    ieqs = third_party._klyachko_qmp_bare_ieqs(dims)
    ieqs += [(list(H[2:4]) + list(H[:2]) + list(H[4:]), c) for (H, c) in ieqs]
    hrepr = HRepr(ieqs=ieqs) & external_tensor_product(
        dims).reduced_positive_weyl_chamber_hrepr

    def proj((H, c)):
        return H[:2 + 2 + 3], c - H[7]

    #bravyi_wrong = klyachko_qmp_hrepr(dims).map(proj).irred()
    bravyi_wrong = hrepr.map(proj).irred()
    assert bravyi_correct.vertices() != bravyi_wrong.vertices()

    # slicing with lambda_C = (0.6, 0.3, 0.1) makes this very apparent!
    slice = HRepr(eqns=[
        ((0, 0, 0, 0, 1, 0, 0), QQ("6/10")),
        ((0, 0, 0, 0, 0, 1, 0), QQ("3/10")),
        ((0, 0, 0, 0, 0, 0, 1), QQ("1/10")),
    ])
    bravyi_slice_correct = (bravyi_correct & slice).irred()
    bravyi_slice_wrong = (bravyi_wrong & slice).irred()

    vertices_correct = {(v[0], v[2]) for v in bravyi_slice_correct.vertices()}
    vertices_wrong = {(v[0], v[2]) for v in bravyi_slice_wrong.vertices()}
    vertices_expected = {(QQ(str(x)), QQ(str(y)))
                         for (x, y) in [
                             ("8/10", "5/10"),
                             ("5/10", "5/10"),
                             ("5/10", "8/10"),
                             ("9/10", "6/10"),
                             ("9/10", "7/10"),
                             ("6/10", "9/10"),
                             ("7/10", "9/10"),
                         ]}
    assert vertices_correct == vertices_expected
    assert vertices_wrong != vertices_expected

    # let's make it painfully clear
    # - klyachko claims that his inequalities hold if lambda_{A,1} >= lambda_{B,1}
    # - this is certainly true for the following test spectrum
    wrong_spec = map(QQ, [
        "95/100", "5/100", "65/100", "35/100", "6/10", "3/10", "1/10", 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0
    ])
    wrong_spec_223 = wrong_spec[:2 + 2 + 3]

    # - check that this spectrum satisfies all of Klyachko's inequalities
    for (H, c) in third_party._klyachko_qmp_bare_ieqs(dims):
        assert H.dot_product(vector(wrong_spec)) >= c

    # - likewise does its truncation to a 2x2x3 spectrum
    assert wrong_spec_223 in bravyi_wrong
    assert wrong_spec_223 in bravyi_slice_wrong

    # - but it does *NOT* satisfy the inequalities by bravyi --- thus there has to be a mistake!
    assert wrong_spec_223 not in bravyi_correct
    assert wrong_spec_223 not in bravyi_slice_correct


@pytest.mark.slow
def test_klyachko_qmp_hrepr_22312_conjecture(algorithm):
    # my conjecture is that klyachko somehow got the extremal edge computation wrong
    dims = (2, 2, 3, 12)
    R = external_tensor_product(dims)
    T = ressayre_tester(R, algorithm=algorithm)

    ieqs = third_party._klyachko_qmp_bare_ieqs(dims)
    for i, ieq in enumerate(ieqs):
        logger.debug("testing inequality %d/%d...", i + 1, len(ieqs))\

        # check that admissibility is just the same as being an extremal edge (in the situation of Klyachko)
        ok = is_admissible(R, ieq)
        assert ok == is_extremal_edge_ieq(
            dims, ieq, assert_dominant=False, assert_primitive=False)

        # check that we can prove all inequalities where the edge is in fact extremal
        if ok:
            assert T.is_ressayre(ieq)
