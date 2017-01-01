from __future__ import absolute_import, print_function
import pytest
from sage.all import QQ, gcd, prod, vector
from moment_polytopes import *


def test_three_qubits(algorithm):
    # extremal edges (up to permutations)
    candidates = qmp.H_AB_dominant(2, 2, include_perms=False)
    assert candidates == {
        ((0, 0), (QQ("1/2"), -QQ("1/2"))),
        ((QQ("1/2"), -QQ("1/2")), (QQ("1/2"), -QQ("1/2"))),
    }
    assert len(qmp.H_AB_dominant(2, 2, include_perms=True)) == 2 + 1

    # candidates for dominant normal vectors
    candidates = sorted(qmp.H_ABC_dominant(2, 2, 2, include_perms=False))
    assert candidates == [
        ((0, 0), (0, 0), (QQ("1/2"), -QQ("1/2"))),
        ((0, 0), (QQ("1/2"), -QQ("1/2")), (QQ("1/2"), -QQ("1/2"))),
        ((QQ("1/2"), -QQ("1/2")), (QQ("1/2"), -QQ("1/2")),
         (QQ("1/2"), -QQ("1/2"))),
    ]
    assert len(qmp.H_ABC_dominant(2, 2, 2, include_perms=True)) == 3 + 3 + 1

    # candidates for corresponding z (easy to visualize for unit cube)
    R = external_tensor_product([2, 2, 2])
    assert c_candidates(R, sum(candidates[0], ())) == {QQ("1/2"), -QQ("1/2")}
    assert c_candidates(R, sum(candidates[1], ())) == {0}
    assert c_candidates(R, sum(candidates[2], ())) == {QQ("1/2"), -QQ("1/2")}

    # candidates for dominant admissible (H, z) (just the ones from above)
    candidates = qmp.H_dominant_admissible((2, 2, 2), include_perms=False)
    assert len(candidates) == 2 + 1 + 2

    # candidates for admissible (H, z)
    candidates = sorted(qmp.H_candidates((2, 2, 2), include_perms=False))
    assert candidates == [
        (((-QQ("1/2"), QQ("1/2")), (-QQ("1/2"), QQ("1/2")),
          (QQ("1/2"), -QQ("1/2"))), -QQ("1/2")),
        (((-QQ("1/2"), QQ("1/2")), (0, 0), (0, 0)), -QQ("1/2")),
        (((0, 0), (QQ("1/2"), -QQ("1/2")), (QQ("1/2"), -QQ("1/2"))), 0),
    ]

    # polytope
    hrepr = qmp.hrepr((2, 2, 2), algorithm=algorithm).irred()
    vs = hrepr.vertices()
    assert len(hrepr.ieqs) == 6
    assert sorted(map(tuple, vs)) == [
        (QQ("1/2"), QQ("1/2"), QQ("1/2"), QQ("1/2"), QQ("1/2"), QQ("1/2")),
        (QQ("1/2"), QQ("1/2"), QQ("1/2"), QQ("1/2"), 1, 0),
        (QQ("1/2"), QQ("1/2"), 1, 0, QQ("1/2"), QQ("1/2")),
        (1, 0, QQ("1/2"), QQ("1/2"), QQ("1/2"), QQ("1/2")),
        (1, 0, 1, 0, 1, 0),
    ]


def test_H_AB_dominant_334_bug():
    # Jonathan's extremal edges were dominant but not primitive, hence we missed some combinations when "assembling" H_ABCs
    H_As_got = {H_A for (H_A, H_B) in qmp.H_AB_dominant(3, 4)}
    assert (1, 0, -1) in H_As_got

    # check that we get at least those that appear in the actual inequalities (make them dominant and primitive in order to compare)
    ieqs = third_party._klyachko_qmp_bare_ieqs((3, 3, 9))
    H_As_expected = {tuple(H[:3])
                     for (H, z) in ieqs} | {tuple(H[3:6])
                                            for (H, z) in ieqs}
    H_As_expected = {tuple(sorted(H_A, reverse=True)) for H_A in H_As_expected}
    H_As_expected = {
        tuple(dual_root_primitive("A2", H_A))
        for H_A in H_As_expected if H_A != (0, 0, 0)
    }
    for H in H_As_expected:
        assert H in H_As_got


@pytest.mark.parametrize("dims", [
    (2, 2),
    (2, 3),
    (2, 4),
    (3, 3),
    (3, 4),
    (4, 4),
])
def test_H_AB_dominant_are_extremal_edges(dims):
    for (H_A, H_B) in qmp.H_AB_dominant(*dims):
        H = tuple(H_A) + tuple(H_B)
        assert is_extremal_edge(dims, H)


@pytest.mark.parametrize("a, b, c, count_wo_perms, count", [
    (2, 2, 2, 3, 7),
    (3, 3, 3, 17, 51),
    (4, 4, 4, 600, 3027),
])
def test_H_ABC_dominant_count(a, b, c, count_wo_perms, count):
    # excluding permutations
    candidates = qmp.H_ABC_dominant(a, b, c, include_perms=False)
    assert len(candidates) == count_wo_perms

    # include permutations
    candidates = qmp.H_ABC_dominant(a, b, c, include_perms=True)
    assert len(candidates) == count


@pytest.mark.parametrize("dims, count_wo_perms, normal_vector_count_wo_perms",
                         [
                             ((2, 2, 2), 5, 3),
                             ((3, 3, 3), 25, 17),
                             ((4, 4, 4), 484, 342),
                             ((2, 2, 4), 3, 3),
                             ((2, 2, 3, 12), 26, 26),
                         ])
def test_H_dominant_admissible(dims, count_wo_perms,
                               normal_vector_count_wo_perms):
    # candidates
    candidates = qmp.H_dominant_admissible(dims, include_perms=False)
    assert len(candidates) == count_wo_perms

    # normal vectors
    normal_vectors = {H for (H, z) in candidates}
    assert len(normal_vectors) == normal_vector_count_wo_perms


@pytest.mark.parametrize("dims, count_wo_perms", [
    ((2, 2, 2), 3),
    ((3, 3, 3), 41),
    ((4, 4, 4), 5633),
    ((2, 2, 4), 9),
    ((2, 2, 3, 12), 25206),
])
def test_H_candidates_count(dims, count_wo_perms):
    candidates = qmp.H_candidates(dims, include_perms=False)
    assert len(candidates) == count_wo_perms


@pytest.mark.parametrize("dims", [
    (2, 2, 2),
    (2, 2, 4),
    (3, 3, 3),
    (4, 4, 4),
    (2, 2, 3, 12),
])
def test_no_duplicates(dims):
    stab = StabilizerGroup(dims)

    # check that there are no duplicates in the set of dominant (H_A, H_B, H_C) candidates
    if len(dims) == 3:
        candidates = set()
        for hs in qmp.H_ABC_dominant(*dims, include_perms=False):
            hs_nf = stab.normal_form(hs)
            assert hs_nf not in candidates
            candidates.add(hs_nf)

    # check that there are no duplicates in the set of (H_A, H_B, H_C, z) candidates
    candidates = set()
    for (hs, z) in qmp.H_candidates(dims, include_perms=False):
        hs_nf = stab.normal_form(hs)
        assert (hs_nf, z) not in candidates
        candidates.add((hs_nf, z))


@pytest.mark.parametrize("dims", [
    (2, 2, 2),
    (3, 3, 3),
    (4, 4, 4),
    (2, 2, 4),
    (2, 3, 6),
    (2, 4, 8),
    (3, 3, 9),
    pytest.mark.slow((2, 2, 3, 12)),
    pytest.mark.slow((2, 2, 2, 2, 16)),
])
def test_H_candidates_traceless_primitive_extremal_edges(dims):
    candidates = qmp.H_candidates(dims, include_perms=False)

    # mixed-state situation?
    if dims[-1] == prod(dims[:-1]):
        # positivity constraint (in normal form)
        stab = StabilizerGroup(dims)
        hs_positivity = [(0, ) * d
                         for d in dims[:-1]] + [(-1, ) * (dims[-1] - 1) +
                                                (dims[-1] - 1, )]
        hs_positivity = stab.normal_form(hs_positivity)
        H_positivity = (hs_positivity, -1)
        assert H_positivity in candidates, 'Forgot positivity constraint?'

        # all other ones should be extremal edges
        for (hs, z) in candidates:
            if (hs, z) != H_positivity:
                ieq = (sum(hs, ()), z)
                assert is_extremal_edge_ieq(
                    dims,
                    ieq,
                    assert_dominant=False,
                    assert_primitive=True,
                    assert_traceless=True)
    else:
        for (hs, z) in candidates:
            # check traceless
            assert all(sum(h) == 0 for h in hs), 'Expect trace to be zero.'

            # check dual root primitive
            diffs = [x - y for h in hs for (x, y) in zip(h, h[1:])]
            c = gcd(diffs)
            assert c in [
                1, -1
            ], 'Expect vector that is primitive in dual of root lattice.'


@pytest.mark.parametrize("dims", [
    (2, 2, 2),
    (2, 2, 4),
    (2, 3, 6),
    (2, 4, 8),
    (3, 3, 3),
    pytest.mark.slow((3, 3, 9)),
    pytest.mark.slow((4, 4, 4)),
    pytest.mark.slow((2, 2, 3, 12)),
    pytest.mark.slow((2, 2, 2, 2, 16)),
])
def test_H_candidates_numerological(dims):
    # check that the "numerological condition" dim n_-(H < 0) == dim M(H < z) holds for all candidates (H, z)
    R = external_tensor_product(dims)
    for (hs, z) in qmp.H_candidates(dims, include_perms=False):
        H = vector(sum(hs, ()))
        dim_n_neg = sum(1 for alpha in R.negative_roots
                        if alpha.dot_product(H) < 0)
        dim_M_neg = sum(1 for omega in R.weights if omega.dot_product(H) < z)
        assert dim_n_neg == dim_M_neg


@pytest.mark.parametrize("dims, count_wo_perms", [
    ((2, 2, 2), 3),
    ((3, 3, 3), 25),
    ((4, 4, 4), 323),
    ((2, 2, 4), 9),
    pytest.mark.slow(((2, 2, 3, 12), 1330)),
    pytest.mark.slow(((2, 2, 2, 2, 16), 535)),
])
def test_H_ressayre(dims, count_wo_perms, algorithm):
    ressayre_wo_perms = qmp.H_ressayre(
        dims, algorithm=algorithm, include_perms=False)
    assert len(ressayre_wo_perms) == count_wo_perms


@pytest.mark.parametrize(
    "dims, num_facets, num_facets_wo_perms",
    [
        ((2, 2, 2), 6, 2),
        ((3, 3, 3), 45, 10),
        ((4, 4, 4), 270, 50),
        ((2, 2, 4), 13,
         9),  # Klyachko: 7 + (1 + 1 + 3 Weyl) + (1 positivity) = 13
        ((2, 3, 6), 50, 50),  # Klyachko:   41 + (1 + 2 + 5) + 1 = 50
        ((2, 4, 8), 246, 246),  # Klyachko:  234 + (1 + 3 + 7) + 1 = 246
        pytest.mark.slow(((3, 3, 9), 400,
                          208)),  # Klyachko:  387 + (2 + 2 + 8) + 1 = 400
        pytest.mark.slow(
            ((2, 2, 3, 12), 599, 322)
        ),  # Klyachko:  442 + (1 + 1 + 2 + 11) + 1 = 457 -- BUT HIS INEQUALITIES ARE NOT CORRECT AS PRINTED
        pytest.mark.slow(
            ((2, 2, 2, 2, 16), 825,
             67)),  # Klyachko:  805 + (1 + 1 + 1 + 1 + 15) + 1 = 825
    ])
def test_qmp_hrepr(dims, num_facets, num_facets_wo_perms, algorithm):
    hrepr = qmp.hrepr(dims, algorithm=algorithm)
    assert len(hrepr.ieqs) == num_facets
    assert len(qmp.ieqs_wo_perms(dims, hrepr.ieqs)) == num_facets_wo_perms


@pytest.mark.parametrize(
    "dims, num_vertices, num_vertices_wo_perms",
    [
        ((2, 2, 2), 5, 3),
        ((3, 3, 3), 33, 11),
        ((4, 4, 4), 328, 65),
        ((2, 2, 4), 10, 8),
        ((2, 3, 6), 56, 56),
        pytest.mark.slow(((2, 4, 8), 248, 248)),  # lrs takes too long?!
        pytest.mark.slow(((3, 3, 9), 561, 297)),  # lrs takes too long?!
        pytest.mark.slow(((2, 2, 3, 12), 0, 0)),  # lrs takes too long?!
        pytest.mark.slow(((2, 2, 2, 2, 16), 0, 0)),  # lrs takes too long?!
    ])
def test_qmp_vrepr(dims, num_vertices, num_vertices_wo_perms, algorithm):
    vrepr = qmp.vrepr(dims, algorithm=algorithm)
    assert not vrepr.rays and not vrepr.lines
    assert len(vrepr.vertices) == num_vertices
    assert len(
        qmp.vertices_wo_perms(dims, vrepr.vertices)) == num_vertices_wo_perms


@pytest.mark.parametrize("dims, expected_fn", [
    ((2, 2, 2), third_party.higuchi_hrepr),
    ((2, 2, 4), third_party.bravyi_hrepr),
    ((3, 3, 3), third_party.franz_hrepr),
])
def test_basic_scenarios(dims, expected_fn, algorithm):
    got = qmp.hrepr(dims, algorithm=algorithm)
    expected = expected_fn()
    assert got.vrepr() == expected.vrepr()


@pytest.mark.parametrize("dims", [
    (2, 2, 2),
    (3, 3, 3),
    (2, 2, 4),
    (2, 3, 6),
    (2, 4, 8),
    pytest.mark.slow((3, 3, 9)),
    pytest.mark.slow(pytest.mark.xfail((2, 2, 3, 12))),
    pytest.mark.slow((2, 2, 2, 2, 16)),
])
def test_vs_klyachko(dims, algorithm):
    # fetch irredundant H-representations
    us = qmp.hrepr(dims, algorithm=algorithm)
    klyachko = third_party.klyachko_qmp_hrepr(dims)

    # compare their normal forms
    def ieqs_nf(l):
        m = set()
        for ieq in l:
            H, c = qmp.facet_normal_form(dims, ieq)
            m.add((tuple(H), c))
        return m

    us = ieqs_nf(us.ieqs)
    klyachko = ieqs_nf(klyachko.ieqs)
    assert len(us) == len(klyachko)
    assert us == klyachko


@pytest.mark.parametrize(
    "dims",
    [
        (2, 2, 2),
        (2, 2, 3),
        (2, 2, 4),
        (2, 3, 6),
        (3, 3, 3),
        (4, 4, 4),
        # it takes a lot of time to verify some of these determinants using Sage:
        pytest.mark.slow((2, 2, 3, 12)),
        pytest.mark.slow((2, 4, 8)),
        pytest.mark.slow((3, 3, 9)),
        pytest.mark.slow((2, 2, 2, 2, 16)),
    ])
def test_compare_against_sage(dims, algorithm):
    if algorithm and algorithm != 'sage':
        ieqs = qmp.H_ressayre(dims, algorithm=algorithm, include_perms=False)
        ieqs_ref = qmp.H_ressayre(dims, include_perms=False)
        assert ieqs == ieqs_ref
