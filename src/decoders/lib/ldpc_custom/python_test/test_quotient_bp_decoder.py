from itertools import product

import numpy as np
import pytest
import scipy.sparse

from ldpc_modified.quotient_bp_decoder import QuotientBpDecoder


def gf2_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (left @ right) % 2


def exact_log_partition(
    D: np.ndarray,
    offset: np.ndarray,
    probabilities: np.ndarray,
) -> float:
    total = 0.0
    for bits in product((0, 1), repeat=D.shape[1]):
        coordinate = np.asarray(bits, dtype=np.uint8)
        error = offset ^ gf2_product(D, coordinate).astype(np.uint8)
        total += float(
            np.prod(np.where(error == 1, probabilities, 1.0 - probabilities))
        )
    return float(np.log(total))


def test_cpp_quotient_model_satisfies_t_k_q_d_invariants():
    H = np.array([[1, 1, 1, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 1, 0]], dtype=np.uint8)
    decoder = QuotientBpDecoder(H, L, error_rate=0.1, max_iter=3)

    T = decoder.T.toarray()
    K = decoder.K.toarray()
    Q = decoder.Q.toarray()
    D = decoder.D.toarray()
    logical_map = decoder.logical_map.toarray()

    np.testing.assert_array_equal(gf2_product(gf2_product(H, T), H), H)
    np.testing.assert_array_equal(
        gf2_product(H, K), np.zeros((H.shape[0], K.shape[1]), dtype=np.uint8)
    )
    np.testing.assert_array_equal(
        gf2_product(H, D), np.zeros((H.shape[0], D.shape[1]), dtype=np.uint8)
    )
    np.testing.assert_array_equal(
        gf2_product(L, D), np.zeros((L.shape[0], D.shape[1]), dtype=np.uint8)
    )
    np.testing.assert_array_equal(
        gf2_product(H, Q), np.zeros((H.shape[0], Q.shape[1]), dtype=np.uint8)
    )
    np.testing.assert_array_equal(gf2_product(L, Q), logical_map)
    assert D.shape[1] + Q.shape[1] == K.shape[1]
    assert decoder.rank_h == 1
    assert decoder.logical_rank == 1
    assert decoder.degeneracy_dimension == 2


def test_h_bp_records_every_iteration_and_unique_cosets_in_first_visit_order():
    H = np.array(
        [
            [0, 0, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=np.uint8,
    )
    L = np.array([[1, 0, 1, 1, 0, 1, 0]], dtype=np.uint8)
    probabilities = np.array(
        [
            0.2628137141616739,
            0.22022352924182012,
            0.21384317109559658,
            0.1215805211501124,
            0.03436378945072717,
            0.10118879327456494,
            0.2860518847262805,
        ]
    )
    syndrome = np.array([1, 0, 0], dtype=np.uint8)
    decoder = QuotientBpDecoder(
        H,
        L,
        error_channel=probabilities,
        max_iter=8,
        bp_method="minimum_sum",
        ms_scaling_factor=0.8,
    )

    decoder.decode(syndrome)

    expected_trajectory = np.array(
        [[0], [1], [1], [0], [0], [0], [1], [0]], dtype=np.uint8
    )
    expected_unique = np.array([[0], [1]], dtype=np.uint8)
    np.testing.assert_array_equal(
        decoder.trajectory_logical_cosets, expected_trajectory
    )
    np.testing.assert_array_equal(
        decoder.candidate_logical_cosets, expected_unique
    )
    assert decoder.h_iter == 8


@pytest.mark.parametrize(
    ("word_schedule", "integer_schedule"),
    [("parallel", 0), ("serial", 1), ("serial_relative", 2)],
)
def test_schedule_integer_alias_matches_existing_bp_decoder_api(
    word_schedule,
    integer_schedule,
):
    H = np.array([[1, 1, 1, 1], [0, 1, 0, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 0, 0]], dtype=np.uint8)
    syndrome = np.array([1, 0], dtype=np.uint8)
    kwargs = {
        "error_channel": np.array([0.08, 0.12, 0.16, 0.20]),
        "max_iter": 5,
    }

    word_decoder = QuotientBpDecoder(H, L, schedule=word_schedule, **kwargs)
    integer_decoder = QuotientBpDecoder(H, L, schedule=integer_schedule, **kwargs)
    word_decoder.decode(syndrome)
    integer_decoder.decode(syndrome)

    np.testing.assert_array_equal(
        integer_decoder.trajectory_logical_cosets,
        word_decoder.trajectory_logical_cosets,
    )


def test_serial_relative_repeated_decode_is_deterministic():
    H = np.array(
        [
            [1, 0, 1, 1, 0, 1],
            [1, 1, 1, 1, 0, 0],
            [1, 1, 0, 1, 0, 1],
        ],
        dtype=np.uint8,
    )
    L = np.array(
        [
            [1, 0, 1, 1, 1, 0],
            [1, 0, 1, 0, 0, 1],
            [1, 0, 0, 1, 1, 0],
        ],
        dtype=np.uint8,
    )
    syndrome = np.array([1, 0, 1], dtype=np.uint8)
    decoder = QuotientBpDecoder(
        H,
        L,
        error_rate=0.1,
        max_iter=7,
        schedule="serial_relative",
        bp_method="minimum_sum",
    )

    first_decoding = decoder.decode(syndrome)
    first_trajectory = decoder.trajectory_logical_cosets
    first_candidates = decoder.candidate_logical_cosets
    first_scores = decoder.candidate_log_partitions
    first_chosen_label = decoder.chosen_logical_coset

    second_decoding = decoder.decode(syndrome)

    np.testing.assert_array_equal(second_decoding, first_decoding)
    np.testing.assert_array_equal(
        decoder.trajectory_logical_cosets, first_trajectory
    )
    np.testing.assert_array_equal(
        decoder.candidate_logical_cosets, first_candidates
    )
    np.testing.assert_allclose(
        decoder.candidate_log_partitions, first_scores, rtol=0, atol=0
    )
    np.testing.assert_array_equal(
        decoder.chosen_logical_coset, first_chosen_label
    )


def test_d_bp_bethe_partition_is_exact_on_a_tree_factor_graph():
    H = np.array([[1, 1, 1, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 1, 0]], dtype=np.uint8)
    probabilities = np.array([0.08, 0.11, 0.17, 0.23])
    syndrome = np.array([1], dtype=np.uint8)
    decoder = QuotientBpDecoder(
        H,
        L,
        error_channel=probabilities,
        max_iter=5,
        bp_method="product_sum",
        d_max_iter=60,
        d_damping=0.2,
        d_tolerance=1e-12,
    )

    correction = decoder.decode(syndrome)
    D = decoder.D.toarray()
    expected_scores = np.array(
        [
            exact_log_partition(D, offset, probabilities)
            for offset in decoder.candidate_representatives
        ]
    )

    np.testing.assert_allclose(
        decoder.candidate_log_partitions, expected_scores, atol=1e-10, rtol=0
    )
    np.testing.assert_array_equal(gf2_product(H, correction), syndrome)
    np.testing.assert_array_equal(
        gf2_product(L, correction), decoder.chosen_logical_coset
    )
    assert np.argmax(decoder.candidate_log_partitions) == np.flatnonzero(
        np.all(
            decoder.candidate_logical_cosets == decoder.chosen_logical_coset,
            axis=1,
        )
    )[0]


def test_each_d_bp_candidate_has_requested_syndrome_and_logical_label():
    H = np.array(
        [
            [0, 0, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=np.uint8,
    )
    L = np.array([[1, 0, 1, 1, 0, 1, 0]], dtype=np.uint8)
    probabilities = np.array([0.26, 0.22, 0.21, 0.12, 0.034, 0.10, 0.29])
    syndrome = np.array([1, 0, 0], dtype=np.uint8)
    decoder = QuotientBpDecoder(
        scipy.sparse.csr_matrix(H),
        scipy.sparse.csr_matrix(L),
        error_channel=probabilities,
        max_iter=8,
        bp_method="minimum_sum",
        ms_scaling_factor=0.8,
    )

    decoder.decode(syndrome)

    for representative, label in zip(
        decoder.candidate_representatives,
        decoder.candidate_logical_cosets,
    ):
        np.testing.assert_array_equal(gf2_product(H, representative), syndrome)
        np.testing.assert_array_equal(gf2_product(L, representative), label)


def test_zero_degeneracy_dimension_reduces_to_representative_probability():
    H = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 0]], dtype=np.uint8)
    probabilities = np.array([0.1, 0.2, 0.3])
    syndrome = np.array([1, 0], dtype=np.uint8)
    decoder = QuotientBpDecoder(
        H,
        L,
        error_channel=probabilities,
        max_iter=4,
        d_max_iter=10,
    )

    decoder.decode(syndrome)

    assert decoder.D.shape == (3, 0)
    expected = []
    for offset in decoder.candidate_representatives:
        expected.append(
            np.log(np.prod(np.where(offset == 1, probabilities, 1 - probabilities)))
        )
    np.testing.assert_allclose(
        decoder.candidate_log_partitions, np.asarray(expected), atol=1e-12, rtol=0
    )


def test_rank_deficient_unreachable_trajectory_label_raises():
    H = np.array([[1, 1]], dtype=np.uint8)
    L = H.copy()
    decoder = QuotientBpDecoder(
        H,
        L,
        error_rate=0.1,
        max_iter=2,
        bp_method="product_sum",
    )

    assert decoder.logical_rank == 0
    with pytest.raises(ValueError, match="logical label is unreachable"):
        decoder.decode(np.array([1], dtype=np.uint8))


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"error_rate": 0.0}, "error_rate"),
        ({"error_rate": 1.0}, "error_rate"),
        ({"error_rate": 0.1, "d_damping": 1.0}, "damping"),
        ({"error_rate": 0.1, "d_tolerance": 0.0}, "tolerance"),
    ],
)
def test_invalid_parameters_raise(kwargs, message):
    H = np.array([[1, 1, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 0]], dtype=np.uint8)

    with pytest.raises((ValueError, TypeError), match=message):
        QuotientBpDecoder(H, L, **kwargs)
