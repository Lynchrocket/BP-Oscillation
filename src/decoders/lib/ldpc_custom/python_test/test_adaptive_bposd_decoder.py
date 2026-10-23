import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp")

import numpy as np

from ldpc_modified.adaptive_bposd_decoder import AdaptiveBpOsdDecoder


def syndrome_for(parity_check_matrix: np.ndarray, error: np.ndarray) -> np.ndarray:
    return ((parity_check_matrix @ error) % 2).astype(np.uint8)


def assert_satisfies_syndrome(
    parity_check_matrix: np.ndarray, decoding: np.ndarray, syndrome: np.ndarray
) -> None:
    np.testing.assert_array_equal(syndrome_for(parity_check_matrix, decoding), syndrome)


def make_decoder(parity_check_matrix: np.ndarray, adaptive_B: int) -> AdaptiveBpOsdDecoder:
    return AdaptiveBpOsdDecoder(
        parity_check_matrix,
        error_rate=0.1,
        max_iter=1,
        bp_method="minimum_sum",
        schedule="parallel",
        osd_method="osd_0",
        adaptive_k=4,
        adaptive_w_max=4,
        adaptive_B=adaptive_B,
    )


def test_adaptive_decode_updates_auxiliary_bp_dimensions_between_runs():
    parity_check_matrix = np.array(
        [
            [1, 1, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 0],
            [0, 0, 1, 1, 0, 0],
            [0, 0, 0, 1, 1, 0],
            [0, 0, 0, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    decoder = make_decoder(parity_check_matrix, adaptive_B=3)

    first_syndrome = syndrome_for(
        parity_check_matrix, np.array([1, 0, 0, 0, 0, 0], dtype=np.uint8)
    )
    first_decoding = decoder.decode(first_syndrome)

    assert_satisfies_syndrome(parity_check_matrix, first_decoding, first_syndrome)
    assert decoder.adaptive_stats["auxiliary_selected_row_count"] == 3
    assert decoder.auxiliary_check_count == parity_check_matrix.shape[0] + 3
    assert decoder.auxiliary_bit_count == parity_check_matrix.shape[1]

    decoder.adaptive_B = 1
    second_syndrome = syndrome_for(
        parity_check_matrix, np.array([0, 0, 0, 0, 0, 1], dtype=np.uint8)
    )
    second_decoding = decoder.decode(second_syndrome)

    assert_satisfies_syndrome(parity_check_matrix, second_decoding, second_syndrome)
    assert decoder.adaptive_stats["auxiliary_selected_row_count"] == 1
    assert decoder.auxiliary_check_count == parity_check_matrix.shape[0] + 1
    assert decoder.auxiliary_bit_count == parity_check_matrix.shape[1]


def test_zero_syndrome_resets_adaptive_osd_stats():
    parity_check_matrix = np.array(
        [
            [1, 1, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 0],
            [0, 0, 1, 1, 0, 0],
            [0, 0, 0, 1, 1, 0],
            [0, 0, 0, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    decoder = make_decoder(parity_check_matrix, adaptive_B=3)

    nonzero_syndrome = syndrome_for(
        parity_check_matrix, np.array([0, 1, 0, 1, 0, 0], dtype=np.uint8)
    )
    decoder.decode(nonzero_syndrome)
    assert decoder.adaptive_stats["auxiliary_selected_row_count"] == 3

    zero_syndrome = np.zeros(parity_check_matrix.shape[0], dtype=np.uint8)
    zero_decoding = decoder.decode(zero_syndrome)

    np.testing.assert_array_equal(
        zero_decoding, np.zeros(parity_check_matrix.shape[1], dtype=np.uint8)
    )
    stats = decoder.adaptive_stats
    assert stats["auxiliary_view_disagreement"] is False
    assert stats["auxiliary_view_disagreement_weight"] == 0
    assert stats["auxiliary_candidate_pair_count"] == 0
    assert stats["auxiliary_selected_row_count"] == 0
    assert stats["auxiliary_bp_converge"] is True
    assert stats["merged_candidate_count"] == 1
    assert np.isclose(
        stats["best_candidate_cost"], -parity_check_matrix.shape[1] * np.log(0.9)
    )
