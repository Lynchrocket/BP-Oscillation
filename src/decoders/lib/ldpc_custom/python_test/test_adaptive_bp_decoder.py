import gc
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp")

import numpy as np
import pytest

import ldpc_modified
from ldpc_modified import AdaptiveBpDecoder
from ldpc_modified.adaptive_bp_decoder._adaptive_bp_decoder import (
    _native_live_instance_count,
)


CHAIN_PCM = np.array(
    [
        [1, 1, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 1, 1],
    ],
    dtype=np.uint8,
)


def syndrome_for(pcm: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    return ((pcm @ candidate) % 2).astype(np.uint8)


def channel_nll(candidate: np.ndarray, probabilities: np.ndarray) -> float:
    return float(
        np.sum(
            np.where(
                candidate,
                -np.log(probabilities),
                -np.log1p(-probabilities),
            )
        )
    )


def make_chain_decoder(**kwargs) -> AdaptiveBpDecoder:
    parameters = {
        "error_rate": 0.1,
        "max_iter": 10,
        "bp_method": "minimum_sum",
        "schedule": "parallel",
        "input_vector_type": "syndrome",
        "adaptive_k": 4,
        "adaptive_max_views": 4,
    }
    parameters.update(kwargs)
    return AdaptiveBpDecoder(CHAIN_PCM, **parameters)


def test_top_level_export_uses_dedicated_adaptive_bp_module():
    assert ldpc_modified.AdaptiveBpDecoder.__module__ == (
        "ldpc_modified.adaptive_bp_decoder._adaptive_bp_decoder"
    )
    assert ldpc_modified.AdaptiveBpOsdDecoder.__module__ == (
        "ldpc_modified.adaptive_bposd_decoder._adaptive_bposd_decoder"
    )


@pytest.mark.parametrize(
    "name",
    ["adaptive_k", "adaptive_max_views"],
)
def test_negative_adaptive_parameter_is_rejected(name: str):
    with pytest.raises(ValueError, match="must be non-negative"):
        make_chain_decoder(**{name: -1})


@pytest.mark.parametrize("name", ["adaptive_k", "adaptive_max_views"])
def test_negative_mutated_adaptive_parameter_is_rejected_before_planning(name: str):
    decoder = make_chain_decoder()
    setattr(decoder, name, -1)

    with pytest.raises(RuntimeError, match="must be non-negative"):
        decoder.decode(np.array([1, 0, 0], dtype=np.uint8))


@pytest.mark.parametrize("name", ["adaptive_w_max", "adaptive_B"])
def test_removed_adaptive_parameter_is_rejected(name: str):
    with pytest.raises(ValueError, match="were removed"):
        make_chain_decoder(**{name: 3})


def test_constructor_errors_after_native_allocation_are_reported_cleanly():
    baseline = _native_live_instance_count()

    with pytest.raises(ValueError, match="length 4"):
        AdaptiveBpDecoder(CHAIN_PCM, error_channel=[0.1, 0.2])
    gc.collect()
    assert _native_live_instance_count() == baseline

    with pytest.raises(ValueError, match="schedule"):
        make_chain_decoder(schedule="not-a-schedule")
    gc.collect()
    assert _native_live_instance_count() == baseline


def test_multiview_decode_covers_feasible_unstable_variables():
    decoder = make_chain_decoder()
    syndrome = np.array([1, 0, 0], dtype=np.uint8)

    decoding = decoder.decode(syndrome)
    stats = decoder.adaptive_stats

    np.testing.assert_array_equal(syndrome_for(CHAIN_PCM, decoding), syndrome)
    assert stats["prepared_view_count"] == 2
    assert stats["executed_view_count"] == 2
    assert stats["selected_rows_per_view"] == [2, 1]
    assert stats["automatic_weight_limit_per_view"] == [2, 4]
    assert stats["selected_row_weights_per_view"] == [[2, 2], [4]]
    assert stats["covered_unstable_variable_count"] == 4
    assert stats["feasible_unstable_variable_count"] == 4
    assert stats["candidate_count"] == 3
    assert stats["valid_candidate_count"] == 3
    assert stats["selected_candidate_valid"] is True
    assert stats["bp_run_count"] == 3


def test_automatic_planner_covers_all_feasible_unstable_variables():
    pcm = np.array(
        [
            [0, 1, 0, 1, 1, 1, 1],
            [0, 0, 0, 1, 0, 1, 0],
            [0, 1, 0, 0, 1, 0, 1],
            [1, 1, 1, 0, 1, 1, 1],
            [0, 1, 1, 0, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    syndrome = np.array([1, 0, 1, 1, 0], dtype=np.uint8)
    decoder = AdaptiveBpDecoder(
        pcm,
        error_rate=0.1,
        max_iter=2,
        input_vector_type="syndrome",
        adaptive_k=7,
        adaptive_max_views=3,
    )

    decoder.decode(syndrome)
    stats = decoder.adaptive_stats

    assert stats["feasible_unstable_variable_count"] == 7
    assert stats["covered_unstable_variable_count"] == 7
    assert all(rows <= 7 for rows in stats["selected_rows_per_view"])
    assert all(
        weight <= limit
        for limit, weights in zip(
            stats["automatic_weight_limit_per_view"],
            stats["selected_row_weights_per_view"],
        )
        for weight in weights
    )


def test_adaptive_max_views_zero_runs_main_bp_only():
    decoder = make_chain_decoder(adaptive_max_views=0)
    syndrome = np.array([1, 0, 0], dtype=np.uint8)

    decoder.decode(syndrome)
    stats = decoder.adaptive_stats

    assert stats["prepared_view_count"] == 0
    assert stats["executed_view_count"] == 0
    assert stats["candidate_count"] == 1
    assert stats["bp_run_count"] == 1
    assert decoder.auxiliary_converge is False
    np.testing.assert_array_equal(decoder.decoding, decoder.main_decoding)


def test_main_only_lazy_instability_matches_multiview_main_bp():
    syndrome = np.array([1, 0, 0], dtype=np.uint8)
    main_only_decoder = make_chain_decoder(adaptive_max_views=0)
    multiview_decoder = make_chain_decoder()

    main_only_decoder.decode(syndrome)
    multiview_decoder.decode(syndrome)

    np.testing.assert_allclose(
        main_only_decoder.instability_scores,
        multiview_decoder.instability_scores,
    )


def test_nonzero_received_vector_runs_main_bp_only():
    decoder = make_chain_decoder(input_vector_type="received_vector")
    received = np.array([1, 0, 0, 0], dtype=np.uint8)

    decoding = decoder.decode(received)
    stats = decoder.adaptive_stats

    assert stats["prepared_view_count"] == 0
    assert stats["executed_view_count"] == 0
    assert stats["candidate_count"] == 1
    assert stats["selected_view"] == 0
    assert stats["bp_run_count"] == 1
    assert decoder.auxiliary_converge is False
    assert decoder.view_candidates.shape == (0, CHAIN_PCM.shape[1])
    assert decoder.view_iterations.size == 0
    np.testing.assert_array_equal(decoding, decoder.main_decoding)


def test_received_vector_clears_previous_auxiliary_diagnostics():
    decoder = make_chain_decoder(input_vector_type="auto")
    decoder.decode(np.array([1, 0, 0], dtype=np.uint8))
    assert decoder.adaptive_stats["executed_view_count"] > 0

    decoder.decode(np.array([1, 0, 0, 0], dtype=np.uint8))

    assert decoder.adaptive_stats["executed_view_count"] == 0
    assert decoder.auxiliary_converge is False
    assert decoder.view_candidates.shape == (0, CHAIN_PCM.shape[1])
    np.testing.assert_array_equal(
        decoder.auxiliary_decoding,
        np.zeros(CHAIN_PCM.shape[1], dtype=int),
    )
    np.testing.assert_allclose(
        decoder.auxiliary_log_prob_ratios,
        np.full(CHAIN_PCM.shape[1], np.log(9.0)),
    )


def test_zero_received_vector_uses_zero_run_fast_path():
    decoder = make_chain_decoder(input_vector_type="received_vector")

    decoding = decoder.decode(np.zeros(CHAIN_PCM.shape[1], dtype=np.uint8))
    stats = decoder.adaptive_stats

    np.testing.assert_array_equal(decoding, np.zeros(CHAIN_PCM.shape[1], dtype=np.uint8))
    assert stats["prepared_view_count"] == 0
    assert stats["executed_view_count"] == 0
    assert stats["candidate_count"] == 0
    assert stats["valid_candidate_count"] == 0
    assert stats["selected_view"] == -1
    assert stats["bp_run_count"] == 0
    assert decoder.auxiliary_converge is False


def test_valid_candidates_are_selected_by_original_channel_nll():
    pcm = np.array(
        [
            [0, 1, 0, 0, 0, 1, 1, 0],
            [0, 0, 0, 0, 0, 1, 0, 0],
            [1, 0, 0, 1, 1, 1, 0, 1],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 0, 0],
        ],
        dtype=np.uint8,
    )
    probabilities = np.array(
        [
            0.23072770,
            0.14634845,
            0.17615726,
            0.10347890,
            0.23981567,
            0.26389831,
            0.11157496,
            0.16439452,
        ]
    )
    syndrome = np.ones(5, dtype=np.uint8)
    decoder = AdaptiveBpDecoder(
        pcm,
        error_channel=probabilities,
        max_iter=2,
        bp_method="minimum_sum",
        schedule="parallel",
        input_vector_type="syndrome",
        adaptive_k=6,
        adaptive_max_views=1,
    )

    decoding = decoder.decode(syndrome)
    main_decoding = decoder.main_decoding
    auxiliary_decoding = decoder.auxiliary_decoding
    stats = decoder.adaptive_stats

    np.testing.assert_array_equal(syndrome_for(pcm, main_decoding), syndrome)
    np.testing.assert_array_equal(syndrome_for(pcm, auxiliary_decoding), syndrome)
    assert channel_nll(auxiliary_decoding, probabilities) < channel_nll(
        main_decoding, probabilities
    )
    np.testing.assert_array_equal(decoding, auxiliary_decoding)
    assert stats["valid_candidate_count"] == 2
    assert stats["selected_view"] == 1
    assert np.isclose(
        stats["best_candidate_cost"],
        channel_nll(auxiliary_decoding, probabilities),
    )
    # Captured before removing the redundant global-gain and row-index keys.
    # These values lock planner ordering across the lossless simplification.
    np.testing.assert_array_equal(
        main_decoding,
        np.array([1, 0, 1, 0, 1, 1, 0, 0]),
    )
    np.testing.assert_array_equal(
        auxiliary_decoding,
        np.array([0, 0, 1, 0, 0, 1, 0, 0]),
    )
    np.testing.assert_array_equal(
        decoder.log_prob_ratios,
        np.array(
            [
                -0.9752880164326114,
                0.7377282320172338,
                -1.7976931348623157e308,
                1.1333496661557159,
                -1.076320854363088,
                -1.7976931348623157e308,
                1.0489491897152243,
                0.600083003691678,
            ]
        ),
    )
    np.testing.assert_array_equal(
        decoder.auxiliary_log_prob_ratios,
        np.array(
            [
                1.2042067757311388,
                1.7635326674150835,
                -1.7976931348623157e308,
                2.2870400229216163,
                1.1536903567659005,
                -1.7976931348623157e308,
                2.074753625113074,
                1.7537733604575785,
            ]
        ),
    )
    np.testing.assert_array_equal(
        decoder.instability_scores,
        np.array(
            [
                3.025338139247065,
                2.3555127167415817,
                2.0,
                1.8823402254937651,
                2.9290909824383364,
                2.0,
                1.9533350230915696,
                2.6664361327456167,
            ]
        ),
    )
    assert stats == {
        "prepared_view_count": 1,
        "executed_view_count": 1,
        "selected_rows_per_view": [3],
        "automatic_weight_limit_per_view": [4],
        "selected_row_weights_per_view": [[4, 4, 4]],
        "unstable_variable_count": 6,
        "feasible_unstable_variable_count": 6,
        "covered_unstable_variable_count": 6,
        "candidate_pair_count": 10,
        "candidate_count": 2,
        "valid_candidate_count": 2,
        "selected_view": 1,
        "best_candidate_cost": 4.170443543956942,
        "fallback_residual_weight": 0,
        "selected_candidate_valid": True,
        "bp_run_count": 2,
    }


def test_zero_syndrome_resets_multiview_state():
    decoder = make_chain_decoder()
    decoder.decode(np.array([1, 0, 0], dtype=np.uint8))

    decoding = decoder.decode(np.zeros(CHAIN_PCM.shape[0], dtype=np.uint8))
    stats = decoder.adaptive_stats

    np.testing.assert_array_equal(decoding, np.zeros(CHAIN_PCM.shape[1], dtype=np.uint8))
    assert stats["prepared_view_count"] == 0
    assert stats["executed_view_count"] == 0
    assert stats["candidate_count"] == 1
    assert stats["valid_candidate_count"] == 1
    assert stats["selected_view"] == 0
    assert stats["selected_candidate_valid"] is True
    assert stats["bp_run_count"] == 0


def test_no_valid_candidate_falls_back_to_main_bp_decoding():
    pcm = np.array([[1, 0], [1, 0]], dtype=np.uint8)
    decoder = AdaptiveBpDecoder(
        pcm,
        error_rate=0.1,
        max_iter=2,
        input_vector_type="syndrome",
        adaptive_k=2,
        adaptive_max_views=2,
    )

    decoding = decoder.decode(np.array([1, 0], dtype=np.uint8))
    stats = decoder.adaptive_stats

    np.testing.assert_array_equal(decoding, decoder.main_decoding)
    assert stats["valid_candidate_count"] == 0
    assert stats["selected_view"] == 0
    assert stats["selected_candidate_valid"] is False
    assert np.isinf(stats["best_candidate_cost"])
    assert stats["fallback_residual_weight"] > 0


def test_repeated_decode_is_deterministic():
    decoder = make_chain_decoder()
    syndrome = np.array([1, 0, 0], dtype=np.uint8)

    first = decoder.decode(syndrome)
    first_stats = decoder.adaptive_stats
    second = decoder.decode(syndrome)
    second_stats = decoder.adaptive_stats

    np.testing.assert_array_equal(first, second)
    assert first_stats == second_stats
