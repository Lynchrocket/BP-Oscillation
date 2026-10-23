from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest
import scipy.sparse


PROJECT_ROOT = Path(__file__).resolve().parents[6]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from experiment.utils.quotient_split import (
    QuotientSplitDecoder as PythonQuotientSplitDecoder,
    logical_neutral_representative,
)
from ldpc_modified.quotient_split_decoder import QuotientSplitDecoder


def gf2_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (left @ right) % 2


def python_reference(decoder, H, L, syndrome, probabilities, **kwargs):
    D = decoder.D.toarray()
    Q = decoder.Q.toarray()
    x_s = logical_neutral_representative(H, L, Q, syndrome)
    return PythonQuotientSplitDecoder(
        D,
        Q,
        x_s,
        probabilities,
        **kwargs,
    ).decode()


@pytest.mark.parametrize("adaptive_split", [False, True])
def test_cpp_split_matches_python_reference(adaptive_split):
    H = np.array([[1, 1, 1, 1, 1]], dtype=np.uint8)
    L = np.array(
        [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]],
        dtype=np.uint8,
    )
    probabilities = np.array([0.06, 0.09, 0.13, 0.18, 0.24])
    syndrome = np.array([1], dtype=np.uint8)
    kwargs = {
        "list_width": 3,
        "split_order": [1, 0],
        "adaptive_split": adaptive_split,
        "internal_iterations": 8,
        "leaf_iterations": 30,
        "damping": 0.2,
        "tolerance": 1e-12,
    }
    decoder = QuotientSplitDecoder(
        scipy.sparse.csr_matrix(H),
        scipy.sparse.csr_matrix(L),
        error_channel=probabilities,
        **kwargs,
    )

    correction = decoder.decode(syndrome)
    reference = python_reference(
        decoder,
        H,
        L,
        syndrome,
        probabilities,
        **kwargs,
    )

    expected_labels = np.asarray(
        [leaf.logical_label for leaf in reference.leaves],
        dtype=np.uint8,
    )
    expected_corrections = np.asarray(
        [leaf.correction for leaf in reference.leaves],
        dtype=np.uint8,
    )
    expected_path_scores = np.asarray(
        [leaf.path_score for leaf in reference.leaves],
        dtype=float,
    )
    expected_log_partitions = np.asarray(
        [leaf.log_z_bethe for leaf in reference.leaves],
        dtype=float,
    )

    np.testing.assert_array_equal(decoder.best_logical_label, reference.best_logical_label)
    np.testing.assert_array_equal(correction, reference.correction)
    np.testing.assert_array_equal(decoder.leaf_logical_labels, expected_labels)
    np.testing.assert_array_equal(decoder.leaf_corrections, expected_corrections)
    np.testing.assert_allclose(
        decoder.leaf_path_scores,
        expected_path_scores,
        atol=1e-12,
        rtol=0,
    )
    np.testing.assert_allclose(
        decoder.leaf_log_partitions,
        expected_log_partitions,
        atol=1e-12,
        rtol=0,
    )
    np.testing.assert_array_equal(
        decoder.surviving_prefix_counts,
        reference.surviving_prefix_counts,
    )


def test_full_width_keeps_every_logical_label_and_valid_representatives():
    H = np.array([[1, 1, 1, 1, 1]], dtype=np.uint8)
    L = np.array(
        [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]],
        dtype=np.uint8,
    )
    probabilities = np.array([0.06, 0.09, 0.13, 0.18, 0.24])
    syndrome = np.array([1], dtype=np.uint8)
    decoder = QuotientSplitDecoder(
        H,
        L,
        error_channel=probabilities,
        list_width=4,
        internal_iterations=5,
        leaf_iterations=20,
    )

    decoder.decode(syndrome)

    assert {
        tuple(map(int, label)) for label in decoder.leaf_logical_labels
    } == set(product((0, 1), repeat=2))
    np.testing.assert_array_equal(decoder.surviving_prefix_counts, [1, 2, 4])
    assert np.all(decoder.leaf_log_partitions[:-1] >= decoder.leaf_log_partitions[1:])
    for correction, label in zip(
        decoder.leaf_corrections,
        decoder.leaf_logical_labels,
    ):
        np.testing.assert_array_equal(gf2_product(H, correction), syndrome)
        np.testing.assert_array_equal(gf2_product(L, correction), label)


def test_zero_logical_dimension_runs_one_d_only_leaf():
    H = np.array([[1, 1, 1]], dtype=np.uint8)
    L = np.zeros((0, 3), dtype=np.uint8)
    syndrome = np.array([1], dtype=np.uint8)
    decoder = QuotientSplitDecoder(
        H,
        L,
        error_channel=np.array([0.1, 0.2, 0.3]),
        list_width=3,
    )

    correction = decoder.decode(syndrome)

    assert decoder.leaf_logical_labels.shape == (1, 0)
    assert decoder.leaf_corrections.shape == (1, 3)
    np.testing.assert_array_equal(decoder.surviving_prefix_counts, [1])
    np.testing.assert_array_equal(gf2_product(H, correction), syndrome)


def test_zero_degeneracy_dimension_and_channel_update_are_deterministic():
    H = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 0]], dtype=np.uint8)
    syndrome = np.array([1, 0], dtype=np.uint8)
    decoder = QuotientSplitDecoder(
        H,
        L,
        error_channel=np.array([0.1, 0.2, 0.3]),
        list_width=2,
        split_order=[0],
    )

    first = decoder.decode(syndrome)
    first_labels = decoder.leaf_logical_labels.copy()
    first_scores = decoder.leaf_log_partitions.copy()
    np.testing.assert_array_equal(decoder.decode(syndrome), first)
    np.testing.assert_array_equal(decoder.leaf_logical_labels, first_labels)
    np.testing.assert_array_equal(decoder.leaf_log_partitions, first_scores)
    assert decoder.D.shape == (3, 0)

    decoder.update_channel_probs(np.array([0.25, 0.15, 0.05]))
    decoder.decode(syndrome)
    assert not np.array_equal(decoder.leaf_log_partitions, first_scores)


def test_invalid_channel_update_does_not_mutate_decoder_state():
    H = np.array([[1, 1, 1, 1]], dtype=np.uint8)
    L = np.array([[1, 0, 0, 0]], dtype=np.uint8)
    syndrome = np.array([1], dtype=np.uint8)
    probabilities = np.array([0.08, 0.12, 0.16, 0.20])
    decoder = QuotientSplitDecoder(
        H,
        L,
        error_channel=probabilities,
        list_width=2,
    )

    expected_correction = decoder.decode(syndrome)
    expected_scores = decoder.leaf_log_partitions.copy()
    with pytest.raises(ValueError, match="error_channel"):
        decoder.update_channel_probs([0.1, 0.2])

    np.testing.assert_array_equal(decoder.error_channel, probabilities)
    np.testing.assert_array_equal(decoder.decode(syndrome), expected_correction)
    np.testing.assert_array_equal(decoder.leaf_log_partitions, expected_scores)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"error_rate": 0.0}, "error_rate"),
        ({"error_rate": 0.1, "list_width": 0}, "list_width"),
        ({"error_rate": 0.1, "split_order": [0, 0]}, "split_order"),
        ({"error_rate": 0.1, "internal_iterations": -1}, "internal_iterations"),
        ({"error_rate": 0.1, "leaf_iterations": -1}, "leaf_iterations"),
        ({"error_rate": 0.1, "damping": 1.0}, "damping"),
        ({"error_rate": 0.1, "tolerance": 0.0}, "tolerance"),
        ({"error_rate": 0.1, "probability_floor": 0.5}, "probability_floor"),
    ],
)
def test_invalid_parameters_raise(kwargs, message):
    H = np.array([[1, 1, 1, 1, 1]], dtype=np.uint8)
    L = np.array(
        [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]],
        dtype=np.uint8,
    )

    with pytest.raises((ValueError, TypeError), match=message):
        QuotientSplitDecoder(H, L, **kwargs)


def test_rank_deficient_logical_map_raises_directly():
    H = np.array([[1, 1]], dtype=np.uint8)
    L = H.copy()

    with pytest.raises(ValueError, match="full logical rank"):
        QuotientSplitDecoder(H, L, error_rate=0.1)
