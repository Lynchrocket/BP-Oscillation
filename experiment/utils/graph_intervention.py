# instability_guided_bp
"""Reusable experiment functions with explicit arguments and dependencies."""
from .bp_trace import BPContext
from itertools import combinations
from dataclasses import dataclass
import numpy as np
import pandas as pd

@dataclass(frozen=True)
class GraphContext(BPContext):
    p: float
    H: np.ndarray
    L: np.ndarray
    priors: np.ndarray
    prior_llr: np.ndarray
    edge_checks: np.ndarray
    edge_variables: np.ndarray
    check_edges: tuple[np.ndarray, ...]
    variable_edges: tuple[np.ndarray, ...]

def graph_context(p, H, L, priors, *, GraphContext):
    H = np.asarray(H, dtype=np.uint8) % 2
    L = np.asarray(L, dtype=np.uint8) % 2
    priors = np.asarray(priors, dtype=float)
    assert H.ndim == 2 and L.shape[1] == H.shape[1] == priors.size
    assert np.all((priors > 0) & (priors < 1))
    checks, variables = np.nonzero(H)
    checks = checks.astype(np.int32)
    variables = variables.astype(np.int32)
    edge_ids = np.arange(len(checks), dtype=np.int32)
    check_edges = tuple((edge_ids[checks == check] for check in range(H.shape[0])))
    variable_edges = tuple((edge_ids[variables == variable] for variable in range(H.shape[1])))
    assert all((len(ids) >= 2 for ids in check_edges))
    return GraphContext(float(p), H, L, priors, np.log((1 - priors) / priors), checks, variables, check_edges, variable_edges)

def context_for(level, p, *, CA_9_2_5, CONTEXT_CACHE, CircuitBuilder, ErrorRate, TriangularCode, detector_error_model_to_check_matrices, graph_context):
    key = (int(level), float(p))
    if key in CONTEXT_CACHE:
        return CONTEXT_CACHE[key]
    code_object = CA_9_2_5() if level == 2 else TriangularCode(level)
    circuit = CircuitBuilder(code_object).get_circuit_level_circuit(ErrorRate(p=float(p), p_meas=0), d_rounds=1)
    matrices = detector_error_model_to_check_matrices(circuit.detector_error_model(), allow_undecomposed_hyperedges=True)
    context = graph_context(p, matrices.check_matrix.toarray(), matrices.observables_matrix.toarray(), matrices.priors)
    distance = 5 if level == 2 else code_object.d
    assert distance in {3, 5, 9}
    CONTEXT_CACHE[key] = (distance, context)
    return CONTEXT_CACHE[key]

def augmented_context(context, extra_rows, *, graph_context):
    extra_rows = np.asarray(extra_rows, dtype=np.uint8).reshape(-1, context.H.shape[1])
    return graph_context(context.p, np.vstack([context.H, extra_rows]), context.L, context.priors)

def logical_label(context, correction, *, bit_text):
    return bit_text(context.L @ np.asarray(correction, dtype=np.uint8) % 2)

def trace_min_sum(context, syndrome, iterations, edge_damping=None, edge_scaling=None, *, CONFIG, MessageTrace):
    trace = context.bp_decoder.trace(
        syndrome, iterations,
        scaling=CONFIG['ms_scaling_factor'] if edge_scaling is None else edge_scaling,
        damping=0.0 if edge_damping is None else edge_damping,
        accumulation='sum_then_prior',
    )
    return MessageTrace(trace.c2v, trace.v2c, trace.llr, trace.hard, trace.valid, trace.first_valid)

def percentile_rank(values):
    values = np.asarray(values, dtype=float)
    order = np.lexsort((np.arange(len(values)), values))
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks / max(1, len(values) - 1)

def outcome_of(context, syndrome, truth, correction, *, logical_label):
    correction = np.asarray(correction, dtype=np.uint8)
    valid = np.array_equal(context.H @ correction % 2, syndrome[:context.H.shape[0]])
    correct = logical_label(context, correction) == truth
    return (valid, correct)

def top_region(scores, count):
    assert 0 < count <= len(scores)
    return np.lexsort((np.arange(len(scores)), -np.asarray(scores)))[:count]

def conditional_outcome_rates(summary):
    summary = summary.copy()
    summary['repair_rate'] = summary.repair_count / summary.baseline_failure_evaluations.replace(0, np.nan)
    summary['regression_rate'] = summary.regression_count / summary.baseline_success_evaluations.replace(0, np.nan)
    return summary.drop(columns=['repair_count', 'regression_count'])

def edge_maps(context, variable_scores, low, high, mode, rng=None, *, percentile_rank):
    normalized = percentile_rank(variable_scores)[context.edge_variables]
    if mode == 'targeted':
        values = low + (high - low) * normalized
    elif mode == 'random':
        assert rng is not None
        values = low + (high - low) * rng.permutation(normalized)
    elif mode == 'global':
        values = np.full(len(normalized), np.mean(low + (high - low) * normalized))
    else:
        raise AssertionError(mode)
    return values

@dataclass(frozen=True)
class CandidateRow:
    candidate_id: int
    source_rows: tuple[int, int]
    row: np.ndarray
    row_weight: int

def generate_redundant_candidates(H, max_row_weight, *, CandidateRow, bit_text):
    existing = {bit_text(row) for row in H}
    unique = {}
    for row_a, row_b in combinations(range(H.shape[0]), 2):
        row = H[row_a] ^ H[row_b]
        text = bit_text(row)
        weight = int(row.sum())
        if weight and weight <= max_row_weight and (text not in existing):
            unique.setdefault(text, (row_a, row_b, row.copy(), weight))
    candidates = [CandidateRow(index, (a, b), row, weight) for index, (a, b, row, weight) in enumerate(unique.values())]
    return candidates

def augmented_syndrome(syndrome, source_rows):
    return np.concatenate([syndrome, np.asarray([int(syndrome[source_rows[0]] ^ syndrome[source_rows[1]])], dtype=np.uint8)])

def candidate_overlap(candidate, variable_scores, region):
    support = np.flatnonzero(candidate.row)
    return (int(np.intersect1d(support, region).size), float(variable_scores[support].sum()))
