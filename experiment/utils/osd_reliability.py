# osd
"""Reusable experiment functions with explicit arguments and dependencies."""
from .bp_trace import BPContext
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from time import perf_counter
from scipy.stats import spearmanr

@dataclass(frozen=True)
class Context(BPContext):
    level: int
    distance: int
    p: float
    H: np.ndarray
    L: np.ndarray
    priors: np.ndarray
    prior_llr: np.ndarray
    dem: object
    edge_checks: np.ndarray
    edge_variables: np.ndarray
    check_edges: tuple
    variable_edges: tuple

@dataclass(frozen=True)
class BPTrace:
    c2v: np.ndarray
    v2c: np.ndarray
    llr: np.ndarray
    hard: np.ndarray
    valid: np.ndarray
    runtime_seconds: np.ndarray
    first_valid: int | None

def gf2_rank(matrix):
    work = np.asarray(matrix, dtype=np.uint8).copy() % 2
    rank = 0
    for column in range(work.shape[1]):
        pivots = np.flatnonzero(work[rank:, column])
        if not len(pivots):
            continue
        pivot = rank + int(pivots[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        for row in np.flatnonzero(work[:, column]):
            if row != rank:
                work[row] ^= work[rank]
        rank += 1
    return rank

def logical_label(context, correction):
    return context.L @ np.asarray(correction, dtype=np.uint8) % 2

def full_order(context, *, CONFIG):
    return int(CONFIG['full_osd_order_by_distance'][str(context.distance)])

def endpoint_of(trace):
    return int(trace.first_valid) if trace.first_valid is not None else len(trace.valid) - 1

def tesseract_predictions(context, syndromes, *, CONFIG, Tesseract):
    compiled = Tesseract(beam=CONFIG['tesseract_beam'], pq_limit=CONFIG['tesseract_pq_limit'], seed=CONFIG['random_seed']).compile_decoder_for_dem(dem=context.dem)
    packed = np.packbits(np.asarray(syndromes, dtype=np.uint8), axis=1, bitorder='little')
    started = perf_counter()
    prediction = compiled.decode_shots_bit_packed(bit_packed_detection_event_data=packed)
    elapsed = perf_counter() - started
    unpacked = np.unpackbits(prediction, axis=1, count=context.L.shape[0], bitorder='little').astype(np.uint8)
    return (unpacked, elapsed / max(1, len(syndromes)))

@dataclass(frozen=True)
class Context_variant_2(BPContext):
    level: int
    distance: int
    p: float
    H: np.ndarray
    L: np.ndarray
    priors: np.ndarray
    prior_llr: np.ndarray
    circuit: object
    dem: object
    edge_checks: np.ndarray
    edge_variables: np.ndarray
    check_edges: tuple[np.ndarray, ...]
    variable_edges: tuple[np.ndarray, ...]

def gf2_rank_variant_2(matrix):
    work = np.asarray(matrix, dtype=np.uint8).copy() % 2
    rank = 0
    for column in range(work.shape[1]):
        pivots = np.flatnonzero(work[rank:, column])
        if not len(pivots):
            continue
        pivot = rank + int(pivots[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        for row in np.flatnonzero(work[:, column]):
            if row != rank:
                work[row] ^= work[rank]
        rank += 1
        if rank == work.shape[0]:
            break
    return rank

def context_for(level, p, *, CA_9_2_5, CONTEXT_CACHE, CircuitBuilder, Context, ErrorRate, TriangularCode, detector_error_model_to_check_matrices, gf2_rank):
    key = (int(level), float(p))
    if key in CONTEXT_CACHE:
        return CONTEXT_CACHE[key]
    code = CA_9_2_5() if level == 2 else TriangularCode(level)
    circuit = CircuitBuilder(code).get_circuit_level_circuit(ErrorRate(p=float(p), p_meas=0), d_rounds=1)
    dem = circuit.detector_error_model()
    matrices = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
    H = matrices.check_matrix.toarray().astype(np.uint8) % 2
    L = matrices.observables_matrix.toarray().astype(np.uint8) % 2
    priors = np.asarray(matrices.priors, dtype=float)
    checks, variables = np.nonzero(H)
    checks, variables = (checks.astype(np.int32), variables.astype(np.int32))
    edge_ids = np.arange(len(checks), dtype=np.int32)
    check_edges = tuple((edge_ids[checks == check] for check in range(H.shape[0])))
    variable_edges = tuple((edge_ids[variables == variable] for variable in range(H.shape[1])))
    distance = 5 if level == 2 else int(code.d)
    assert distance in {5, 9} and H.shape[1] == L.shape[1] == priors.size
    assert gf2_rank(H) == H.shape[0] and all((len(ids) >= 2 for ids in check_edges))
    assert np.all((priors > 0) & (priors < 1))
    result = Context(level, distance, float(p), H, L, priors, np.log((1 - priors) / priors), circuit, dem, checks, variables, check_edges, variable_edges)
    CONTEXT_CACHE[key] = result
    return result

def trace_min_sum(context, syndrome, iterations, scaling=1.0, damping=0.0, schedule='parallel', early_stop=False, *, BPTrace):
    trace = context.bp_decoder.trace(
        syndrome, iterations, scaling=scaling, damping=damping, schedule=schedule,
        stop_early=early_stop, serial_accumulation='numpy', combine_scaling_damping=True,
        accumulation='sum_then_prior',
    )
    return BPTrace(trace.c2v, trace.v2c, trace.llr, trace.hard, trace.valid, trace.runtime_seconds, trace.first_valid)

def outcome(context, correction, syndrome, truth, error=None, *, logical_label):
    correction = np.asarray(correction, dtype=np.uint8)
    syndrome_consistent = np.array_equal(context.H @ correction % 2, syndrome)
    logical_correct = syndrome_consistent and np.array_equal(logical_label(context, correction), truth)
    exact = np.nan if error is None else int(np.array_equal(correction, error))
    return (int(syndrome_consistent), int(logical_correct), exact)

def ordering_proxy(reliability):
    reliability = np.asarray(reliability, dtype=float)
    assert reliability.ndim == 1 and np.isfinite(reliability).all()
    order = np.lexsort((np.arange(len(reliability)), reliability))
    proxy = np.empty(len(reliability), dtype=np.float64)
    proxy[order] = np.arange(len(reliability), dtype=np.float64)
    return proxy

def osd_from_reliability(context, syndrome, reliability, order, *, ordered_statistics_postprocess, ordering_proxy):
    correction = ordered_statistics_postprocess(np.ascontiguousarray(ordering_proxy(reliability)), np.ascontiguousarray(syndrome, dtype=np.uint8), np.ascontiguousarray(context.H, dtype=np.uint8), int(order), np.ascontiguousarray(context.priors, dtype=np.float64))
    assert np.array_equal(context.H @ correction % 2, syndrome)
    return correction

def ordinal_ranks(reliability):
    order = np.lexsort((np.arange(len(reliability)), np.asarray(reliability, dtype=float)))
    ranks = np.empty(len(order), dtype=int)
    ranks[order] = np.arange(1, len(order) + 1)
    return ranks

def safe_corr(x, y, kind):
    x, y = (np.asarray(x, float), np.asarray(y, float))
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float((spearmanr if kind == 'spearman' else pearsonr)(x, y).statistic)
