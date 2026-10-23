# gap_mechanisms, osd
"""Reusable experiment functions with explicit arguments and dependencies."""
from .bp_trace import BPContext
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.sparse import csc_matrix
from dataclasses import dataclass
import json
from scipy.optimize import milp
import numpy as np
import pandas as pd
from time import perf_counter
from dataclasses import replace
from hashlib import sha256
from tqdm.auto import tqdm

def split_for(shot_id, *, CONFIG):
    value = int.from_bytes(sha256(f"{CONFIG['random_seed']}:{shot_id}".encode()).digest()[:8], 'little') / 2 ** 64
    return 'development' if value < CONFIG['development_fraction'] else 'test'

@dataclass(frozen=True)
class CodeContext(BPContext):
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
    detector_coordinates: dict

@dataclass(frozen=True)
class BPTrace:
    c2v: np.ndarray
    v2c: np.ndarray
    llr: np.ndarray
    hard: np.ndarray
    valid: np.ndarray
    runtime: np.ndarray
    first_valid: int | None

def gf2_rank(matrix):
    work = np.asarray(matrix, dtype=np.uint8).copy() % 2
    rank = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        for row in np.flatnonzero(work[:, column]):
            if row != rank:
                work[row] ^= work[rank]
        rank += 1
        if rank == work.shape[0]:
            break
    return rank

def graph_adjacency(H):
    checks, variables = np.nonzero(H)
    checks, variables = (checks.astype(np.int32), variables.astype(np.int32))
    edge_ids = np.arange(len(checks), dtype=np.int32)
    check_edges = tuple((edge_ids[checks == check][::-1] for check in range(H.shape[0])))
    variable_edges = tuple((edge_ids[variables == variable][::-1] for variable in range(H.shape[1])))
    return (checks, variables, check_edges, variable_edges)

def context_for(level, p, *, CA_9_2_5, CONTEXTS, CircuitBuilder, CodeContext, ErrorRate, TriangularCode, detector_error_model_to_check_matrices, gf2_rank, graph_adjacency):
    key = (int(level), float(p))
    if key in CONTEXTS:
        return CONTEXTS[key]
    code = CA_9_2_5() if level == 2 else TriangularCode(level)
    circuit = CircuitBuilder(code).get_circuit_level_circuit(ErrorRate(p=float(p), p_meas=0), d_rounds=1)
    dem = circuit.detector_error_model()
    matrices = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
    H = matrices.check_matrix.toarray().astype(np.uint8) % 2
    L = matrices.observables_matrix.toarray().astype(np.uint8) % 2
    priors = np.asarray(matrices.priors, dtype=float)
    assert np.all((priors > 0) & (priors < 1)) and gf2_rank(H) == H.shape[0]
    checks, variables, check_edges, variable_edges = graph_adjacency(H)
    distance = 5 if level == 2 else int(code.d)
    assert distance in {5, 9}
    context = CodeContext(level, distance, float(p), H, L, priors, np.log((1 - priors) / priors), dem, checks, variables, check_edges, variable_edges, dem.get_detector_coordinates())
    CONTEXTS[key] = context
    return context

def transformed_context(context, transform, *, graph_adjacency):
    transform = np.asarray(transform, dtype=np.uint8) % 2
    H = transform @ context.H % 2
    assert H.shape[1] == context.H.shape[1] and (not np.any(H.sum(axis=1) == 0))
    checks, variables, check_edges, variable_edges = graph_adjacency(H)
    return replace(context, H=H, edge_checks=checks, edge_variables=variables, check_edges=check_edges, variable_edges=variable_edges)

def trace_bp(context, syndrome, iterations, update_rule='min_sum', scaling=1.0, offset=0.0, damping=0.0, schedule='parallel', rng=None, prior_llr=None, initial_message_sigma=0.0, *, BPTrace):
    trace = context.bp_decoder.trace(
        syndrome, iterations, update_rule=update_rule,
        scaling=scaling, offset=offset, damping=damping, schedule=schedule,
        rng=rng, prior_llr=prior_llr, initial_message_sigma=initial_message_sigma,
    )
    return BPTrace(trace.c2v, trace.v2c, trace.llr, trace.hard, trace.valid, trace.runtime_seconds, trace.first_valid)

def logical_label(context, correction, *, bits):
    return bits(context.L @ np.asarray(correction, dtype=np.uint8) % 2)

def outcome(context, correction, syndrome, truth, error=None, *, bits):
    correction = np.asarray(correction, dtype=np.uint8)
    valid = np.array_equal(context.H @ correction % 2, syndrome)
    sector = context.L @ correction % 2
    correct = valid and np.array_equal(sector, truth)
    exact = np.nan if error is None else int(np.array_equal(correction, error))
    return (int(valid), int(correct), exact, bits(sector))

def prior_cost(context, correction, prior_llr=None):
    weights = context.prior_llr if prior_llr is None else np.asarray(prior_llr)
    return float(np.asarray(correction, dtype=np.uint8) @ weights)

def full_osd_order(context, *, CONFIG):
    return int(CONFIG['configured_full_osd_order_by_distance'][str(context.distance)])

def osd_decode(context, syndrome, llr, order, *, ordered_statistics_postprocess):
    llr = np.ascontiguousarray(np.asarray(llr, dtype=np.float64))
    syndrome = np.ascontiguousarray(np.asarray(syndrome, dtype=np.uint8))
    started = perf_counter()
    correction = ordered_statistics_postprocess(llr, syndrome, np.ascontiguousarray(context.H), int(order), np.ascontiguousarray(context.priors, dtype=np.float64))
    runtime = perf_counter() - started
    assert np.array_equal(context.H @ correction % 2, syndrome)
    return (correction, runtime)

def standard_bposd(context, syndrome, trace, *, CONFIG, full_osd_order, osd_decode):
    valid_iterations = np.flatnonzero(trace.valid[1:CONFIG['bp_standard_iterations'] + 1]) + 1
    if len(valid_iterations):
        endpoint = int(valid_iterations[0])
        return (trace.hard[endpoint].copy(), endpoint, 0, float(trace.runtime[endpoint]))
    endpoint = CONFIG['bp_standard_iterations']
    correction, osd_runtime = osd_decode(context, syndrome, trace.llr[endpoint], full_osd_order(context))
    return (correction, endpoint, 1, float(trace.runtime[endpoint] + osd_runtime))

def trace_metrics(context, trace, endpoint=None):
    endpoint = len(trace.valid) - 1 if endpoint is None else int(endpoint)
    if endpoint == 0:
        return {'message_instability': 0.0, 'message_instability_max': 0.0, 'belief_instability': 0.0, 'sign_flip': 0.0, 'extrinsic_conflict': 0.0, 'posterior_abs_mean': float(np.mean(np.abs(trace.llr[0]))), 'reliability_rank': json.dumps(np.argsort(np.abs(trace.llr[0])).tolist())}
    start = max(1, endpoint - max(2, endpoint // 5) + 1)
    belief_move = np.mean(np.abs(np.diff(trace.llr[start - 1:endpoint + 1], axis=0)))
    sign_flip = np.mean(np.diff(trace.hard[start - 1:endpoint + 1], axis=0) != 0)
    message_delta = np.abs(np.diff(trace.c2v[start - 1:endpoint + 1], axis=0))
    message_move = np.mean(message_delta)
    conflict_values = []
    for variable, ids in enumerate(context.variable_edges):
        incoming = trace.c2v[start:endpoint + 1, ids]
        denom = np.abs(incoming).sum(axis=1)
        conflict_values.extend((1 - np.abs(incoming.sum(axis=1)) / (denom + 1e-12)).tolist())
    rank_order = np.argsort(np.abs(trace.llr[endpoint]))
    return {'message_instability': float(message_move), 'message_instability_max': float(np.max(message_delta)), 'belief_instability': float(belief_move), 'sign_flip': float(sign_flip), 'extrinsic_conflict': float(np.mean(conflict_values)), 'posterior_abs_mean': float(np.mean(np.abs(trace.llr[endpoint]))), 'reliability_rank': json.dumps(rank_order.tolist())}

def tesseract_predictions(context, syndromes, *, CONFIG, Tesseract):
    compile_started = perf_counter()
    decoder = Tesseract(beam=CONFIG['tesseract_beam'], pq_limit=CONFIG['tesseract_pq_limit'], seed=CONFIG['random_seed']).compile_decoder_for_dem(dem=context.dem)
    compile_runtime = perf_counter() - compile_started
    packed = np.packbits(np.asarray(syndromes, dtype=np.uint8), axis=1, bitorder='little')
    predictions = decoder.decode_shots_bit_packed(bit_packed_detection_event_data=packed)
    per_shot_runtime = np.full(len(syndromes), np.nan)
    timing_count = min(len(syndromes), CONFIG['tesseract_timing_shots_per_batch'])
    compile_amortized = compile_runtime / max(1, len(syndromes))
    for index in tqdm(range(timing_count), desc='Tesseract timing', leave=False, disable=CONFIG['print_level'] == 1):
        started = perf_counter()
        timed = decoder.decode_shots_bit_packed(bit_packed_detection_event_data=packed[index:index + 1])
        per_shot_runtime[index] = perf_counter() - started + compile_amortized
        assert np.array_equal(timed, predictions[index:index + 1])
    unpacked = np.unpackbits(predictions, axis=1, count=context.L.shape[0], bitorder='little').astype(np.uint8)
    return (unpacked, per_shot_runtime, compile_amortized)

def exact_sector_costs(context, syndrome, *, CONFIG, bits):
    n, m, k = (context.H.shape[1], context.H.shape[0], context.L.shape[0])
    A = np.zeros((m + k, n + m + k))
    A[:m, :n] = context.H
    A[:m, n:n + m] = -2 * np.eye(m)
    A[m:, :n] = context.L
    A[m:, n + m:] = -2 * np.eye(k)
    lower = np.zeros(n + m + k)
    upper = np.concatenate([np.ones(n), context.H.sum(axis=1) // 2, context.L.sum(axis=1) // 2])
    objective = np.concatenate([context.prior_llr, np.zeros(m + k)])
    results = []
    for sector_index in range(1 << k):
        sector = np.asarray([sector_index >> k - 1 - i & 1 for i in range(k)], dtype=float)
        rhs = np.concatenate([syndrome.astype(float), sector])
        started = perf_counter()
        result = milp(objective, integrality=np.ones(n + m + k), bounds=Bounds(lower, upper), constraints=LinearConstraint(csc_matrix(A), rhs, rhs), options={'time_limit': CONFIG['milp_time_limit_seconds']})
        assert result.status in {0, 1, 2}, result.message
        results.append({'sector': bits(sector), 'cost': np.nan if result.fun is None else float(result.fun), 'exact': int(result.status == 0), 'solver_status': int(result.status), 'runtime': perf_counter() - started})
    return results
