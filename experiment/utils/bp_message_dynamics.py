# bp_failure_dynamics, instability_guided_bp, quotient, trapping_sets
"""Reusable experiment functions with explicit arguments and dependencies."""
from .bp_trace import BPContext
from dataclasses import dataclass
import numpy as np
import pandas as pd
import stim

@dataclass(frozen=True)
class Context(BPContext):
    p: float
    H: np.ndarray
    L: np.ndarray
    priors: np.ndarray
    prior_llr: np.ndarray
    edge_checks: np.ndarray
    edge_variables: np.ndarray
    check_edges: tuple[np.ndarray, ...]
    variable_edges: tuple[np.ndarray, ...]

@dataclass(frozen=True)
class MessageTrace:
    c2v: np.ndarray
    v2c: np.ndarray
    llr: np.ndarray
    hard: np.ndarray
    syndrome_consistent: np.ndarray
    first_valid_iteration: int | None

def logical(context, correction):
    return (context.L @ correction % 2).astype(np.uint8)

def context_for(p, *, CACHE, CA_9_2_5, CircuitBuilder, Context, ErrorRate, detector_error_model_to_check_matrices):
    p = float(p)
    if p in CACHE:
        return CACHE[p]
    circuit = CircuitBuilder(CA_9_2_5()).get_circuit_level_circuit(ErrorRate(p=p, p_meas=0), d_rounds=1)
    matrices = detector_error_model_to_check_matrices(circuit.detector_error_model(), allow_undecomposed_hyperedges=True)
    H = matrices.check_matrix.toarray().astype(np.uint8) % 2
    L = matrices.observables_matrix.toarray().astype(np.uint8) % 2
    priors = np.asarray(matrices.priors, float)
    assert H.shape[1] == priors.size and np.all((priors > 0) & (priors < 1))
    checks, variables = np.nonzero(H)
    checks = checks.astype(np.int32)
    variables = variables.astype(np.int32)
    edge_ids = np.arange(checks.size, dtype=np.int32)
    check_edges = tuple((edge_ids[checks == c] for c in range(H.shape[0])))
    variable_edges = tuple((edge_ids[variables == v] for v in range(H.shape[1])))
    assert all((ids.size >= 2 for ids in check_edges))
    result = Context(p, H, L, priors, np.log((1 - priors) / priors), checks, variables, check_edges, variable_edges)
    CACHE[p] = result
    return result

def trace_min_sum(context, syndrome, iterations, stop_early=False, *, CONFIG, MessageTrace, rng=None, initial_message_sigma=0.0):
    trace = context.bp_decoder.trace(
        syndrome, iterations, scaling=CONFIG['ms_scaling_factor'],
        stop_early=stop_early, accumulation='sum_then_prior',
        rng=rng, initial_message_sigma=initial_message_sigma,
    )
    return MessageTrace(trace.c2v, trace.v2c, trace.llr, trace.hard, trace.valid, trace.first_valid)

def edge_metrics(context, trace, *, CONFIG):
    messages = trace.c2v[1:]
    delta = np.abs(np.diff(messages, axis=0))
    assert delta.shape[0] >= 2
    late_n = max(1, int(np.ceil(CONFIG['late_fraction'] * delta.shape[0])))
    late = delta[-late_n:]
    tol = CONFIG['message_zero_tolerance']
    signs = np.where(np.abs(messages) <= tol, 0, np.where(messages > 0, 1, -1))
    comparable = (signs[1:] != 0) & (signs[:-1] != 0)
    flips = (signs[1:] != signs[:-1]) & comparable
    d1 = delta.mean(0)
    d2 = np.abs(messages[2:] - messages[:-2]).mean(0)
    prior_scale = np.abs(context.prior_llr[context.edge_variables]) + np.finfo(float).eps
    return pd.DataFrame({'edge_id': np.arange(context.edge_checks.size), 'variable_id': context.edge_variables, 'late_residual': late.mean(0), 'normalized_late_residual': late.mean(0) / prior_scale, 'sign_flip_rate': flips.sum(0) / max(1, flips.shape[0]), 'period2_strength': d1 * d1 / (d1 + d2 + np.finfo(float).eps)})

def variable_metrics(context, trace, edges, *, CONFIG):
    incoming = trace.c2v[1:]
    beliefs = trace.llr[1:]
    movement = np.abs(np.diff(beliefs, axis=0))
    late_n = max(1, int(np.ceil(CONFIG['late_fraction'] * movement.shape[0])))
    rows = []
    for variable, ids in enumerate(context.variable_edges):
        messages = incoming[:, ids]
        ext_abs = np.abs(messages).sum(1)
        conflict = 1 - np.abs(messages.sum(1)) / (ext_abs + np.finfo(float).eps)
        rows.append({'belief_instability': movement[-late_n:, variable].mean(), 'message_instability': edges.loc[ids, 'late_residual'].mean(), 'extrinsic_conflict': conflict[-late_n:].mean()})
    return pd.DataFrame(rows)

def trace_prefix(trace, iterations, *, MessageTrace):
    assert 0 < iterations < len(trace.llr)
    end = iterations + 1
    valid_indices = np.flatnonzero(trace.syndrome_consistent[:end])
    first = None if valid_indices.size == 0 else int(valid_indices[0])
    return MessageTrace(trace.c2v[:end], trace.v2c[:end], trace.llr[:end], trace.hard[:end], trace.syndrome_consistent[:end], first)

def benchmark_state(trace, *, CONFIG):
    index = trace.first_valid_iteration if trace.first_valid_iteration is not None else CONFIG['benchmark_iterations']
    assert index <= CONFIG['benchmark_iterations'] and index < len(trace.hard)
    return (index, trace.hard[index], bool(trace.syndrome_consistent[index]))

def classify_outcome(syndrome_consistent, logical_failure):
    labels = {(False, True): 'invalid-logical-mismatch', (False, False): 'invalid-logical-match', (True, False): 'valid-correct', (True, True): 'valid-logical-wrong'}
    return labels[bool(syndrome_consistent), bool(logical_failure)]

def normalized_late_residual(context, trace, endpoint, *, CONFIG):
    messages = trace.c2v[1:endpoint + 1]
    delta = np.abs(np.diff(messages, axis=0))
    late_n = max(1, int(np.ceil(CONFIG['late_fraction'] * len(delta))))
    scale = np.abs(context.prior_llr[context.edge_variables]) + np.finfo(float).eps
    return float(np.max(delta[-late_n:] / scale))

@dataclass
class Matrices:
    H: np.ndarray
    L: np.ndarray
    p: np.ndarray

    def to_dem(self) -> stim.DetectorErrorModel:
        num_detectors, num_errors = self.H.shape
        num_observables, _ = self.L.shape
        dem = stim.DetectorErrorModel()
        for d in range(num_detectors):
            dem.append('detector', [], [stim.DemTarget.relative_detector_id(int(d))])
        for ell in range(num_observables):
            dem.append('logical_observable', [], [stim.DemTarget.logical_observable_id(int(ell))])
        for err in range(num_errors):
            targets = []
            for d in np.flatnonzero(self.H[:, err]):
                targets.append(stim.DemTarget.relative_detector_id(int(d)))
            for ell in np.flatnonzero(self.L[:, err]):
                targets.append(stim.DemTarget.logical_observable_id(int(ell)))
            if len(targets) == 0:
                continue
            dem.append('error', float(self.p[err]), targets)
        return dem

def message_edge_metrics(context, trace, direction, *, CONFIG):
    assert direction in {'c2v', 'v2c'}
    x = getattr(trace, direction)[1:]
    delta = np.abs(np.diff(x, axis=0))
    late_n = max(1, int(np.ceil(CONFIG['late_fraction'] * delta.shape[0])))
    late = delta[-late_n:]
    tol = CONFIG['message_zero_tolerance']
    signs = np.where(np.abs(x) <= tol, 0, np.where(x > 0, 1, -1))
    comparable = (signs[1:] != 0) & (signs[:-1] != 0)
    flips = (signs[1:] != signs[:-1]) & comparable
    d1 = delta.mean(0)
    d2 = np.abs(x[2:] - x[:-2]).mean(0)
    return pd.DataFrame({'edge_id': np.arange(context.edge_checks.size), 'check_id': context.edge_checks, 'variable_id': context.edge_variables, 'direction': direction, 'mean_residual': d1, 'max_residual': delta.max(0), 'late_residual': late.mean(0), 'normalized_late_residual': late.mean(0) / (np.abs(context.prior_llr[context.edge_variables]) + np.finfo(float).eps), 'sign_flip_rate': flips.sum(0) / max(1, flips.shape[0]), 'period2_strength': d1 * d1 / (d1 + d2 + np.finfo(float).eps)})

def message_variable_metrics(context, trace, edges, *, CONFIG):
    incoming = trace.c2v[1:]
    beliefs = trace.llr[1:]
    movement = np.abs(np.diff(beliefs, axis=0))
    late_n = max(1, int(np.ceil(CONFIG['late_fraction'] * movement.shape[0])))
    rows = []
    for v, ids in enumerate(context.variable_edges):
        msgs = incoming[:, ids]
        ext_abs = np.abs(msgs).sum(1)
        cancellation = 1 - np.abs(beliefs[:, v]) / (np.abs(context.prior_llr[v]) + ext_abs + np.finfo(float).eps)
        conflict = 1 - np.abs(msgs.sum(1)) / (ext_abs + np.finfo(float).eps)
        belief_instability = movement[-late_n:, v].mean()
        message_instability = edges.loc[ids, 'late_residual'].mean()
        rows.append({'variable_id': v, 'belief_instability': belief_instability, 'normalized_belief_instability': belief_instability / (abs(context.prior_llr[v]) + np.finfo(float).eps), 'message_instability': message_instability, 'normalized_message_instability': edges.loc[ids, 'normalized_late_residual'].mean(), 'late_cancellation': cancellation[-late_n:].mean(), 'late_extrinsic_conflict': conflict[-late_n:].mean()})
    return pd.DataFrame(rows)
