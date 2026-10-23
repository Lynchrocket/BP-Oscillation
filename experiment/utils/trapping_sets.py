# trapping_sets
"""Shared collection and output helpers for trapping-set structural studies."""
from pathlib import Path
from collections import defaultdict
from itertools import combinations

import numpy as np
import networkx as nx
from .misc import id2str, str2id, pack_bits_little, unpack_bits_little
from ast import literal_eval

import pandas as pd
from tqdm.auto import tqdm

from src.noises_ca import CircuitBuilder, ErrorRate
from experiment.utils.graph_struct import (
    build_tanner_graph_from_dem, err_subgraph_struct, det_subgraph_struct,
)

EXPANSION_COLUMNS = ['variables', 'checks', 'odd_checks', 'a', 'b', 'critical_number', 'strength',
                     'minimal_failures', 'seed_cycle', 'expanded_variables', 'expanded_checks']


def parse_node_ids(value: str):
    """Read saved integer node collections without evaluating an expression."""
    nodes = literal_eval(value)
    assert isinstance(nodes, (list, tuple, set)), nodes
    assert all(type(node) is int for node in nodes), nodes
    return nodes


def collect_low_weight_structure(code, cases: pd.DataFrame, *, weight: int, print_level: int):
    """Preserve row-aligned error- and detector-induced feature populations."""
    assert weight in {2, 3, 4}
    assert type(print_level) is int and print_level in {1, 2, 3}
    assert {'p', 'err', 'det', 'true', 'pred_b', 'pred_t'} <= set(cases.columns)
    builder = CircuitBuilder(code)
    error_features, detector_features = [], []
    for _, row in tqdm(cases.iterrows(), total=len(cases), desc=f"weight-{weight} structures",
                       disable=print_level == 1):
        assert isinstance(row['err'], str) and set(row['err']) <= {'0', '1'}
        assert row['err'].count('1') == weight
        circuit = builder.get_circuit_level_circuit(ErrorRate(p=float(row['p']), p_meas=0.0), d_rounds=1)
        graph, _ = build_tanner_graph_from_dem(circuit.detector_error_model())
        error_features.append(err_subgraph_struct(graph, row['err']))
        detector_features.append(det_subgraph_struct(graph, row['det']))
    assert len(error_features) == len(detector_features) == len(cases)
    return pd.DataFrame(error_features), pd.DataFrame(detector_features)


def low_weight_summary(cases: pd.DataFrame, features: pd.DataFrame, *, weight: int, population: str):
    """Aggregate structures separately for each noise value and outcome label."""
    assert population in {'error_induced', 'detector_induced'}
    assert len(cases) == len(features) and 'WL_hash' in features
    frame = cases[['p', 'true', 'pred_b', 'pred_t']].reset_index(drop=True)
    frame['label'] = ((frame['pred_b'] == frame['true']).astype(int).astype(str)
                      + (frame['pred_t'] == frame['true']).astype(int).astype(str))
    frame['WL_hash'] = features['WL_hash'].reset_index(drop=True)
    summary = frame.groupby(['p', 'label', 'WL_hash'], dropna=False).size().rename('case_count').reset_index()
    summary['population_count'] = summary.groupby(['p', 'label'])['case_count'].transform('sum')
    summary['fraction'] = summary['case_count'] / summary['population_count']
    summary.insert(0, 'population', population)
    summary.insert(0, 'weight', weight)
    assert summary['case_count'].sum() == len(cases)
    return summary


def save_plotly(figure, path: Path):
    """Persist a returned Plotly figure locally while retaining notebook display."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(path)
    return figure


def save_matplotlib(axes, path: Path):
    """Save a Pandas plot's Matplotlib figure and retain the returned axes."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    axes.get_figure().savefig(path, bbox_inches='tight')
    return axes


def save_expansion_records(records, *, trace_dir: Path, stats_dir: Path):
    """Keep constructive records separate from aggregate (a,b) statistics."""
    frame = pd.DataFrame([record.to_dict() for record in records], columns=EXPANSION_COLUMNS)
    summary = frame.groupby(['a', 'b', 'critical_number'], dropna=False).agg(
        record_count=('strength', 'size'), strength_sum=('strength', 'sum')).reset_index()
    assert summary['record_count'].sum() == len(records)
    trace_dir, stats_dir = Path(trace_dir), Path(stats_dir)
    trace_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(trace_dir / 'expansion_contraction_records.csv', index=False)
    summary.to_csv(stats_dir / 'expansion_contraction_summary.csv', index=False)
    return frame, summary


def expansion_probability_summary(frame: pd.DataFrame):
    """Preserve the existing per-probability contraction-summary fields."""
    assert {'p', 'variables', 'a', 'b', 'critical_number', 'strength'} <= set(frame.columns)
    return frame.groupby('p').agg(
        num_records=('p', 'size'), num_unique_ts=('variables', 'nunique'),
        mean_a=('a', 'mean'), mean_b=('b', 'mean'), min_critical=('critical_number', 'min'),
        mean_critical=('critical_number', 'mean'), max_strength=('strength', 'max'),
        mean_strength=('strength', 'mean')).reset_index()


def collect_expansion_sweep(code, p_values, *, expansion_contraction, decoder_fails,
                            depth: int, max_error_weight: int, trace_path: Path,
                            stats_path: Path, print_level: int):
    """Run the four existing sweeps with common output/schema handling."""
    assert len(p_values) and type(print_level) is int and print_level in {1, 2, 3}
    builder = CircuitBuilder(code)
    trace_path, stats_path = Path(trace_path), Path(stats_path)
    trace_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    columns = EXPANSION_COLUMNS + ['p']
    # An empty completed cohort must retain a parseable schema, not a blank file.
    pd.DataFrame(columns=columns).to_csv(trace_path, index=False)
    for p in tqdm(p_values, desc='Expansion/contraction probabilities', disable=print_level == 1):
        if print_level >= 2:
            print(f'Processing p = {p:.6f} ...')
        circuit = builder.get_circuit_level_circuit(ErrorRate(p=float(p), p_meas=0.0), d_rounds=1)
        dem = circuit.detector_error_model()
        records = expansion_contraction(
            dem, decoder_fails=decoder_fails, max_cycle_len=6, expand_depth=depth,
            max_error_weight=max_error_weight, max_cycles=5000)
        if print_level >= 2:
            print(f'Found {len(records)} trapping sets for p = {p:.6f}')
        frame = pd.DataFrame([{**r.to_dict(), 'p': float(p)} for r in records], columns=columns)
        frame.to_csv(trace_path, mode='a', header=False, index=False)
    table = pd.read_csv(trace_path)
    expansion_probability_summary(table).to_csv(stats_path, index=False)
    if print_level >= 2:
        print(f'saved: {trace_path}')
    return table, dem

def construct_fault_in_weight(candidate, num_vars, weight):
    faults = []
    for var_ids in combinations(candidate, weight):
        fault = np.zeros(num_vars, dtype=np.uint8)
        fault[list(var_ids)] = 1
        faults.append(fault)
    
    return np.array(faults)

def trapping_set_candidates(G: nx.Graph, a_max: int, b_max: int):
    """
    Finds candidate trapping sets in the graph G with parameters (a, b) s.t. a <= a_max and b <= b_max.
    
    (a, b) classical-type TS: a variable nodes, b odd-degree check nodes
    (a, 0) symmetric stabilizer TS: a variable nodes, 0 odd-degree check nodes, and a is even.

    Usually, we are interested in small trapping sets (small a and b) since they are more likely to cause decoding failures.
    
    Parameters
    ----------
    G: nx.Graph
        the tanner graph built from the DEM.
    a_max: int
        maximum number of variable nodes in the trapping set.
    b_max: int
        maximum number of odd-degree check nodes in the trapping set.
    
    Returns
    ----------
    ct_candidates: dict
        dictionary mapping (a, b) to list of candidate classical-type trapping sets, where each candidate is represented as a set of node ids (including only variable nodes).
    ss_candidates: dict
        dictionary mapping (a, b) to list of candidate symmetric stabilizer trapping sets, where each candidate is represented as a set of node ids (including only variable nodes).
    """
    # (a, b) -> list of candidate sets with those parameters
    ct_candidates = defaultdict(list)
    ss_candidates = defaultdict(list)
    error_nodes = [n for n, d in G.nodes(data=True) if d['kind'] == 'error']
    detector_nodes = [n for n, d in G.nodes(data=True) if d['kind'] == 'detector']
    num_dets = len(detector_nodes)
    a_bound = min(a_max, len(error_nodes))
    for a in range(1, a_bound + 1):
        for es in combinations(error_nodes, a):
            neighbor_checks = set()
            for e in es:
                neighbor_checks.update(G.neighbors(e))
            SG = G.subgraph(neighbor_checks.union(es))
            odd_degree_checks = [c for c in neighbor_checks if SG.degree(c) % 2 == 1]
            b = len(odd_degree_checks)
            if b <= b_max:
                ct_candidates[(a, b)].append(set(map(lambda x: x - num_dets, es)))
            if b == 0 and a % 2 == 0:
                ss_candidates[(a, b)].append(set(map(lambda x: x - num_dets, es)))

    return ct_candidates, ss_candidates

def candidates_to_fault_det(ct_candidates, ss_candidates, H_dem, num_vars, weight):
    """
    Converts the candidate trapping sets to their corresponding faults (and then detector events by applying the detector error matrix).

    Parameters
    ----------
    ct_candidates: dict
        dictionary mapping (a, b) to list of candidate classical-type trapping sets, where each candidate is represented as a set of node ids (including only variable nodes).
    ss_candidates: dict
        dictionary mapping (a, b) to list of candidate symmetric stabilizer trapping sets, where each candidate is represented as a set of node ids (including only variable nodes).
    H_dem: np.ndarray
        the detector error matrix.
    num_vars: int
        the number of variable nodes.
    weight: int
        the weight of faults to consider.

    Returns
    ----------
    ct_faults: dict
    ss_faults: dict
    ct_dets: dict
    ss_dets: dict
    """
    ct_faults, ss_faults = {}, {}
    ct_dets, ss_dets = {}, {}
    if len(ct_candidates) > 0:
        for (a, b), candidates in ct_candidates.items():
            if a < weight:
                continue
            faults = np.vstack([construct_fault_in_weight(candidate, num_vars, weight) for candidate in candidates], dtype=np.uint8)
            unique_faults = np.unique(faults, axis=0)
            ct_faults[(a, b)] = unique_faults
            dets = (unique_faults @ H_dem.T) % 2
            ct_dets[(a, b)] = dets
    
    if len(ss_candidates) > 0:
        for (a, b), candidates in ss_candidates.items():
            if a < weight:
                continue
            faults = np.vstack([construct_fault_in_weight(candidate, num_vars, weight) for candidate in candidates], dtype=np.uint8)
            unique_faults = np.unique(faults, axis=0)
            ss_faults[(a, b)] = unique_faults
            dets = (unique_faults @ H_dem.T) % 2
            ss_dets[(a, b)] = dets
    
    return ct_faults, ss_faults, ct_dets, ss_dets

def critical_number(ts):
    ...

def strength(ts):
    ...

def decode_on_candidates(dem, compiled_decoder, ct_candidates, ss_candidates, num_vars, save_folder: Path):
    from beliefmatching import detector_error_model_to_check_matrices

    demMatrices = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
    H_dem = demMatrices.check_matrix
    L_dem = demMatrices.observables_matrix
    num_obs = L_dem.shape[0]

    overview = []
    failure_cases = []
    max_weight = max(max((a for (a, b) in ct_candidates.keys()), default=0), max((a for (a, b) in ss_candidates.keys()), default=0))
    for weight in range(1, max_weight + 1):
        print(f'Analyzing candidates with fault weight {weight}...')
        ct_faults, ss_faults, ct_dets, ss_dets = candidates_to_fault_det(ct_candidates, ss_candidates, H_dem, num_vars, weight)

        for (a, b), dets in ct_dets.items():
            faults = ct_faults[(a, b)]
            obs = (faults @ L_dem.T) % 2
            dets_packed = pack_bits_little(dets)
            pred_packed = compiled_decoder.decode_shots_bit_packed(bit_packed_detection_event_data=dets_packed)
            pred = unpack_bits_little(pred_packed, num_obs)
            succ = np.all(pred == obs, axis=1)
            print(f'    Classical-type TS (a={a}, b={b}): total {len(succ)}, {np.sum(succ)} succeeded, {len(succ) - np.sum(succ)} failed, success rate {np.sum(succ) / len(succ) * 100:.2f}%')
            overview.append({
                'type': 'classical-type',
                'a': a,
                'b': b,
                'weight': weight,
                'total': len(succ),
                'succeeded': np.sum(succ),
                'failed': len(succ) - np.sum(succ),
            })
            failure_cases.extend([{
                'type': 'classical-type',
                'a': a,
                'b': b,
                'weight': weight,
                'fault': ''.join(map(str, faults[i])),
                'det': ''.join(map(str, dets[i])),
                'obs': ''.join(map(str, obs[i])),
                'pred': ''.join(map(str, pred[i])),
            } for i in range(len(succ)) if not succ[i]])

        print('---------------------------------------------------------------------------')

        for (a, b), dets in ss_dets.items():
            faults = ss_faults[(a, b)]
            obs = (faults @ L_dem.T) % 2
            dets_packed = pack_bits_little(dets)
            pred_packed = compiled_decoder.decode_shots_bit_packed(bit_packed_detection_event_data=dets_packed)
            pred = unpack_bits_little(pred_packed, num_obs)
            succ = np.all(pred == obs, axis=1)
            print(f'    Symmetric stabilizer TS (a={a}, b={b}): total {len(succ)}, {np.sum(succ)} succeeded, {len(succ) - np.sum(succ)} failed, success rate {np.sum(succ) / len(succ) * 100:.2f}%')
            overview.append({
                'type': 'symmetric stabilizer',
                'a': a,
                'b': b,
                'weight': weight,
                'total': len(succ),
                'succeeded': np.sum(succ),
                'failed': len(succ) - np.sum(succ),
            })
            failure_cases.extend([{
                'type': 'symmetric stabilizer',
                'a': a,
                'b': b,
                'weight': weight,
                'fault': ''.join(map(str, faults[i])),
                'det': ''.join(map(str, dets[i])),
                'obs': ''.join(map(str, obs[i])),
                'pred': ''.join(map(str, pred[i])),
            } for i in range(len(succ)) if not succ[i]])
        print('===========================================================================')

    overview, failure_cases = pd.DataFrame(overview), pd.DataFrame(failure_cases)

    overview.to_csv(save_folder / 'TS_overview.csv', index=False)
    failure_cases.to_feather(save_folder / 'TS_failure_cases.feather')

    return overview, failure_cases
