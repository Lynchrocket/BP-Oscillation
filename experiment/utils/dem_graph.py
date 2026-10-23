# failure_regions, trapping_sets
"""Detector-error-model parsing and graph construction helpers."""
from collections import Counter
import stim
import numpy as np
import networkx as nx
import pandas as pd


def extract_dem_info(dem: stim.DetectorErrorModel) -> tuple[pd.DataFrame, Counter, Counter]:
    """
    Parse only error(...) terms and extract structure tokens D*, L*

    Returns
    --------
    terms_df: DataFrame
        columns ['p_error', 'n_det', 'n_obs', 'det_ids', 'obs_ids'] for each error(...) term in the DEM.
    det_degree: Counter
        mapping detector_id -> count of error(...) terms it appears in. -> a det influenced by how many error terms
    obs_degree: Counter
        mapping observable_id -> count of error(...) terms it appears in. -> an obs influenced by how many error terms
    """
    term_rows = []
    det_degree = Counter()
    obs_degree = Counter()
    error_idx = 0
    for inst in dem:
        if inst.type != 'error':
            continue
        p_val = inst.args_copy()[0]
        targets = inst.targets_copy()
        det_ids = []
        obs_ids = []
        for t in targets:
            if t.is_relative_detector_id():
                di = t.val
                det_ids.append(di)
                det_degree[di] += 1
            elif t.is_logical_observable_id():
                li = t.val
                obs_ids.append(li)
                obs_degree[li] += 1
        term_rows.append({'error_idx': error_idx, 'p_error': float(p_val), 'n_det': len(det_ids), 'n_obs': len(obs_ids), 'det_ids': det_ids, 'obs_ids': obs_ids})
        error_idx += 1
    terms_df = pd.DataFrame(term_rows)
    return (terms_df, det_degree, obs_degree)

def extract_coord_from_graph(G: nx.Graph, dim: int=3) -> dict[int, tuple[float, float]]:
    """
    Extracts the 'coord' attribute from graph nodes and returns a mapping from node_id to (x, y, t) coordinates.

    If a node does not have a 'coord' attribute, it defaults to (-1, -1, -1).
    """
    default_coord = (-2, -2, -2)[:dim]
    coord = {}
    for node, data in G.nodes(data=True):
        coord[node] = tuple(data.get('coord', default_coord))[:dim]
    return coord

def build_tanner_graph_from_dem(dem: stim.DetectorErrorModel, *, extract_dem_info):
    terms_df, _, _ = extract_dem_info(dem)
    coords = dem.get_detector_coordinates()
    x_min, x_max, y_min, y_max = (min((c[0] for c in coords.values())), max((c[0] for c in coords.values())), min((c[1] for c in coords.values())), max((c[1] for c in coords.values())))
    x_center, y_center = ((x_min + x_max) / 2, (y_min + y_max) / 2)
    G = nx.Graph()
    detector_nodes = [(c[0], dict(coord=c[1])) for c in coords.items()]
    G.add_nodes_from(detector_nodes, kind='detector')
    for _, row in terms_df.iterrows():
        error_id = int(row['error_idx']) + dem.num_detectors
        det_ids = list(row['det_ids'])
        if len(det_ids) == 1:
            d = int(det_ids[0])
            coord = [coords[d][0] - 1, coords[d][1] + 0.5, *coords[d][2:]] if coords[d][0] < x_center else [coords[d][0] + 1, coords[d][1] + 0.5, *coords[d][2:]]
            G.add_node(error_id, coord=coord, kind='error')
            G.add_edge(error_id, d)
        elif len(det_ids) > 1:
            center_coord = np.mean([coords[int(d)] for d in det_ids], axis=0).tolist()
            G.add_node(error_id, coord=center_coord, kind='error')
            for d in det_ids:
                G.add_edge(error_id, int(d))
    return (G, terms_df)

def event_to_active_ids(detector_event_bits: str):
    """
    Converts a detector event (string of '0' and '1') into a list of active detector ids.

    '00101001' => [2, 4, 7]
    """
    return [i for i, ch in enumerate(detector_event_bits) if ch == '1']
