# failure_regions, force_osd_failures, low_weight_failures, quotient, trapping_sets
import stim
import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter
from itertools import combinations

from .misc import str2id

# =====================================
# Errors
# =====================================

def error_cross_time(err_times: list[float], time_window: float = 0.0) -> bool:
    """
    Checks if all error times are within a certain time window, which indicates that they are likely to be part of the same error event.
    
    Parameters
    ----------
    err_times: list[float]
        List of error times corresponding to the error events in the DEM.
    time_window: float
        The maximum allowed time difference between error events for them to be considered as part of the same event.
    
    Returns
    -------
    bool
        True if all error times are across the specified time window, False otherwise.
    """
    if len(err_times) <= 1:
        return True
    min_time, max_time = min(err_times), max(err_times)
    return (max_time - min_time) > time_window

def get_error_time(dem: stim.DetectorErrorModel, error_event: str) -> list[float]:
    errs = np.array([inst for inst in dem.flattened() if inst.type == "error"])
    err_ids = str2id(error_event)
    selected_errs = errs[err_ids]
    
    dets_coord = dem.get_detector_coordinates()
    err_times = []
    for err in selected_errs:
        targets = err.targets_copy()
        time_coords = [dets_coord[t.val][-1] for t in targets if t.is_relative_detector_id()]
        err_times.append(np.mean(time_coords))
    return err_times

def error_llr_weight(dem: stim.DetectorErrorModel, error_event: str) -> float:
    err_probs = np.array([inst.args_copy()[0] for inst in dem.flattened() if inst.type == "error"])
    err_ids = str2id(error_event)
    probs = err_probs[err_ids]
    llr_weight = sum(np.log((1 - probs) / probs))
    return llr_weight

def error_diameter(G: nx.Graph, error_event: str) -> float:
    """
    Computes the diameter of the errors, i.e, the maximum shortest path length between any two error nodes
    """
    num_dets = sum(1 for n, data in G.nodes(data=True) if data.get("kind") == "detector")
    err_ids = str2id(error_event, offset=num_dets)
    error_nodes = [n for n, data in G.nodes(data=True) if data.get("kind") == "error" and n in err_ids]

    dist = [nx.shortest_path_length(G, u, v) for u, v in combinations(error_nodes, 2)]
    max_dist = max(dist) if len(dist) > 0 else 0
    return max_dist

def dem_error_to_explanation(tag_circ: stim.Circuit, dem: stim.DetectorErrorModel) -> list[stim.ExplainedError]:
    """
    Converts a DEM error event to the real circuit errors.

    Parameters:
    ------------
    tag_circ: stim.Circuit
        The circuit with tags.
    dem: stim.DetectorErrorModel
        The detector error model containing the error events to be explained. Note that this is the merged DEM from circuit instead of tag_circ.
    """
    dem = dem.flattened()
    explanations = []
    for inst in dem:
        if inst.type != 'error':
            continue
        dem_filter = stim.DetectorErrorModel()
        dem_filter.append(inst)
        exp = tag_circ.explain_detector_error_model_errors(
            dem_filter=dem_filter,
        )
        explanations.append(exp[0])
    return explanations

def most_occurring_circ_error(explanations: list[stim.ExplainedError], err_id: int) -> str:
    """
    Given a mapping from DEM error indices to lists of explanations (circuit error locations), 
    returns the most frequently occurring circuit error type.
    """
    exp = explanations[err_id]
    from collections import Counter
    exps = [eloc.noise_tag for eloc in exp.circuit_error_locations]
    explanation_counter = Counter(exps)
    most_common_explanation, _ = explanation_counter.most_common(1)[0]
    return most_common_explanation

# =====================================
# Graph construction
# =====================================

def _extract_errs(dem: stim.DetectorErrorModel) -> list[dict]:
    """
    Parse only error(...) terms and extract structure tokens D*, L*
    
    Returns
    --------
    terms: list[dict]
        [{'p'=..., 'n_det'=..., 'n_obs'=..., 'det_ids'=..., 'obs_ids'=...}] for each error(...) term in the DEM.
    """
    terms = []
    error_idx = 0
    for inst in dem.flattened():
        if inst.type != 'error':
            continue
        
        p = inst.args_copy()[0] # error rate inside error(...)
        targets = inst.targets_copy()
        
        det_ids, obs_ids = [], []
        for t in targets:
            if t.is_relative_detector_id():
                det_ids.append(t.val)
            elif t.is_logical_observable_id():
                obs_ids.append(t.val)

        terms.append({
            'error_idx': error_idx,
            'p': float(p),
            'n_det': len(det_ids),
            'n_obs': len(obs_ids),
            'det_ids': det_ids,
            'obs_ids': obs_ids,
        })
        error_idx += 1
        
    return terms

def build_matching_graph_from_dem(dem: stim.DetectorErrorModel) -> tuple[nx.Graph, list[dict]]:
    """
    Builds a matching graph from the Detector Error Model (DEM).

    Returns
    ----------
        G: networkx.Graph 
            nodes -> detectors & edges -> error mechanism. 
            Edge attributes:
            * weight: minimum log((1-p)/p) across shared terms.
        terms: list[dict]
            [{'p'=..., 'n_det'=..., 'n_obs'=..., 'det_ids'=..., 'obs_ids'=...}] for each error(...) term in the DEM.
    """
    terms = _extract_errs(dem)

    coords = dem.get_detector_coordinates()
    x_min, x_max, y_min, y_max = min(c[0] for c in coords.values()), max(c[0] for c in coords.values()), min(c[1] for c in coords.values()), max(c[1] for c in coords.values())
    x_center, y_center = (x_min + x_max) / 2, (y_min + y_max) / 2

    G = nx.Graph()
    detector_nodes = [(c[0], dict(coord=c[1])) for c in coords.items()]
    G.add_nodes_from(detector_nodes)
    # virtual nodes for boundary (single-detector) terms
    G.graph['virtual_boundary_nodes'] = []
    # virtual nodes for hyperedges (multi-detector) terms
    G.graph['virtual_hyperedge_nodes'] = []
    vnode_id = 0
    hyperedge_id = 0

    for term in terms:
        det_ids = list(term['det_ids'])
        p = term['p']
        
        if len(det_ids) == 1: # boundary
            d = int(det_ids[0])
            vnode = f'B{vnode_id}'
            vnode_id += 1
            coord = [coords[d][0]-1, coords[d][1]+0.5, *coords[d][2:]] if coords[d][0] < x_center else [coords[d][0]+1, coords[d][1]+0.5, *coords[d][2:]]
            G.add_node(vnode, coord=coord) # virtual node for this boundary term, positioned at the near coordinate of the involved detector
            G.add_edge(d, vnode, p=p, is_boundary=True)
            G.graph['virtual_boundary_nodes'].append(vnode)
        elif len(det_ids) == 2:
            a, b = int(det_ids[0]), int(det_ids[1])
            if G.has_edge(a, b):
                G[a][b]['p'] = G[a][b]['p']*(1-p) + (1-G[a][b]['p'])*p
            else:
                G.add_edge(a, b, p=p)
        elif len(det_ids) > 2:
            # add a virtual node for this hyperedge and connect to all involved detectors
            center_coord = np.mean([coords[d] for d in det_ids], axis=0).tolist() # position the hyperedge node at the centroid of involved detectors
            vnode = f'H{hyperedge_id}'
            hyperedge_id += 1
            G.add_node(vnode, coord=center_coord)
            for d in det_ids:
                G.add_edge(vnode, d, p=p, is_hyperedge=True)
            G.graph['virtual_hyperedge_nodes'].append(vnode)

    for u, v, data in G.edges(data=True):
        p = data['p']
        data['weight'] = np.log((1 - p) / p)

    return G, terms

def build_tanner_graph_from_dem(dem: stim.DetectorErrorModel):
    terms = _extract_errs(dem)

    coords = dem.get_detector_coordinates()
    has_coords = any(len(v) > 0 for k, v in coords.items())
    if has_coords:
        x_min, x_max, y_min, y_max = min(c[0] for c in coords.values()), max(c[0] for c in coords.values()), min(c[1] for c in coords.values()), max(c[1] for c in coords.values())
        x_center, y_center = (x_min + x_max) / 2, (y_min + y_max) / 2
    
    G = nx.Graph()
    detector_nodes = [(c[0], dict(coord=c[1])) for c in coords.items()]
    G.add_nodes_from(detector_nodes, kind='detector')

    # edges between detectors and errors
    for term in terms:
        err_id = int(term['error_idx']) + dem.num_detectors # error nodes start after detector nodes
        det_ids = list(term['det_ids'])
        p = term['p']

        if len(det_ids) == 1:
            d = int(det_ids[0])
            if has_coords:
                coord = [coords[d][0]-1, coords[d][1]+0.5, *coords[d][2:]] if coords[d][0] < x_center else [coords[d][0]+1, coords[d][1]+0.5, *coords[d][2:]]
                G.add_node(err_id, coord=coord, kind='error', p=p)
            else:
                G.add_node(err_id, kind='error', p=p)
            G.add_edge(err_id, d)
        elif len(det_ids) > 1:
            if has_coords:
                center_coord = np.mean([coords[d] for d in det_ids], axis=0).tolist()
                G.add_node(err_id, coord=center_coord, kind='error', p=p)
            else:
                G.add_node(err_id, kind='error', p=p)
            for d in det_ids:
                G.add_edge(err_id, d)

    if not has_coords:
        coords = nx.spring_layout(G, dim=3, seed=42)
        nx.set_node_attributes(G, coords, 'coord')

    return G, terms

def build_tanner_graph_from_pcm(H: np.ndarray, stab_coords: list[tuple[float, float]], qubit_coords: list[tuple[float, float]]) -> nx.Graph:
    G = nx.Graph()
    num_dets, num_errs = H.shape
    for d in range(num_dets):
        G.add_node(d, coord=stab_coords[d], kind='detector')
    for e in range(num_errs):
        err_id = num_dets + e
        G.add_node(err_id, coord=qubit_coords[e], kind='error')
        for d in range(num_dets):
            if H[d, e] == 1:
                G.add_edge(err_id, d)
    return G

def convert_tanner_to_matrix(G: nx.Graph) -> tuple[np.ndarray, np.ndarray]:
    """
    Converts the tanner graph into a binary parity-check matrix H and a coordinate matrix C.

    Returns
    --------
    H: np.ndarray
        binary parity-check matrix where rows correspond to error nodes and columns correspond to detector nodes. H[i, j] = 1 if error node i is connected to detector node j.
    C: np.ndarray
        coordinate matrix where rows correspond to nodes in G and columns correspond to (x, y, t) coordinates. The order of rows follows the node ordering in G.nodes().
    """
    det_nodes = [n for n, d in G.nodes(data=True) if d.get('kind') == 'detector']
    error_nodes = [n for n, d in G.nodes(data=True) if d.get('kind') == 'error']
    
    det_index = {d: i for i, d in enumerate(det_nodes)}
    error_index = {e: i for i, e in enumerate(error_nodes)}

    H = np.zeros((len(error_nodes), len(det_nodes)), dtype=int)
    for e in error_nodes:
        for d in G.neighbors(e):
            if d in det_index:
                H[error_index[e], det_index[d]] = 1

    C = np.zeros((G.number_of_nodes(), 3), dtype=float)
    for i, (n, data) in enumerate(G.nodes(data=True)):
        C[i] = data.get('coord', (-1, -1, -1))[:3]

    return H, C

def extract_coord_from_graph(G: nx.Graph, dim: int = 3) -> dict[int, tuple[float, float]]:
    """
    Extracts the 'coord' attribute from graph nodes and returns a mapping from node_id to (x, y, t) coordinates.

    If a node does not have a 'coord' attribute, it defaults to (-1, -1, -1).
    """
    default_coord = (-2, -2, -2)
    coord = {}
    for node, data in G.nodes(data=True):
        coord[node] = tuple(data.get('coord', default_coord))[:dim]
    return coord

def det_events_to_coords(dem: stim.DetectorErrorModel, det_events: str) -> list[tuple[float, float, float]]:
    """
    Converts a string of detector events (e.g., '00101001') into a list of coordinates corresponding to the detectors coordinates.

    Returns
    --------
    coords: list[tuple[float, float, float]]
        In the format of [(x, y, t) ...]
    """
    coords = dem.get_detector_coordinates()
    det_ids = str2id(det_events)
    return [coords[i] for i in det_ids]

# =====================================
# Whole Graph Structures
# =====================================

def bipartite_graph_density(G: nx.Graph) -> float:
    num_dets = sum(1 for n, d in G.nodes(data=True) if d.get('kind') == 'detector')
    num_errs = sum(1 for n, d in G.nodes(data=True) if d.get('kind') == 'error')
    density = G.number_of_edges() / (num_dets * num_errs) if num_dets > 0 and num_errs > 0 else 0.0
    return density

def graph_struct(G: nx.Graph):
    """
    Extracts graph features from the whole Tanner graph.
    
    Parameters
    ----------
    G: nx.Graph
        the Tanner graph built from the DEM.
    
    Returns
    ----------
    features: dict
        * cycle_rank             (int): the number of independent cycles in the graph (Cyclomatic number, first Betti number), defined as |E| - |V| + |C|
        * inverse_girth        (float): the inverse of the length of the shortest cycle in the graph (0 if no cycle)
        * density              (float): density of the graph, defined as #edges over max possible #edges: |E|/(|V_e|*|V_d|)
        * num_boundary_err       (int): the number of error nodes that are connected to only one detector node
        * num_graphlike_err      (int): the number of error nodes that are connected to exactly two detector nodes
        * num_hyperedge_err      (int): the number of error nodes that are connected to more than two detector nodes
        * err_degree_dist    (list[int]): the error node degree distribution of the graph, i.e., a list of degrees for each error node in the graph
        * det_degree_dist    (list[int]): the detector node degree distribution of the graph, i.e., a list of degrees for each detector node in the graph
    """
    cycle_rank = G.number_of_edges() - G.number_of_nodes() + nx.number_connected_components(G)
    girth = nx.girth(G)

    density = bipartite_graph_density(G)
    
    errs = [n for n, d in G.nodes(data=True) if d.get('kind') == 'error']
    err_degree_dist = [d for n, d in G.degree() if G.nodes[n].get('kind') == 'error']
    det_degree_dist = [d for n, d in G.degree() if G.nodes[n].get('kind') == 'detector']

    return {
        'cycle_rank': cycle_rank,
        'inverse_girth': 1 / girth if np.isfinite(girth) else 0.0,
        'density': density,
        'num_boundary_err': sum(1 for n in errs if G.degree(n) == 1),
        'num_graphlike_err': sum(1 for n in errs if G.degree(n) == 2),
        'num_hyperedge_err': sum(1 for n in errs if G.degree(n) > 2),
        'err_degree_dist': err_degree_dist,
        'det_degree_dist': det_degree_dist,
    }

# =====================================
# Subgraph Structures
# =====================================

def local_incident_subgraph(G: nx.Graph, nodes: set[int]) -> nx.Graph:
    """
    Returns the subgraph of G induced by the given nodes and their incident edges.
    """
    affected_edges = [(u, v) for u, v in G.edges if u in nodes or v in nodes]
    SG = G.edge_subgraph(affected_edges).copy()
    missing_nodes = set(nodes) - set(SG.nodes)
    for n in missing_nodes:
        SG.add_node(n, **G.nodes[n])
    return SG

def err_subgraph_struct(G: nx.Graph, err_event: str):
    """
    Extract features from the error-support induced local subgraph.

    Given true error set E, construct:
        G[E cup N(E)]

    This subgraph asks whether the true error itself forms a local structure that may confuse BP.

    Parameters
    ----------
    G: nx.Graph
        the whole graph
    err_event: str
        string of '0' and '1' representing error event.

    Returns
    ----------
    features: dict
        * WL_hash                   (str): Weisfeiler-Lehman graph hash of the subgraph
        * cycle_rank                (int): the number of independent cycles in the subgraph, defined as |E| - |V| + |C|
        * inverse_girth           (float): the inverse of the length of the shortest cycle in the subgraph (0 if acyclic)
        * density                 (float): density of the graph, defined as #edges over max possible #edges: |E|/(|V_e|*|V_d|)
        * avg_degree              (float): average degree of the nodes in the subgraph
        * num_components            (int): number of connected components in the subgraph
        * largest_component_size    (int): size of the largest connected component in the subgraph
        * avg_comp_dets             (float): average number of detector nodes in each connected component of the subgraph
        * avg_comp_errs             (float): average number of error nodes in each connected component of the subgraph
        # * num_comp_dets_dist  (list[int]): the detector node degree distribution of the subgraph, i.e., a list of degrees for each detector node in the subgraph
        # * num_comp_errs_dist  (list[int]): the error node degree distribution of the subgraph, i.e., a list of degrees for each error node in the subgraph
        * num_annotated_errs        (int): number of annotated errors in the subgraph
    """
    dets = [n for n, d in G.nodes(data=True) if d.get('kind') == 'detector']
    errs = [n for n, d in G.nodes(data=True) if d.get('kind') == 'error']

    if len(errs) == 0:
        return {
            'WL_hash': '',
            'cycle_rank': 0,
            'inverse_girth': 0.0,
            'density': 0.0,
            'avg_degree': 0.0,
            'num_components': 0,
            'largest_component_size': 0,
            'avg_comp_dets': 0.0,
            'avg_comp_errs': 0.0,
            'num_annotated_errs': 0,
        }

    annotated_errs = [e for e, b in zip(errs, err_event) if b == '1']
    SG = local_incident_subgraph(G, annotated_errs)

    cycle_rank = SG.number_of_edges() - SG.number_of_nodes() + nx.number_connected_components(SG)
    girth = nx.girth(SG)

    density = bipartite_graph_density(SG)
    avg_degree = np.mean([d for n, d in SG.degree()]) if SG.number_of_nodes() > 0 else 0.0

    comps = list(nx.connected_components(SG))
    largest_comp_size = max(len(c) for c in comps) if comps else 0
    n_comp = len(comps)

    num_comp_dets_dist = [sum(1 for n in c if G.nodes[n].get('kind') == 'detector') for c in comps]
    num_comp_errs_dist = [sum(1 for n in c if G.nodes[n].get('kind') == 'error') for c in comps]

    return {
        'WL_hash': nx.weisfeiler_lehman_graph_hash(SG, node_attr='kind'),
        # graph features
        'cycle_rank': cycle_rank,
        'inverse_girth': 1 / girth if np.isfinite(girth) else 0.0,
        'density': density,
        'avg_degree': avg_degree,
        # component features
        'num_components': n_comp,
        'largest_component_size': largest_comp_size,
        'avg_comp_dets': np.mean(num_comp_dets_dist) if len(num_comp_dets_dist) > 0 else 0.0,
        'avg_comp_errs': np.mean(num_comp_errs_dist) if len(num_comp_errs_dist) > 0 else 0.0,
        # 'num_comp_dets_dist': num_comp_dets_dist,
        # 'num_comp_errs_dist': num_comp_errs_dist,
        # error-induced features
        'num_annotated_errs': len(annotated_errs),
    }

def det_subgraph_struct(G: nx.Graph, det_event: str):
    """
    Extracts graph features from the subgraph.
    
    Parameters
    ----------
    G: nx.Graph
        the whole graph
    det_event: str
        string of '0' and '1' representing detection event.

    Returns
    ----------
    features: dict
        * WL_hash                            (str): Weisfeiler-Lehman graph hash of the subgraph
        * cycle_rank                         (int): the number of independent cycles in the subgraph, defined as |E| - |V| + |C|
        * inverse_girth                    (float): the inverse of the length of the shortest cycle in the subgraph (0 if acyclic)
        * density                          (float): density of the graph, defined as #edges over max possible #edges: |E|/(|V_e|*|V_d|)
        * avg_degree                       (float): average degree of the nodes in the subgraph
        * num_components                     (int): number of connected components in the subgraph
        * largest_component_size             (int): size of the largest connected component in the subgraph
        * avg_comp_dets                      (float): average number of detector nodes in each connected component of the subgraph
        * avg_comp_errs                      (float): average number of error nodes in each connected component of the subgraph
        # * num_comp_dets_dist           (list[int]): the detector node degree distribution of the subgraph, i.e., a list of degrees for each detector node in the subgraph
        # * num_comp_errs_dist           (list[int]): the error node degree distribution of the subgraph, i.e., a list of degrees for each error node in the subgraph
        * num_activated_dets                 (int): number of activated detectors
        * num_candidate_errs                 (int): number of candidate error nodes that are connected to at least one activated detector node
        # * candidate_errs               (list[int]): list of error nodes that are connected to at least one activated detector node
        * candidate_errs_per_activated_det (float): number of candidate error nodes per activated detector node
        * num_boundary_dets                  (int): number of detector nodes that are connected to only one error node
    """
    dets = [n for n, d in G.nodes(data=True) if d.get('kind') == 'detector']
    errs = [n for n, d in G.nodes(data=True) if d.get('kind') == 'error']

    if len(dets) == 0:
        return {
            'WL_hash': '',
            'cycle_rank': 0,
            'inverse_girth': 0.0,
            'density': 0.0,
            'avg_degree': 0.0,
            'num_components': 0,
            'largest_component_size': 0,
            'avg_comp_dets': 0.0,
            'avg_comp_errs': 0.0,
            'num_activated_dets': 0,
            'num_candidate_errs': 0,
            'candidate_errs_per_activated_det': 0.0,
            'num_boundary_dets': 0,
        }

    activated_dets = [d for d, b in zip(dets, det_event) if b == '1']
    SG = local_incident_subgraph(G, activated_dets)

    cycle_rank = SG.number_of_edges() - SG.number_of_nodes() + nx.number_connected_components(SG)
    girth = nx.girth(SG)

    density = bipartite_graph_density(SG)
    avg_degree = np.mean([d for n, d in SG.degree()]) if SG.number_of_nodes() > 0 else 0.0

    comps = list(nx.connected_components(SG))
    largest_comp_size = max(len(c) for c in comps) if comps else 0
    n_comp = len(comps)

    num_comp_dets_dist = [sum(1 for n in c if G.nodes[n].get('kind') == 'detector') for c in comps]
    num_comp_errs_dist = [sum(1 for n in c if G.nodes[n].get('kind') == 'error') for c in comps]

    candidate_errs = [n for n, d in SG.nodes(data=True) if d.get('kind') == 'error']
    candidate_errs_per_activated_det = len(candidate_errs) / len(activated_dets) if len(activated_dets) > 0 else 0.0

    num_boundary_dets = sum(1 for n in candidate_errs if G.degree(n) == 1)

    return {
        'WL_hash': nx.weisfeiler_lehman_graph_hash(SG, node_attr='kind'),
        # graph features
        'cycle_rank': cycle_rank,
        'inverse_girth': 1 / girth if np.isfinite(girth) else 0.0,
        'density': density,
        'avg_degree': avg_degree,
        # component features
        'num_components': n_comp,
        'largest_component_size': largest_comp_size,
        'avg_comp_dets': np.mean(num_comp_dets_dist) if len(num_comp_dets_dist) > 0 else 0.0,
        'avg_comp_errs': np.mean(num_comp_errs_dist) if len(num_comp_errs_dist) > 0 else 0.0,
        # 'num_comp_dets_dist': num_comp_dets_dist,
        # 'num_comp_errs_dist': num_comp_errs_dist,
        # syndrome-induced features
        'num_activated_dets': len(activated_dets),
        'num_candidate_errs': len(candidate_errs),
        # 'candidate_errs': candidate_errs,
        'candidate_errs_per_activated_det': candidate_errs_per_activated_det,
        'num_boundary_dets': num_boundary_dets,
    }





def induced_subgraph(G: nx.Graph, active_dets: set[int]) -> nx.Graph:
    """
    Returns the subgraph of G induced by the active detectors.

    Notice that this is not just the subgraph of active detectors, but also includes any edges that are connected to at least one active detector.
    """
    affected_edges = [(u, v) for u, v in G.edges if u in active_dets or v in active_dets]
    SG = G.edge_subgraph(affected_edges).copy()
    for u in SG.nodes:
        SG.nodes[u]['active'] = True if u in active_dets else False
    return SG
