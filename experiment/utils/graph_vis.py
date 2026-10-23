# failure_regions, force_osd_failures, low_weight_failures, quotient, trapping_sets
from typing import Literal

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from .misc import str2id

def extract_coord_from_graph(G: nx.Graph, dim: int = 3) -> dict[int, tuple[float, float]]:
    """
    Extracts the 'coord' attribute from graph nodes and returns a mapping from node_id to (x, y, t) coordinates.

    If a node does not have a 'coord' attribute, it defaults to (-2, -2, -2).
    """
    default_coord = (-2, -2, -2)
    coord = {}
    for node, data in G.nodes(data=True):
        coord[node] = tuple(data.get('coord', default_coord))
        if len(coord[node]) > dim:
            coord[node] = coord[node][:dim]
        if len(coord[node]) < dim:
            coord[node] = coord[node] + default_coord[len(coord[node]):dim]
    return coord

def extract_p_from_graph(G: nx.Graph) -> dict[int, float]:
    """
    Extracts the 'p' attribute from graph nodes and returns a mapping from node_id to (x, y, t) coordinates.

    If a node does not have a 'p' attribute, it defaults to -1.
    """
    default_p = -1
    p_values = {}
    for node, data in G.nodes(data=True):
        p_values[node] = data.get('p', default_p)
    return p_values

def flatten_graph(G: nx.Graph) -> nx.Graph:
    """
    Flatten graph node coordinates into 2D.

    Returns
    -------
    nx.Graph
        A copied graph with flattened 3D coordinates.
    """
    flattened = G.copy()
    pos = nx.kamada_kawai_layout(flattened)
    for node, (x, y) in pos.items():
        flattened.nodes[node]["coord"] = (x, y, 0.0)
        
    return flattened

# ==================================================
# Visualization of detection events on the graph
# ==================================================

from matplotlib import animation
from IPython.display import HTML

def visualize_det_on_graph_matplotlib(
    G: nx.Graph,
    Gtype: Literal['matching', 'tanner'],
    det_events: str,
):
    active = str2id(det_events)
    residual = set(G.nodes) - set(active) 
    if Gtype == 'matching':
        residual -= set(G.graph['virtual_boundary_nodes'])

    pos = extract_coord_from_graph(G, dim=3)
    
    active_nodes = np.array([pos[d] for d in active])
    if Gtype == 'matching':
        virtual_boundary_nodes = np.array([pos[d] for d in G.graph['virtual_boundary_nodes']])
    residual_nodes = np.array([pos[d] for d in residual])
    if Gtype == 'tanner':
        residual_check = np.array([pos[d] for d in residual if G.nodes[d].get('kind') == 'detector'])
        residual_variable = np.array([pos[d] for d in residual if G.nodes[d].get('kind') == 'error'])

    ae = [(u, v) for u, v in G.edges if u in active or v in active]
    affected_edges = np.array([(pos[u], pos[v]) for u, v in ae])
    residual_edges = np.array([(pos[u], pos[v]) for u, v in G.edges - ae])

    def _init():
        ax.clear()
        if Gtype == 'matching':
            ax.scatter(*active_nodes.T, alpha=0.2, s=100, color="red")
            ax.scatter(*virtual_boundary_nodes.T, alpha=0.2, s=100, color="green")
            ax.scatter(*residual_nodes.T, alpha=0.2, s=100, color="blue")
        else:
            ax.scatter(*active_nodes.T, alpha=0.2, s=100, color="red", marker='s')
            ax.scatter(*residual_check.T, alpha=0.2, s=100, color="blue", marker='s')
            ax.scatter(*residual_variable.T, alpha=0.2, s=100, color="blue", marker='o')
        for vizedge in affected_edges:
            ax.plot(*vizedge.T, color="red")
        for vizedge in residual_edges:
            ax.plot(*vizedge.T, color="gray")
        ax.grid(False)
        ax.set_axis_off()

    def _frame_update(index):
        ax.view_init(index * 0.2, index * 0.5)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    fig.tight_layout()

    ani = animation.FuncAnimation(
        fig,
        _frame_update,
        init_func=_init,
        interval=50,
        cache_frame_data=False,
        frames=100,
    )

    return ani

import plotly
import plotly.graph_objects as go
import numpy as np

def visualize_graph_plotly(
    G: nx.Graph,
):
    pos = extract_coord_from_graph(G, dim=3)
    p_values = extract_p_from_graph(G)

    dets = [n for n in G.nodes if G.nodes[n].get('kind') == 'detector']
    errs = [n for n in G.nodes if G.nodes[n].get('kind') == 'error']

    def edge_trace(edges, color):
        x, y, z = [], [], []
        for u, v in edges: # None to separate edges in the plot
            x.extend([pos[u][0], pos[v][0], None])
            y.extend([pos[u][1], pos[v][1], None])
            z.extend([pos[u][2], pos[v][2], None])
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=2, color=color),
            hoverinfo='none'
        )
    
    def node_trace(nodes, color, symbol, name):
        if len(nodes) == 0:
            x, y, z = [], [], []
            ps = []
        else:
            arr = np.array([pos[d] for d in nodes])            
            ps = [p_values[d] for d in nodes]
            x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=5, color=color, symbol=symbol),
            name=name,
            customdata=nodes,
            hovertext=ps,
            hovertemplate=(
                'x:%{x}<br>'
                'y:%{y}<br>'
                'z:%{z}<br>'
                '<extra>%{fullData.name}<br>id=%{customdata}<br>p=%{hovertext:.2f}</extra>'
            )
        )
    
    traces = []
    traces.append(node_trace(dets, 'blue', 'square', 'detector'))
    traces.append(node_trace(errs, 'red', 'circle', 'error'))
    traces.append(edge_trace(G.edges, 'gray'))

    fig = go.Figure(data=traces)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        showlegend=True
    )

    return fig

def visualize_det_on_tanner_graph(
    G: nx.Graph,
    det_events: str,
    other_events: dict[tuple[str, str, str], str] = None,
):
    """
    Visualization of detection events on the Tanner graph using plotly. The coloring scheme is as follows:

    - active check (red), residual check (blue square), residual variable (blue circle)
    - affected edges (red), residual edges (gray)
    
    Parameters
    ----------
    G: nx.Graph
        The graph to visualize. Nodes should have a 'coord' attribute for their 3D coordinates. For matching graphs, there should be a 'virtual_boundary_nodes' attribute in G.graph.
    det_events: str
        A string of '0' and '1' indicating which nodes are active (detection events). The length of this string should match the number of nodes in G.
    other_events: dict[tuple[str, str, str], str], optional
        A dictionary mapping event (names, colors, shapes) to strings of '0' and '1' for other types of events to visualize.
    """
    pos = extract_coord_from_graph(G, dim=3)
    p_values = extract_p_from_graph(G)

    def edge_trace(edges, color):
        x, y, z = [], [], []
        for u, v in edges: # None to separate edges in the plot
            x.extend([pos[u][0], pos[v][0], None])
            y.extend([pos[u][1], pos[v][1], None])
            z.extend([pos[u][2], pos[v][2], None])
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=2, color=color),
            hoverinfo='none'
        )
    
    def node_trace(nodes, color, symbol, name):
        if len(nodes) == 0:
            x, y, z = [], [], []
            ps = []
        else:
            arr = np.array([pos[d] for d in nodes])
            ps = [p_values[d] for d in nodes]
            x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=5, color=color, symbol=symbol),
            name=name,
            customdata=nodes,
            hovertext=ps,
            hovertemplate=(
                'x:%{x}<br>'
                'y:%{y}<br>'
                'z:%{z}<br>'
                '<extra>%{fullData.name}<br>id=%{customdata}<br>p=%{hovertext:.2f}</extra>'
            )
        )
    
    active = str2id(det_events)
    residual = set(G.nodes) - set(active)

    traces = []

    residual_check = [d for d in residual if G.nodes[d].get('kind') == 'detector']
    residual_variable = [d for d in residual if G.nodes[d].get('kind') == 'error']
    traces.append(node_trace(active, 'red', 'square', 'active'))
    traces.append(node_trace(residual_check, 'blue', 'square', 'check'))
    traces.append(node_trace(residual_variable, 'blue', 'circle', 'variable'))

    if other_events is not None:
        det_cnt = sum(1 for n in G.nodes if G.nodes[n].get('kind') == 'detector')
        for (name, color, shape), events in other_events.items():
            other = str2id(events, offset=det_cnt)
            traces.append(node_trace(other, color, shape, name))

    ae = [(u, v) for u, v in G.edges if u in active or v in active]
    re = [(u, v) for u, v in G.edges if (u, v) not in ae]

    traces.append(edge_trace(ae, 'red'))
    traces.append(edge_trace(re, 'gray'))

    fig = go.Figure(data=traces)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        showlegend=True
    )

    return fig

def visualize_det_on_matching_graph(
    G: nx.Graph,
    det_events: str,
    other_events: dict[tuple[str, str, str], str] = None,
):
    """
    Visualization of detection events on the matching graph using plotly. The coloring scheme is as follows:

    - active (red), virtual boundary (green), residual (blue)
    - affected edges (red), residual edges (gray)
    
    Parameters
    ----------
    G: nx.Graph
        The graph to visualize. Nodes should have a 'coord' attribute for their 3D coordinates. For matching graphs, there should be a 'virtual_boundary_nodes' attribute in G.graph.
    det_events: str
        A string of '0' and '1' indicating which nodes are active (detection events). The length of this string should match the number of nodes in G.
    other_events: dict[tuple[str, str, str], str], optional
        A dictionary mapping event (names, colors, shapes) to strings of '0' and '1' for other types of events to visualize.
    """
    pos = extract_coord_from_graph(G, dim=3)
    p_values = extract_p_from_graph(G)

    def edge_trace(edges, color):
        x, y, z = [], [], []
        for u, v in edges: # None to separate edges in the plot
            x.extend([pos[u][0], pos[v][0], None])
            y.extend([pos[u][1], pos[v][1], None])
            z.extend([pos[u][2], pos[v][2], None])
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=2, color=color),
            hoverinfo='none'
        )
    
    def node_trace(nodes, color, symbol, name):
        if len(nodes) == 0:
            x, y, z = [], [], []
            ps = []
        else:
            arr = np.array([pos[d] for d in nodes])
            ps = [p_values[d] for d in nodes]
            x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=5, color=color, symbol=symbol),
            name=name,
            customdata=nodes,
            hovertext=ps,
            hovertemplate=(
                'x:%{x}<br>'
                'y:%{y}<br>'
                'z:%{z}<br>'
                '<extra>%{fullData.name}<br>id=%{customdata}<br>p=%{hovertext:.2f}</extra>'
            )
        )
    
    active = str2id(det_events)
    residual = set(G.nodes) - set(active) - set(G.graph['virtual_boundary_nodes'])

    traces = []

    traces.append(node_trace(active, 'red', 'circle', 'active'))
    traces.append(node_trace(residual, 'blue', 'circle', 'inactive'))
    traces.append(node_trace(G.graph['virtual_boundary_nodes'], 'green', 'circle', 'boundary'))
    
    if other_events is not None:
        det_cnt = sum(1 for n in G.nodes if G.nodes[n].get('kind') == 'detector')
        for (name, color, shape), events in other_events.items():
            other = str2id(events, offset=det_cnt)
            traces.append(node_trace(other, color, shape, name))

    ae = [(u, v) for u, v in G.edges if u in active or v in active]
    re = [(u, v) for u, v in G.edges if (u, v) not in ae]

    traces.append(edge_trace(ae, 'red'))
    traces.append(edge_trace(re, 'gray'))

    fig = go.Figure(data=traces)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        showlegend=True
    )

    return fig

def animate_hard_decisions_on_graph(
    G: nx.Graph,
    det_events: str,
    hard_decisions_list: list[str],
    err_events: str = None,
):
    pos = extract_coord_from_graph(G, dim=3)
    p_values = extract_p_from_graph(G)

    detector_nodes = sorted([n for n in G.nodes if G.nodes[n].get('kind') == 'detector'])
    error_nodes = sorted([n for n in G.nodes if G.nodes[n].get('kind') == 'error'])
    active_dets = str2id(det_events)
    active_errs = str2id(err_events, offset=len(detector_nodes)) if err_events else []

    def edge_trace(edges, color, name):
        x, y, z = [], [], []
        for u, v in edges:
            x.extend([pos[u][0], pos[v][0], None])
            y.extend([pos[u][1], pos[v][1], None])
            z.extend([pos[u][2], pos[v][2], None])
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=2, color=color),
            hoverinfo='none',
            name=name,
        )

    def node_trace(nodes, color, symbol, name):
        if len(nodes) == 0:
            x, y, z = [], [], []
            ps = []
        else:
            arr = np.array([pos[d] for d in nodes])
            ps = [p_values[d] for d in nodes]
            x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=5, color=color, symbol=symbol),
            name=name,
            customdata=nodes,
            hovertext=ps,
            hovertemplate=(
                'x:%{x}<br>'
                'y:%{y}<br>'
                'z:%{z}<br>'
                '<extra>%{fullData.name}<br>id=%{customdata}<br>p=%{hovertext:.2f}</extra>'
            )
        )

    def frame_traces(hard_decisions):
        decided_errs = [error_node for error_node, decision in zip(error_nodes, hard_decisions) if decision == '1']
        undecided_errs = [node for node in error_nodes if node not in decided_errs and node not in active_errs]
        untouched_errs = [node for node in active_errs if node not in decided_errs]

        return [
            node_trace(detector_nodes, 'blue', 'square', 'check'),
            node_trace(active_dets, 'skyblue', 'square', 'activated check'),
            node_trace(undecided_errs, 'green', 'circle', '0'),
            node_trace(decided_errs, 'red', 'circle', '1'),
            node_trace(untouched_errs, 'red', 'circle-open', 'untouched'),
            edge_trace(G.edges, 'gray', 'edge'),
        ]

    frames = [
        go.Frame(data=frame_traces(hard_decisions), name=str(time_step))
        for time_step, hard_decisions in enumerate(hard_decisions_list)
    ]

    fig = go.Figure(data=frames[0].data, frames=frames)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        showlegend=True,
        updatemenus=[
            dict(
                type='buttons',
                showactive=False,
                buttons=[
                    dict(
                        label='Play',
                        method='animate',
                        args=[
                            None,
                            dict(
                                frame=dict(duration=500, redraw=True),
                                fromcurrent=True,
                                transition=dict(duration=0),
                            ),
                        ],
                    ),
                    dict(
                        label='Pause',
                        method='animate',
                        args=[
                            [None],
                            dict(
                                frame=dict(duration=0, redraw=True),
                                mode='immediate',
                                transition=dict(duration=0),
                            ),
                        ],
                    ),
                ],
            )
        ],
        sliders=[
            dict(
                active=0,
                currentvalue=dict(prefix='time step: '),
                steps=[
                    dict(
                        label=str(time_step),
                        method='animate',
                        args=[
                            [str(time_step)],
                            dict(
                                frame=dict(duration=0, redraw=True),
                                mode='immediate',
                                transition=dict(duration=0),
                            ),
                        ],
                    )
                    for time_step in range(len(hard_decisions_list))
                ],
            )
        ],
    )

    return fig


def visualize_det_on_graph_plotly(
    G: nx.Graph, 
    Gtype: Literal['matching', 'tanner'],
    det_events: str,
    other_events: dict | None = None,
):
    assert Gtype in {"matching", "tanner"}, Gtype
    pos = extract_coord_from_graph(G, dim=3)

    def edge_trace(edges, color):
        x, y, z = [], [], []
        for u, v in edges: # None to separate edges in the plot
            x.extend([pos[u][0], pos[v][0], None])
            y.extend([pos[u][1], pos[v][1], None])
            z.extend([pos[u][2], pos[v][2], None])
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=2, color=color),
            hoverinfo='none'
        )

    def node_trace(nodes, color, symbol, name):
        arr = np.array([pos[d] for d in nodes]).reshape(-1, 3)
        x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        return go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=5, color=color, symbol=symbol),
            name=name
        )

    active = str2id(det_events)
    assert set(active) <= set(G.nodes), active
    residual = set(G.nodes) - set(active)

    if Gtype == 'matching':
        residual -= set(G.graph['virtual_boundary_nodes'])

    traces = []

    if Gtype == 'matching':
        traces.append(node_trace(active, 'red', 'circle', 'active'))
        traces.append(node_trace(residual, 'blue', 'circle', 'inactive'))
        traces.append(node_trace(G.graph['virtual_boundary_nodes'], 'green', 'circle', 'boundary'))
    else:
        residual_check = [d for d in residual if G.nodes[d].get('kind') == 'detector']
        residual_variable = [d for d in residual if G.nodes[d].get('kind') == 'error']
        traces.append(node_trace(active, 'red', 'square', 'active'))
        traces.append(node_trace(residual_check, 'blue', 'square', 'check'))
        traces.append(node_trace(residual_variable, 'blue', 'circle', 'variable'))

    if other_events is not None:
        for (label, color, symbol), bits in other_events.items():
            assert label == "fault" and Gtype == "tanner", label
            detectors = [n for n in G.nodes if G.nodes[n].get("kind") == "detector"]
            nodes = [len(detectors) + i for i in str2id(bits)]
            assert all(n in G and G.nodes[n].get("kind") == "error" for n in nodes), nodes
            traces.append(node_trace(nodes, color, symbol, label))

    ae = [(u, v) for u, v in G.edges if u in active or v in active]
    re = [(u, v) for u, v in G.edges if (u, v) not in ae]

    traces.append(edge_trace(ae, 'red'))
    traces.append(edge_trace(re, 'gray'))

    fig = go.Figure(data=traces)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        showlegend=True
    )

    return fig
