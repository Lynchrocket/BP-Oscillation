"""
1. Code-capacity
    DEPOLARIZE1(p) on data qubits

2. Phenomenological
    extends code-capacity with d_rounds rounds of 
      * DEPOLARIZE1(p) on data qubits
      * BIT_FLIP(p)    on syndrome measurement outcomes (classical readout error)
    
3. Circuit-level
    circuit noise model:
      * DEPOLARIZE1(p) on idle/data qubits
      * DEPOLARIZE2(p) after every two-qubit CX gate
      * DEPOLARIZE1(p) after ancilla resets (preparation noise)
      * BIT_FLIP(p)    on syndrome measurements outcomes (classical readout noise)
"""

import numpy as np
import scipy.sparse as sp
import networkx as nx
import stim
from collections import defaultdict
from typing import Literal
from dataclasses import dataclass
import random

# ---------------------------------------------------------------------------
# Error Configuration
# ---------------------------------------------------------------------------

@dataclass
class ErrorRate:
    p: float                   # data qubit (and the default value for other error rates)
    p_meas: float = None       # measurement readout
    p_cx: float = None         # two-qubit gate
    p_prep: float = None       # preparation noise (reset)
    p_idle: float = None       # idle qubit

    def __post_init__(self):   # set default values for unspecified error rates
        if self.p_meas is None:
            self.p_meas = self.p
        if self.p_cx is None:
            self.p_cx = self.p
        if self.p_prep is None:
            self.p_prep = self.p
        if self.p_idle is None:
            self.p_idle = self.p

# ---------------------------------------------------------------------------
# CNOT Schedule
# ---------------------------------------------------------------------------

def _edge_coloring(H: np.ndarray) -> tuple[defaultdict[int, list[tuple[int, int]]], defaultdict[int, list[int]]]:
    """
    Parameters
    ----------
    H : np.ndarray
        Binary parity-check matrix.

    Returns
    -------
    color_to_edges : defaultdict[int, list[tuple[int, int]]]
        Colour index -> list of (check_idx, data_idx).
    idle_qubits : defaultdict[int, list[int]]
        Colour index -> list of idle qubit indices.
    """
    m, n = H.shape
    all_qubits = set(range(n + m))
    # Build bipartite graph: check nodes are offset by cols to avoid index collisions with data nodes.
    G = nx.Graph()
    edges = np.where(H == 1)
    for i, j in zip(*edges):
        G.add_edge(i + n, j)

    coloring = nx.coloring.greedy_color(nx.line_graph(G), strategy="largest_first")

    # color -> list[(check_idx, data_idx)]
    color_to_edges: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
    idle_qubits: defaultdict[int, list[int]] = defaultdict(list)
    for (u, v), color in coloring.items():
        row, col = (v - n, u) if u < n else (u - n, v)
        color_to_edges[color].append((row, col))
    for color, edges in color_to_edges.items():
        used_qubits = set()
        for e in edges:
            used_qubits.add(e[0] + n) # check qubit
            used_qubits.add(e[1]) # data qubit
        idle_qubits[color].extend(all_qubits - used_qubits)
    return color_to_edges, idle_qubits

def _random_schedule(H: np.ndarray, seed=None):
    rng = random.Random(seed)
    m, n = H.shape
    all_qubits = set(range(n + m))
    edges = list(zip(*np.where(H == 1)))

    color_to_edges: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
    idle_qubits: defaultdict[int, list[int]] = defaultdict(list)

    color = 0
    while len(edges) > 0:
        rng.shuffle(edges)
        used_qubits = set()
        layer = []
        next_edges = []
        for ec, ev in edges:
            if (ec + n) not in used_qubits and ev not in used_qubits:
                layer.append((ec, ev))
                used_qubits.add(ec + n) # check qubit
                used_qubits.add(ev)     # data qubit
            else:
                next_edges.append((ec, ev))
        color_to_edges[color].extend(layer)
        idle_qubits[color].extend(all_qubits - used_qubits)
        edges = next_edges
        color += 1
    return color_to_edges, idle_qubits

def _reversed_schedule(edges, idles):
    colors = sorted(edges.keys())

    reversed_edges: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
    reversed_idles: defaultdict[int, list[int]] = defaultdict(list)

    for new_color, old_color in enumerate(reversed(colors)):
        reversed_edges[new_color].extend(edges[old_color])
        reversed_idles[new_color].extend(idles[old_color])

    return reversed_edges, reversed_idles

_schedule_strategies = {
    0: lambda H, seed=None: _edge_coloring(H),
    1: lambda H, seed=None: _random_schedule(H, seed=seed),
    2: lambda H, seed=None: _reversed_schedule(*_edge_coloring(H)),
}

# ---------------------------------------------------------------------------
# CircuitBuilder
# ---------------------------------------------------------------------------

class CircuitBuilder:
    """Syndrome extraction circuit builder for CA codes."""
    def __init__(self, code, need_tags = False, schedule: int = 0, seed: int = 3407):
        self.code = code
        self.need_tags = need_tags
        self.schedule = schedule
        self.seed = seed

        self.H = np.asarray(
            code.H.toarray() if sp.issparse(code.H) else code.H, dtype=np.uint8
        )
        self.L_X = np.asarray(
            code.L_X.toarray() if sp.issparse(code.L_X) else code.L_X, dtype=np.uint8
        )
        self.L_Z = np.asarray(
            code.L_Z.toarray() if sp.issparse(code.L_Z) else code.L_Z, dtype=np.uint8
        )

        self.m = self.H.shape[0]
        self.n = self.H.shape[1]
        self.k = self.L_Z.shape[0]

        self.datas = list(range(self.n))
        self.ancillas  = list(range(self.n, self.n + self.m))

        self.coloring, self.idle_qubits = _schedule_strategies[self.schedule](self.H, self.seed)
        self.stab_coords = code.stabilizer_coords() if hasattr(code, "stabilizer_coords") else None
        self.qubit_coords = code.qubit_coords() if hasattr(code, "qubit_coords") else None

    def get_circuit(
        self, 
        noise_model: Literal["code_capacity", "phenomenological", "circuit_level"],
        error_rate: ErrorRate,
        d_rounds: int = None
    ) -> stim.Circuit:
        """
        Parameters
        ----------
        noise_model : str
            One of "code_capacity", "phenomenological", "circuit_level".
        error_rate : ErrorRate
            Error rates for different types of noise.
        d_rounds : int, optional
            Number of syndrome measurement rounds (only for phenomenological and circuit-level).
        """
        match noise_model:
            case "code_capacity":
                return self.get_code_capacity_circuit(error_rate)
            case "phenomenological":
                return self.get_phenomenological_circuit(error_rate, d_rounds)
            case "circuit_level":
                return self.get_circuit_level_circuit(error_rate, d_rounds)
            case _:
                raise ValueError(f"Unknown noise model: {noise_model}")

    def get_code_capacity_circuit(self, error_rate: ErrorRate) -> stim.Circuit:
        """
        Noise model: DEPOLARIZE1(p) on data qubits, perfect measurements.
        
        Parameters
        ----------
        error_rate : ErrorRate
            Error rates for different types of noise. Only `p` is used for code-capacity.
        """
        return self._build_code_capacity(error_rate.p)

    def get_phenomenological_circuit(self, error_rate: ErrorRate, d_rounds: int = None) -> stim.Circuit:
        """Multi-round phenomenological circuit.

        Parameters
        ----------
        error_rate : ErrorRate
            Error rates for different types of noise. Only `p` and `p_meas` are used for phenomenological.
        d_rounds : int, optional
        """
        if d_rounds is None:
            d_rounds = max(
                int(min(self.L_Z.sum(axis=1).min(), self.L_X.sum(axis=1).min())), 1
            )
        return self._build_phenomenological(error_rate.p, error_rate.p_meas, d_rounds)

    def get_circuit_level_circuit(self, error_rate: ErrorRate, d_rounds: int = None) -> stim.Circuit:
        """
        Parameters
        ----------
        error_rate : ErrorRate
            Error rates for different types of noise.
        d_rounds : int, optional
        """
        if d_rounds is None:
            d_rounds = max(
                int(min(self.L_Z.sum(axis=1).min(), self.L_X.sum(axis=1).min())), 1
            )
        return self._build_circuit_level(error_rate, d_rounds)

    # ------------------------------------------------------------------
    # MPP Z-stabiliser measurements
    # ------------------------------------------------------------------
    
    def _mpp_z_stabilisers(self, circuit: stim.Circuit, flip: float = 0.0) -> None:
        """Append one MPP Z-stabiliser measurement per row of H_Z.

        Parameters
        ----------
        circuit : stim.Circuit
        flip : float
            If > 0, each measurement outcome is independently flipped with
            this probability (readout noise).
        """
        for row in self.H:
            targets = []
            for j in np.where(row)[0]:
                targets.append(stim.target_z(int(j)))
                targets.append(stim.target_combiner())
            targets.pop()  # remove trailing combiner
            if flip > 0:
                # BIT_FLIP on the measurement outcome
                circuit.append("MPP", targets, flip, tag='measurement' if self.need_tags else '')
            else:
                circuit.append("MPP", targets)

    # ------------------------------------------------------------------
    # Code-capacity
    # ------------------------------------------------------------------

    def _build_code_capacity(self, p: float) -> stim.Circuit:
        """
        Parameters
        ----------
        p : float
            Physical error rate for DEPOLARIZE1.
        """
        n, m = self.n, self.m
        circ  = stim.Circuit()

        # Noiseless reference round: project onto +1 eigenspace of all Z-stabilisers.
        self._mpp_z_stabilisers(circ, flip=0.0)

        # DEPOLARIZE1(p) on data qubits
        circ.append("DEPOLARIZE1", self.datas, p, tag='data' if self.need_tags else '')

        # Noiseless measurement round
        self._mpp_z_stabilisers(circ, flip=0.0)

        for i in range(m):
            coords = [*self.stab_coords[i], 0] if self.stab_coords is not None else [i, 0]
            circ.append(
                "DETECTOR",
                [stim.target_rec(-m + i), stim.target_rec(-2 * m + i)],
                coords,
            )
        
        circ.append("M", self.datas)

        for i, row in enumerate(self.L_Z):
            targets = [stim.target_rec(-n + j) for j, v in enumerate(row) if v == 1]
            circ.append("OBSERVABLE_INCLUDE", targets, i)

        return circ

    # ------------------------------------------------------------------
    # Phenomenological
    # ------------------------------------------------------------------
    
    def _build_phenomenological(self, p: float, p_meas: float, d_rounds: int) -> stim.Circuit:
        """
        Parameters
        ----------
        p : float
            Physical error rate for DEPOLARIZE1.
        p_meas : float
            Physical error rate for Measurement.
        d_rounds : int
            Number of noisy data-error rounds.
        """
        n, m = self.n, self.m
        circ  = stim.Circuit()

        # Noiseless reference round: project onto +1 eigenspace of all Z-stabilisers.
        self._mpp_z_stabilisers(circ, flip=0.0)

        for r in range(d_rounds):
            # DEPOLARIZE1(p) on data qubits
            circ.append("DEPOLARIZE1", self.datas, p, tag='data' if self.need_tags else '')

            self._mpp_z_stabilisers(circ, flip=p_meas)

            for i in range(m):
                coords = [*self.stab_coords[i], 3*r] if self.stab_coords is not None else [i, r]
                circ.append(
                    "DETECTOR",
                    [stim.target_rec(-m + i), stim.target_rec(-2 * m + i)],
                    coords,
                )
        
        if p_meas > 0:
            circ.append("M", self.datas, p_meas, tag='measurement' if self.need_tags else '')
        else:
            circ.append("M", self.datas)

        for i, row in enumerate(self.H):
            coords = [*self.stab_coords[i], 3*d_rounds] if self.stab_coords is not None else [i, d_rounds]
            targets = [stim.target_rec(-n + j) for j, v in enumerate(row) if v == 1]
            targets.append(stim.target_rec(-n - m + i))
            circ.append("DETECTOR", targets, coords)

        for i, row in enumerate(self.L_Z):
            targets = [stim.target_rec(-n + j) for j, v in enumerate(row) if v == 1]
            circ.append("OBSERVABLE_INCLUDE", targets, i)

        return circ

    # ------------------------------------------------------------------
    # Circuit-level
    # ------------------------------------------------------------------

    def _one_syndrome_round(self, circ: stim.Circuit, error_rate: ErrorRate, round: int) -> None:
        """
        Parameters
        ----------
        circ : stim.Circuit
        error_rate : ErrorRate
        round: int
        """
        for color, edges in sorted(self.coloring.items()):
            targets = []
            for check_idx, data_idx in edges:
                targets.extend([self.datas[data_idx], self.ancillas[check_idx]])
            circ.append("CX", targets)
            # DEPOLARIZE2(p) after every CX gate
            circ.append("DEPOLARIZE2", targets, error_rate.p_cx, tag=f'cx({round},0)' if self.need_tags else '')
            circ.append("TICK")
            # DPOLARIZE1(p) on idle qubits
            circ.append("DEPOLARIZE1", self.idle_qubits[color], error_rate.p_idle, tag=f'idle({round},1)' if self.need_tags else '')
            circ.append("TICK")

        # BIT_FLIP(p) on measurement outcomes
        if error_rate.p_meas > 0:
            circ.append("M", self.ancillas, error_rate.p_meas, tag=f'measurement({round},2)' if self.need_tags else '')
        else:
            circ.append("M", self.ancillas, tag=f'measurement({round},2)' if self.need_tags else '')
        # DEPOLARIZE1(p) on idle qubits during measurement and reset
        circ.append("DEPOLARIZE1", self.datas, error_rate.p_idle, tag=f'idle({round},2)' if self.need_tags else '')
        circ.append("TICK")
        circ.append("R", self.ancillas)
        circ.append("DEPOLARIZE1", self.datas, error_rate.p_idle, tag=f'idle({round},3)' if self.need_tags else '')
        circ.append("TICK")
        # DEPOLARIZE1(p) on ancilla reset
        circ.append("DEPOLARIZE1", self.ancillas, error_rate.p_prep, tag=f'reset({round},4)' if self.need_tags else '')
        circ.append("TICK")

    def _build_circuit_level(self, error_rate: ErrorRate, d_rounds: int) -> stim.Circuit:
        """
        Parameters
        ----------
        error_rate : ErrorRate
        d_rounds : int
        """
        n, m = self.n, self.m
        circ = stim.Circuit()

        circ.append("R", self.datas + self.ancillas)
        # DEPOLARIZE1(p) on preparation noise
        circ.append("DEPOLARIZE1", self.datas + self.ancillas, error_rate.p_prep, tag='reset(-1,0)' if self.need_tags else '')
        circ.append("TICK")

        for r in range(d_rounds):
            self._one_syndrome_round(circ, error_rate, r)

            if r == 0:
                for k in range(m):
                    coords = [*self.stab_coords[k], 3*r] if self.stab_coords is not None else [k, 0]
                    circ.append("DETECTOR", [stim.target_rec(-m + k)], coords)
            else:
                for k in range(m):
                    coords = [*self.stab_coords[k], 3*r] if self.stab_coords is not None else [k, r]
                    circ.append(
                        "DETECTOR",
                        [stim.target_rec(-m + k), stim.target_rec(-2 * m + k)],
                        coords,
                    ) 

        if error_rate.p_meas > 0:
            circ.append("M", self.datas, error_rate.p_meas, tag=f'measurement({d_rounds},0)' if self.need_tags else '')
        else:
            circ.append("M", self.datas, tag=f'measurement({d_rounds},0)' if self.need_tags else '')

        for i, row in enumerate(self.H):
            coords = [*self.stab_coords[i], 3*d_rounds] if self.stab_coords is not None else [i, d_rounds]
            targets = [stim.target_rec(-n + j) for j, v in enumerate(row) if v == 1]
            targets.append(stim.target_rec(-n - m + i))
            circ.append("DETECTOR", targets, coords)

        for i, row in enumerate(self.L_Z):
            targets = [stim.target_rec(-n + j) for j, v in enumerate(row) if v == 1]
            circ.append("OBSERVABLE_INCLUDE", targets, i)

        return circ