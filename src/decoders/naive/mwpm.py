from itertools import combinations

import numpy as np

from ..base import Decoder

try:
    import networkx as nx
except ImportError:  # pragma: no cover - project already depends on networkx elsewhere
    nx = None


def _ensure_binary_matrix(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(matrix, dtype=np.uint8) % 2


def _build_css_check_matrix(code) -> np.ndarray:
    hz = _ensure_binary_matrix(code.H_Z)
    hx = _ensure_binary_matrix(code.H_X)
    mz, nz = hz.shape
    mx, nx_ = hx.shape
    if nz != nx_:
        raise ValueError("H_Z and H_X must act on the same number of qubits")

    top = np.hstack([hz, np.zeros((mz, nz), dtype=np.uint8)])
    bottom = np.hstack([np.zeros((mx, nz), dtype=np.uint8), hx])
    return np.vstack([top, bottom]).astype(np.uint8)


def _validate_code(code) -> None:
    if code is None:
        raise ValueError("Decoder requires a code instance")
    for attr in ("H_X", "H_Z", "L_X", "L_Z", "n", "k", "c_x", "c_z"):
        if not hasattr(code, attr):
            raise ValueError(f"Code wrapper missing required attribute {attr!r}")


def _rref_binary(matrix: np.ndarray) -> tuple[np.ndarray, list[int]]:
    work = _ensure_binary_matrix(matrix).copy()
    rows, cols = work.shape
    pivot_cols: list[int] = []
    pivot_row = 0

    for col in range(cols):
        pivot = None
        for row in range(pivot_row, rows):
            if work[row, col]:
                pivot = row
                break
        if pivot is None:
            continue

        if pivot != pivot_row:
            work[[pivot_row, pivot]] = work[[pivot, pivot_row]]

        for row in range(rows):
            if row != pivot_row and work[row, col]:
                work[row] ^= work[pivot_row]

        pivot_cols.append(col)
        pivot_row += 1
        if pivot_row == rows:
            break

    return work, pivot_cols


def _solve_gf2(H: np.ndarray, syndrome: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    H = _ensure_binary_matrix(H)
    syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
    rows, cols = H.shape

    augmented = np.hstack([H.copy(), syndrome.reshape(rows, 1)])
    rref, pivot_cols = _rref_binary(augmented)

    for row in range(rows):
        if not rref[row, :cols].any() and rref[row, cols]:
            raise ValueError("Syndrome is inconsistent with the supplied check matrix")

    pivot_set = set(pivot_cols)
    free_cols = [col for col in range(cols) if col not in pivot_set]

    particular = np.zeros(cols, dtype=np.uint8)
    for row, pivot_col in enumerate(pivot_cols):
        particular[pivot_col] = rref[row, cols]

    nullspace_basis: list[np.ndarray] = []
    for free_col in free_cols:
        basis = np.zeros(cols, dtype=np.uint8)
        basis[free_col] = 1
        for row, pivot_col in enumerate(pivot_cols):
            basis[pivot_col] = rref[row, free_col]
        nullspace_basis.append(basis)

    return particular, nullspace_basis


def _min_weight_solution(H: np.ndarray, syndrome: np.ndarray, max_nullity: int = 20) -> np.ndarray:
    particular, nullspace_basis = _solve_gf2(H, syndrome)
    nullity = len(nullspace_basis)

    if nullity == 0:
        return particular
    if nullity > max_nullity:
        return particular

    best = particular.copy()
    best_weight = int(best.sum())

    for mask in range(1, 1 << nullity):
        candidate = particular.copy()
        for basis_idx, vector in enumerate(nullspace_basis):
            if (mask >> basis_idx) & 1:
                candidate ^= vector

        weight = int(candidate.sum())
        if weight < best_weight:
            best = candidate
            best_weight = weight

    return best


class MWPM(Decoder):
    """MWPM decoder for graph-like CSS sectors with an algebraic fallback elsewhere."""

    _BOUNDARY = "__boundary__"

    def __init__(self, code=None, error_rate: float = 0.001, **kwargs):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        if self.code is not None:
            self._init_backend()

    def _init_backend(self) -> None:
        _validate_code(self.code)
        self.H = _build_css_check_matrix(self.code)
        self.m, self.total_bits = self.H.shape
        self._graphlike_x = self._build_detector_graph(_ensure_binary_matrix(self.code.H_Z))
        self._graphlike_z = self._build_detector_graph(_ensure_binary_matrix(self.code.H_X))

    def _build_detector_graph(self, check_matrix: np.ndarray):
        if nx is None:
            return None

        checks, bits = check_matrix.shape
        graph = nx.Graph()
        graph.add_nodes_from(range(checks))
        graph.add_node(self._BOUNDARY)

        for bit_idx in range(bits):
            support = np.flatnonzero(check_matrix[:, bit_idx]).tolist()
            if len(support) == 0:
                continue
            if len(support) > 2:
                return None

            if len(support) == 1:
                u, v = support[0], self._BOUNDARY
            else:
                u, v = support

            if graph.has_edge(u, v):
                continue
            graph.add_edge(u, v, weight=1, qubit=bit_idx)

        return graph

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if not hasattr(self, "H"):
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        if syndrome.size != self.m:
            raise ValueError(f"syndrome has length {syndrome.size}, expected {self.m}")

        z_syndrome = syndrome[: self.code.c_z]
        x_syndrome = syndrome[self.code.c_z :]

        x_correction = self._decode_sector(_ensure_binary_matrix(self.code.H_Z), z_syndrome, self._graphlike_x)
        z_correction = self._decode_sector(_ensure_binary_matrix(self.code.H_X), x_syndrome, self._graphlike_z)
        return np.concatenate([x_correction, z_correction]).astype(np.uint8)

    def _decode_sector(
        self, check_matrix: np.ndarray, syndrome: np.ndarray, detector_graph
    ) -> np.ndarray:
        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        if syndrome.size != check_matrix.shape[0]:
            raise ValueError(
                f"Sector syndrome has length {syndrome.size}, expected {check_matrix.shape[0]}"
            )
        if syndrome.sum() == 0:
            return np.zeros(check_matrix.shape[1], dtype=np.uint8)

        if detector_graph is None or nx is None:
            return _min_weight_solution(check_matrix, syndrome)

        defects = np.flatnonzero(syndrome).tolist()
        if len(defects) == 0:
            return np.zeros(check_matrix.shape[1], dtype=np.uint8)

        complete = nx.Graph()
        defect_nodes = [f"d{i}" for i in range(len(defects))]
        copy_nodes = [f"c{i}" for i in range(len(defects))]
        pair_paths: dict[tuple[str, str], list[tuple[object, object]]] = {}
        boundary_paths: dict[str, list[tuple[object, object]]] = {}

        for i, j in combinations(range(len(defects)), 2):
            source = defects[i]
            target = defects[j]
            try:
                path = nx.shortest_path(detector_graph, source, target, weight="weight")
            except nx.NetworkXNoPath:
                continue

            edges = list(zip(path[:-1], path[1:]))
            weight = len(edges)
            u, v = defect_nodes[i], defect_nodes[j]
            complete.add_edge(u, v, weight=weight)
            pair_paths[(u, v)] = edges
            pair_paths[(v, u)] = list(reversed([(b, a) for a, b in edges]))

        for i, defect in enumerate(defects):
            try:
                path = nx.shortest_path(detector_graph, defect, self._BOUNDARY, weight="weight")
            except nx.NetworkXNoPath:
                continue

            edges = list(zip(path[:-1], path[1:]))
            defect_node = defect_nodes[i]
            copy_node = copy_nodes[i]
            complete.add_edge(defect_node, copy_node, weight=len(edges))
            boundary_paths[defect_node] = edges

        for i, j in combinations(range(len(copy_nodes)), 2):
            complete.add_edge(copy_nodes[i], copy_nodes[j], weight=0)

        matching = nx.algorithms.matching.min_weight_matching(
            complete, maxcardinality=True, weight="weight"
        )

        if len(matching) * 2 != len(complete.nodes):
            return _min_weight_solution(check_matrix, syndrome)

        correction = np.zeros(check_matrix.shape[1], dtype=np.uint8)
        for left, right in matching:
            if left.startswith("c") and right.startswith("c"):
                continue

            if left.startswith("d") and right.startswith("d"):
                edges = pair_paths.get((left, right), [])
            elif left.startswith("d") and right.startswith("c"):
                edges = boundary_paths.get(left, [])
            elif left.startswith("c") and right.startswith("d"):
                edges = boundary_paths.get(right, [])
            else:
                continue

            for u, v in edges:
                qubit = detector_graph[u][v]["qubit"]
                correction[qubit] ^= 1

        if np.array_equal((check_matrix @ correction) % 2, syndrome):
            return correction
        return _min_weight_solution(check_matrix, syndrome)
