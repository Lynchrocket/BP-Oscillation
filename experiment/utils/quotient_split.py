# quotient
"""
quotient_split.py

Quotient Successive-Cancellation List BP
========================================

Implements the decomposition

    x = x_s XOR D u XOR Q ell,

and a beam/list decoder that starts from [D, Q], successively conditions
logical coordinates, keeps the best affine logical blocks, and finally ranks
surviving fixed-coset leaves using a D-only Bethe partition estimate.

Only NumPy is required.
"""

from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

Array = np.ndarray
EdgeKey = Tuple[int, int]


def as_binary_array(x, *, ndim: Optional[int] = None) -> Array:
    arr = np.asarray(x, dtype=np.uint8) & 1
    if ndim is not None and arr.ndim != ndim:
        raise ValueError(f"Expected ndim={ndim}, got shape {arr.shape}.")
    return arr


def gf2_rank(a: Array) -> int:
    m = as_binary_array(a, ndim=2).copy()
    rows, cols = m.shape
    rank = 0
    for col in range(cols):
        candidates = np.flatnonzero(m[rank:, col])
        if candidates.size == 0:
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            m[[rank, pivot]] = m[[pivot, rank]]
        for row in range(rows):
            if row != rank and m[row, col]:
                m[row] ^= m[rank]
        rank += 1
        if rank == rows:
            break
    return rank


def gf2_solve(a: Array, b: Array) -> Array:
    """Return one solution of A x = b over GF(2), with free variables set to 0."""
    a = as_binary_array(a, ndim=2)
    b = as_binary_array(b, ndim=1)
    if a.shape[0] != b.size:
        raise ValueError("A and b have incompatible dimensions.")

    aug = np.concatenate([a.copy(), b[:, None]], axis=1)
    rows, cols = a.shape
    pivot_cols: List[int] = []
    pivot_row = 0

    for col in range(cols):
        candidates = np.flatnonzero(aug[pivot_row:, col])
        if candidates.size == 0:
            continue
        row = pivot_row + int(candidates[0])
        if row != pivot_row:
            aug[[pivot_row, row]] = aug[[row, pivot_row]]
        for other in range(rows):
            if other != pivot_row and aug[other, col]:
                aug[other] ^= aug[pivot_row]
        pivot_cols.append(col)
        pivot_row += 1
        if pivot_row == rows:
            break

    for row in range(pivot_row, rows):
        if not aug[row, :cols].any() and aug[row, cols]:
            raise ValueError("The GF(2) system is inconsistent.")

    x = np.zeros(cols, dtype=np.uint8)
    for row, col in enumerate(pivot_cols):
        x[col] = aug[row, cols]
    return x


def logical_neutral_representative(h: Array, l: Array, q: Array, syndrome: Array) -> Array:
    """Construct x_s with H x_s = syndrome and L x_s = 0."""
    h = as_binary_array(h, ndim=2)
    l = as_binary_array(l, ndim=2)
    q = as_binary_array(q, ndim=2)
    syndrome = as_binary_array(syndrome, ndim=1)

    if h.shape[1] != l.shape[1] or q.shape != (h.shape[1], l.shape[0]):
        raise ValueError("H, L, and Q have incompatible shapes.")
    if np.any((h @ q) & 1):
        raise ValueError("Q does not satisfy H Q = 0.")
    if not np.array_equal((l @ q) & 1, np.eye(l.shape[0], dtype=np.uint8)):
        raise ValueError("Q does not satisfy L Q = I.")

    x0 = gf2_solve(h, syndrome)
    x_s = x0 ^ ((q @ ((l @ x0) & 1)) & 1)
    if not np.array_equal((h @ x_s) & 1, syndrome) or np.any((l @ x_s) & 1):
        raise RuntimeError("Failed to construct a logical-neutral representative.")
    return x_s


def validate_quotient_decomposition(
    d_matrix: Array,
    q_matrix: Array,
    *,
    h: Optional[Array] = None,
    l: Optional[Array] = None,
) -> None:
    d_matrix = as_binary_array(d_matrix, ndim=2)
    q_matrix = as_binary_array(q_matrix, ndim=2)
    if d_matrix.shape[0] != q_matrix.shape[0]:
        raise ValueError("D and Q must have the same number of rows.")
    g = np.concatenate([d_matrix, q_matrix], axis=1)
    if gf2_rank(g) != g.shape[1]:
        raise ValueError("[D, Q] is not full-column-rank.")

    if h is None and l is None:
        return
    if h is None or l is None:
        raise ValueError("Supply both H and L, or neither.")

    h = as_binary_array(h, ndim=2)
    l = as_binary_array(l, ndim=2)
    n, d_count = d_matrix.shape
    k = q_matrix.shape[1]
    if h.shape[1] != n or l.shape != (k, n):
        raise ValueError("H or L has an incompatible shape.")
    if np.any((h @ d_matrix) & 1):
        raise ValueError("D does not satisfy H D = 0.")
    if np.any((l @ d_matrix) & 1):
        raise ValueError("D does not satisfy L D = 0.")
    if np.any((h @ q_matrix) & 1):
        raise ValueError("Q does not satisfy H Q = 0.")
    if not np.array_equal((l @ q_matrix) & 1, np.eye(k, dtype=np.uint8)):
        raise ValueError("Q does not satisfy L Q = I.")


@dataclass
class BPResult:
    converged: bool
    iterations: int
    variable_prob_1: Array
    log_z_bethe: float
    max_message_change: float
    message_state: Dict[EdgeKey, Tuple[Array, Array]]


class GeneratorGraphBP:
    """Sum-product BP for Z = sum_w product_i P_i(offset_i XOR (G w)_i)."""

    def __init__(
        self,
        generator: Array,
        offset: Array,
        p_error: Array,
        *,
        coordinate_ids: Optional[Sequence[int]] = None,
        epsilon: float = 1e-15,
    ) -> None:
        self.g = as_binary_array(generator, ndim=2)
        self.offset = as_binary_array(offset, ndim=1)
        self.p_error = np.asarray(p_error, dtype=float)
        self.epsilon = float(epsilon)

        n, r = self.g.shape
        if self.offset.shape != (n,) or self.p_error.shape != (n,):
            raise ValueError("generator, offset, and p_error have incompatible dimensions.")
        if np.any(~np.isfinite(self.p_error)) or np.any((self.p_error <= 0) | (self.p_error >= 1)):
            raise ValueError("Every p_error[i] must lie strictly between 0 and 1.")

        if coordinate_ids is None:
            self.coordinate_ids = tuple(range(r))
        else:
            if len(coordinate_ids) != r or len(set(coordinate_ids)) != r:
                raise ValueError("coordinate_ids must be unique and match generator columns.")
            self.coordinate_ids = tuple(int(x) for x in coordinate_ids)

        self.row_edges: List[List[int]] = [[] for _ in range(n)]
        self.col_edges: List[List[int]] = [[] for _ in range(r)]
        edge_rows: List[int] = []
        edge_cols: List[int] = []
        for i, j in np.argwhere(self.g == 1):
            edge = len(edge_rows)
            edge_rows.append(int(i))
            edge_cols.append(int(j))
            self.row_edges[int(i)].append(edge)
            self.col_edges[int(j)].append(edge)

        self.edge_rows = np.asarray(edge_rows, dtype=np.int64)
        self.edge_cols = np.asarray(edge_cols, dtype=np.int64)
        self.edge_keys = [
            (int(i), self.coordinate_ids[int(j)])
            for i, j in zip(self.edge_rows, self.edge_cols)
        ]

    def _normalize(self, values: Array) -> Array:
        out = np.maximum(np.asarray(values, dtype=float), self.epsilon)
        total = float(out.sum())
        if not np.isfinite(total) or total <= 0:
            return np.array([0.5, 0.5], dtype=float)
        return out / total

    def _initialize(self, initial_message_state):
        e_count = len(self.edge_keys)
        vf = np.full((e_count, 2), 0.5, dtype=float)
        fv = np.full((e_count, 2), 0.5, dtype=float)
        if initial_message_state is not None:
            for edge, key in enumerate(self.edge_keys):
                stored = initial_message_state.get(key)
                if stored is not None:
                    vf[edge] = self._normalize(stored[0])
                    fv[edge] = self._normalize(stored[1])
        return vf, fv

    def _factor_update(self, vf: Array) -> Array:
        out = np.empty_like(vf)
        for i, edges in enumerate(self.row_edges):
            degree = len(edges)
            if degree == 0:
                continue
            biases = vf[edges, 0] - vf[edges, 1]
            prefix = np.ones(degree + 1)
            suffix = np.ones(degree + 1)
            for t in range(degree):
                prefix[t + 1] = prefix[t] * biases[t]
            for t in range(degree - 1, -1, -1):
                suffix[t] = suffix[t + 1] * biases[t]

            likelihood = np.array([1.0 - self.p_error[i], self.p_error[i]])
            off = int(self.offset[i])
            for local, edge in enumerate(edges):
                parity_bias = prefix[local] * suffix[local + 1]
                p_even = 0.5 * (1.0 + parity_bias)
                p_odd = 0.5 * (1.0 - parity_bias)
                msg = np.empty(2)
                for bit in (0, 1):
                    x_even = off ^ bit
                    msg[bit] = p_even * likelihood[x_even] + p_odd * likelihood[x_even ^ 1]
                out[edge] = self._normalize(msg)
        return out

    def _variable_update(self, fv: Array) -> Array:
        out = np.empty_like(fv)
        for edges in self.col_edges:
            if not edges:
                continue
            logs = np.log(np.maximum(fv[edges], self.epsilon))
            total = logs.sum(axis=0)
            for local, edge in enumerate(edges):
                raw = total - logs[local]
                raw -= np.max(raw)
                out[edge] = self._normalize(np.exp(raw))
        return out

    def _variable_beliefs(self, fv: Array) -> Array:
        result = np.empty(self.g.shape[1], dtype=float)
        for j, edges in enumerate(self.col_edges):
            if not edges:
                result[j] = 0.5
                continue
            log_b = np.log(np.maximum(fv[edges], self.epsilon)).sum(axis=0)
            log_b -= np.max(log_b)
            result[j] = self._normalize(np.exp(log_b))[1]
        return result

    def _bethe_log_partition(self, vf: Array, fv: Array) -> float:
        log_z = 0.0

        for i, edges in enumerate(self.row_edges):
            likelihood = np.array([1.0 - self.p_error[i], self.p_error[i]])
            off = int(self.offset[i])
            if not edges:
                z_i = likelihood[off]
            else:
                parity_bias = float(np.prod(vf[edges, 0] - vf[edges, 1]))
                p_even = 0.5 * (1.0 + parity_bias)
                p_odd = 0.5 * (1.0 - parity_bias)
                z_i = p_even * likelihood[off] + p_odd * likelihood[off ^ 1]
            log_z += np.log(max(z_i, self.epsilon))

        for edges in self.col_edges:
            if not edges:
                z_j = 2.0
            else:
                z_j = float(np.prod(fv[edges, 0]) + np.prod(fv[edges, 1]))
            log_z += np.log(max(z_j, self.epsilon))

        for edge in range(len(self.edge_keys)):
            z_ij = float(np.dot(vf[edge], fv[edge]))
            log_z -= np.log(max(z_ij, self.epsilon))

        return float(log_z)

    def run(
        self,
        *,
        max_iterations: int = 50,
        damping: float = 0.0,
        tolerance: float = 1e-9,
        initial_message_state: Optional[Mapping[EdgeKey, Tuple[Array, Array]]] = None,
    ) -> BPResult:
        if max_iterations < 0:
            raise ValueError("max_iterations must be non-negative.")
        if not (0.0 <= damping < 1.0):
            raise ValueError("damping must lie in [0, 1).")

        vf, fv = self._initialize(initial_message_state)
        converged = False
        max_change = np.inf
        iterations_done = 0

        for iteration in range(1, max_iterations + 1):
            proposed_fv = self._factor_update(vf)
            new_fv = damping * fv + (1.0 - damping) * proposed_fv
            for edge in range(new_fv.shape[0]):
                new_fv[edge] = self._normalize(new_fv[edge])

            proposed_vf = self._variable_update(new_fv)
            new_vf = damping * vf + (1.0 - damping) * proposed_vf
            for edge in range(new_vf.shape[0]):
                new_vf[edge] = self._normalize(new_vf[edge])

            change_fv = float(np.max(np.abs(new_fv - fv))) if new_fv.size else 0.0
            change_vf = float(np.max(np.abs(new_vf - vf))) if new_vf.size else 0.0
            max_change = max(change_fv, change_vf)
            vf, fv = new_vf, new_fv
            iterations_done = iteration
            if max_change < tolerance:
                converged = True
                break

        state = {
            key: (vf[edge].copy(), fv[edge].copy())
            for edge, key in enumerate(self.edge_keys)
        }
        return BPResult(
            converged=converged,
            iterations=iterations_done,
            variable_prob_1=self._variable_beliefs(fv),
            log_z_bethe=self._bethe_log_partition(vf, fv),
            max_message_change=float(max_change),
            message_state=state,
        )


@dataclass
class SplitState:
    logical_bits: Array
    remaining: Tuple[int, ...]
    offset: Array
    path_score: float
    bp_result: Optional[BPResult]


@dataclass
class LeafResult:
    logical_label: Array
    correction: Array
    path_score: float
    log_z_bethe: float
    bp_converged: bool
    bp_iterations: int


@dataclass
class QuotientSplitResult:
    best_logical_label: Array
    correction: Array
    leaves: List[LeafResult]
    surviving_prefix_counts: List[int]


class QuotientSplitDecoder:
    """Quotient successive-cancellation list BP decoder."""

    UNFIXED = np.uint8(255)

    def __init__(
        self,
        d_matrix: Array,
        q_matrix: Array,
        x_s: Array,
        p_error: Array,
        *,
        list_width: int = 2,
        split_order: Optional[Sequence[int]] = None,
        adaptive_split: bool = False,
        internal_iterations: int = 5,
        leaf_iterations: int = 50,
        damping: float = 0.2,
        tolerance: float = 1e-9,
        probability_floor: float = 1e-300,
        validate: bool = True,
    ) -> None:
        self.d = as_binary_array(d_matrix, ndim=2)
        self.q = as_binary_array(q_matrix, ndim=2)
        self.x_s = as_binary_array(x_s, ndim=1)
        self.p_error = np.asarray(p_error, dtype=float)

        self.n, self.d_count = self.d.shape
        if self.q.shape[0] != self.n:
            raise ValueError("D and Q must have the same number of rows.")
        self.k = self.q.shape[1]
        if self.x_s.shape != (self.n,) or self.p_error.shape != (self.n,):
            raise ValueError("x_s and p_error have incompatible shapes.")
        if list_width < 1:
            raise ValueError("list_width must be at least 1.")

        self.list_width = int(list_width)
        self.internal_iterations = int(internal_iterations)
        self.leaf_iterations = int(leaf_iterations)
        self.damping = float(damping)
        self.tolerance = float(tolerance)
        self.probability_floor = float(probability_floor)
        self.adaptive_split = bool(adaptive_split)

        if split_order is None:
            self.split_order = tuple(range(self.k))
        else:
            order = tuple(int(j) for j in split_order)
            if sorted(order) != list(range(self.k)):
                raise ValueError("split_order must be a permutation of range(k).")
            self.split_order = order

        if validate:
            validate_quotient_decomposition(self.d, self.q)

        self.d_ids = tuple(range(self.d_count))
        self.q_ids = tuple(self.d_count + j for j in range(self.k))

    def _build_graph(self, state: SplitState):
        q_remaining = (
            self.q[:, state.remaining]
            if state.remaining
            else np.zeros((self.n, 0), dtype=np.uint8)
        )
        generator = np.concatenate([self.d, q_remaining], axis=1)
        coordinate_ids = list(self.d_ids) + [self.q_ids[j] for j in state.remaining]
        graph = GeneratorGraphBP(
            generator,
            state.offset,
            self.p_error,
            coordinate_ids=coordinate_ids,
        )
        return graph

    def _run_internal_bp(self, state: SplitState) -> SplitState:
        graph = self._build_graph(state)
        initial = state.bp_result.message_state if state.bp_result else None
        result = graph.run(
            max_iterations=self.internal_iterations,
            damping=self.damping,
            tolerance=self.tolerance,
            initial_message_state=initial,
        )
        return SplitState(
            logical_bits=state.logical_bits,
            remaining=state.remaining,
            offset=state.offset,
            path_score=state.path_score,
            bp_result=result,
        )

    def _choose_split_coordinate(self, active: Sequence[SplitState], depth: int) -> int:
        if not self.adaptive_split:
            return self.split_order[depth]

        common = set(active[0].remaining)
        for state in active[1:]:
            common.intersection_update(state.remaining)
        if not common:
            raise RuntimeError("No common remaining logical coordinate.")

        scores = np.array([state.path_score for state in active], dtype=float)
        scores -= np.max(scores)
        weights = np.exp(scores)
        weights /= weights.sum()

        best_j = None
        best_entropy = np.inf
        for j in sorted(common):
            avg_entropy = 0.0
            for weight, state in zip(weights, active):
                if state.bp_result is None:
                    raise RuntimeError("Missing BP result.")
                local = self.d_count + state.remaining.index(j)
                p1 = float(np.clip(
                    state.bp_result.variable_prob_1[local],
                    self.probability_floor,
                    1.0 - self.probability_floor,
                ))
                p0 = 1.0 - p1
                avg_entropy += weight * (-(p0 * np.log(p0) + p1 * np.log(p1)))
            if avg_entropy < best_entropy:
                best_entropy = avg_entropy
                best_j = j
        return int(best_j)

    def _branch(self, parent: SplitState, logical_index: int) -> List[SplitState]:
        if parent.bp_result is None:
            raise RuntimeError("Cannot branch without BP results.")
        local = self.d_count + parent.remaining.index(logical_index)
        p1 = float(np.clip(
            parent.bp_result.variable_prob_1[local],
            self.probability_floor,
            1.0 - self.probability_floor,
        ))
        probabilities = (1.0 - p1, p1)
        remaining = tuple(j for j in parent.remaining if j != logical_index)

        children = []
        for bit in (0, 1):
            bits = parent.logical_bits.copy()
            bits[logical_index] = bit
            offset = parent.offset.copy()
            if bit:
                offset ^= self.q[:, logical_index]
            children.append(SplitState(
                logical_bits=bits,
                remaining=remaining,
                offset=offset,
                path_score=parent.path_score + np.log(max(probabilities[bit], self.probability_floor)),
                bp_result=parent.bp_result,
            ))
        return children

    def _evaluate_leaf(self, state: SplitState) -> LeafResult:
        graph = GeneratorGraphBP(
            self.d,
            state.offset,
            self.p_error,
            coordinate_ids=self.d_ids,
        )
        initial = state.bp_result.message_state if state.bp_result else None
        result = graph.run(
            max_iterations=self.leaf_iterations,
            damping=self.damping,
            tolerance=self.tolerance,
            initial_message_state=initial,
        )
        return LeafResult(
            logical_label=state.logical_bits.copy(),
            correction=state.offset.copy(),
            path_score=float(state.path_score),
            log_z_bethe=float(result.log_z_bethe),
            bp_converged=result.converged,
            bp_iterations=result.iterations,
        )

    def decode(self) -> QuotientSplitResult:
        root = SplitState(
            logical_bits=np.full(self.k, self.UNFIXED, dtype=np.uint8),
            remaining=tuple(range(self.k)),
            offset=self.x_s.copy(),
            path_score=0.0,
            bp_result=None,
        )
        active = [self._run_internal_bp(root)]
        surviving_counts = [1]

        for depth in range(self.k):
            split_index = self._choose_split_coordinate(active, depth)
            candidates: List[SplitState] = []
            for state in active:
                candidates.extend(self._branch(state, split_index))
            candidates.sort(key=lambda state: state.path_score, reverse=True)
            active = candidates[: self.list_width]
            if depth < self.k - 1:
                active = [self._run_internal_bp(state) for state in active]
            surviving_counts.append(len(active))

        leaves = [self._evaluate_leaf(state) for state in active]
        leaves.sort(key=lambda leaf: leaf.log_z_bethe, reverse=True)
        best = leaves[0]
        return QuotientSplitResult(
            best_logical_label=best.logical_label.copy(),
            correction=best.correction.copy(),
            leaves=leaves,
            surviving_prefix_counts=surviving_counts,
        )


def exact_affine_partition(
    generator: Array,
    offset: Array,
    p_error: Array,
    *,
    max_latent_variables: int = 24,
) -> float:
    """Exact enumeration helper for small tests."""
    generator = as_binary_array(generator, ndim=2)
    offset = as_binary_array(offset, ndim=1)
    p_error = np.asarray(p_error, dtype=float)
    n, r = generator.shape
    if r > max_latent_variables:
        raise ValueError("Too many latent variables for exact enumeration.")

    total = 0.0
    for bits in product((0, 1), repeat=r):
        latent = np.fromiter(bits, dtype=np.uint8, count=r)
        x = offset ^ ((generator @ latent) & 1)
        total += float(np.prod(np.where(x, p_error, 1.0 - p_error)))
    return total


def exact_coset_masses(
    d_matrix: Array,
    q_matrix: Array,
    x_s: Array,
    p_error: Array,
    *,
    max_total_variables: int = 24,
):
    """Exact Z_ell values for every logical label on a small instance."""
    d_matrix = as_binary_array(d_matrix, ndim=2)
    q_matrix = as_binary_array(q_matrix, ndim=2)
    x_s = as_binary_array(x_s, ndim=1)
    k = q_matrix.shape[1]
    if d_matrix.shape[1] + k > max_total_variables:
        raise ValueError("Too many variables for exact enumeration.")

    masses = {}
    for label_bits in product((0, 1), repeat=k):
        label = np.fromiter(label_bits, dtype=np.uint8, count=k)
        offset = x_s ^ ((q_matrix @ label) & 1)
        masses[tuple(int(v) for v in label)] = exact_affine_partition(
            d_matrix,
            offset,
            p_error,
            max_latent_variables=max_total_variables,
        )
    return masses


if __name__ == "__main__":
    D = np.array([[1], [1], [0]], dtype=np.uint8)
    Q = np.array([[0], [1], [1]], dtype=np.uint8)
    x_s = np.array([1, 0, 0], dtype=np.uint8)
    p_error = np.array([0.08, 0.12, 0.05], dtype=float)

    decoder = QuotientSplitDecoder(
        D,
        Q,
        x_s,
        p_error,
        list_width=2,
        internal_iterations=10,
        leaf_iterations=50,
        damping=0.2,
    )
    result = decoder.decode()
    print("Best logical label:", result.best_logical_label.tolist())
    print("Canonical correction:", result.correction.tolist())
    for leaf in result.leaves:
        print(
            "leaf",
            leaf.logical_label.tolist(),
            "log Z_B =",
            leaf.log_z_bethe,
            "path score =",
            leaf.path_score,
        )
    print("Exact masses:", exact_coset_masses(D, Q, x_s, p_error))
