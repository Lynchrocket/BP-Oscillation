# quotient
"""Trajectory-screened quotient BP decoder.

This module extends the original quotient_bp.py with the two-stage algorithm:

    one H-BP run -> unique logical labels along the hard-decision trajectory
    -> D-BP only on those shortlisted logical cosets.

It also provides a reliability-frozen variant that keeps the final values of
reliable H-BP variables fixed and projects all syndrome-compatible flips of the
remaining variables onto their reachable logical cosets.

The exhaustive all-coset baseline is retained as
``decode_all_logical_cosets_bp``.  The main entry point
``decode_logical_coset_bp`` now implements the trajectory-screened algorithm.
"""

from dataclasses import dataclass
from itertools import product
from math import comb, log
from typing import Any

import numpy as np


def _gf2(A: np.ndarray) -> np.ndarray:
    return np.asarray(A, dtype=np.uint8) & 1


def gf2_rref(
    A: np.ndarray,
    *,
    return_transform: bool = False,
) -> tuple[np.ndarray, list[int]] | tuple[np.ndarray, list[int], np.ndarray]:
    """Reduced row-echelon form over GF(2).

    If return_transform=True, also return P such that R = P @ A (mod 2).
    """
    R = _gf2(A).copy()
    m, n = R.shape
    P = np.eye(m, dtype=np.uint8)
    pivots: list[int] = []
    row = 0

    for col in range(n):
        candidates = np.flatnonzero(R[row:, col])
        if candidates.size == 0:
            continue

        pivot_row = row + int(candidates[0])
        if pivot_row != row:
            R[[row, pivot_row]] = R[[pivot_row, row]]
            P[[row, pivot_row]] = P[[pivot_row, row]]

        for r in range(m):
            if r != row and R[r, col]:
                R[r] ^= R[row]
                P[r] ^= P[row]

        pivots.append(col)
        row += 1
        if row == m:
            break

    if return_transform:
        return R, pivots, P
    return R, pivots


def gf2_rank(A: np.ndarray) -> int:
    return len(gf2_rref(A)[1])


def gf2_nullspace(A: np.ndarray) -> np.ndarray:
    """Return a matrix whose columns form a basis of ker(A) over GF(2)."""
    R, pivots = gf2_rref(A)
    _, n = R.shape
    pivot_set = set(pivots)
    free_cols = [j for j in range(n) if j not in pivot_set]

    N = np.zeros((n, len(free_cols)), dtype=np.uint8)
    for basis_index, free_col in enumerate(free_cols):
        N[free_col, basis_index] = 1
        for pivot_row, pivot_col in enumerate(pivots):
            # In GF(2), subtraction and addition are identical.
            N[pivot_col, basis_index] = R[pivot_row, free_col]
    return N


def gf2_inverse(A: np.ndarray) -> np.ndarray:
    """Inverse of a square nonsingular matrix over GF(2)."""
    A = _gf2(A)
    n, m = A.shape
    if n != m:
        raise ValueError("A must be square.")

    aug = np.concatenate([A.copy(), np.eye(n, dtype=np.uint8)], axis=1)
    row = 0
    for col in range(n):
        candidates = np.flatnonzero(aug[row:, col])
        if candidates.size == 0:
            raise ValueError("A is singular over GF(2).")

        pivot_row = row + int(candidates[0])
        if pivot_row != row:
            aug[[row, pivot_row]] = aug[[pivot_row, row]]

        for r in range(n):
            if r != row and aug[r, col]:
                aug[r] ^= aug[row]
        row += 1

    return aug[:, n:]


@dataclass(frozen=True)
class QuotientModel:
    H: np.ndarray
    L: np.ndarray

    # Particular-solution map: x_s = T s and H T H = H.
    T: np.ndarray

    # Basis of ker(H).
    K: np.ndarray

    # Basis of ker(H) ∩ ker(L): H D = 0 and L D = 0.
    D: np.ndarray

    # Syndrome-null logical directions: H Q = 0.
    Q: np.ndarray

    # L Q. If all observable directions are reachable, this is I.
    logical_map: np.ndarray

    rank_H: int
    logical_rank: int

    @property
    def num_error_mechanisms(self) -> int:
        return self.H.shape[1]

    @property
    def num_observables(self) -> int:
        return self.L.shape[0]

    @property
    def degeneracy_dimension(self) -> int:
        return self.D.shape[1]

    @property
    def canonical_logical_coordinates(self) -> bool:
        return (
            self.logical_rank == self.num_observables
            and np.array_equal(
                self.logical_map,
                np.eye(self.num_observables, dtype=np.uint8),
            )
        )


def build_quotient_model(H: np.ndarray, L: np.ndarray) -> QuotientModel:
    """Construct T, K, D and Q from detector matrix H and observable matrix L.

    Every syndrome-compatible configuration can be represented as

        x = x_s ⊕ D u ⊕ Q λ,

    when the returned columns of [D Q] span ker(H).
    """
    H = _gf2(H)
    L = _gf2(L)

    m, n = H.shape
    k, n_L = L.shape
    if n != n_L:
        raise ValueError(
            f"H and L must have the same number of columns; got {n} and {n_L}."
        )

    # Build a generalized right inverse T using the row operations of H.
    # For every legal syndrome s in im(H), H(Ts) = s.
    _, pivot_columns_H, row_transform = gf2_rref(
        H,
        return_transform=True,
    )
    rank_H = len(pivot_columns_H)

    T = np.zeros((n, m), dtype=np.uint8)
    for pivot_row, pivot_col in enumerate(pivot_columns_H):
        T[pivot_col, :] = row_transform[pivot_row, :]

    # Kernel coordinates.
    K = gf2_nullspace(H)
    M = (L @ K) & 1

    # Degeneracy directions: z = K a with M a = 0.
    A = gf2_nullspace(M)
    D = (K @ A) & 1

    # Select independent logical directions from columns of K.
    _, logical_pivots = gf2_rref(M)
    logical_rank = len(logical_pivots)

    if logical_rank == 0:
        Q = np.zeros((n, 0), dtype=np.uint8)
        logical_map = np.zeros((k, 0), dtype=np.uint8)
    else:
        Q0 = K[:, logical_pivots]
        R = (L @ Q0) & 1

        # If every observable bit is independently reachable inside ker(H),
        # change basis so that LQ = I.
        if logical_rank == k:
            Q = (Q0 @ gf2_inverse(R)) & 1
            logical_map = np.eye(k, dtype=np.uint8)
        else:
            Q = Q0
            logical_map = R

    model = QuotientModel(
        H=H,
        L=L,
        T=T,
        K=K,
        D=D,
        Q=Q,
        logical_map=logical_map,
        rank_H=rank_H,
        logical_rank=logical_rank,
    )

    # Internal consistency checks.
    assert np.array_equal((H @ T @ H) & 1, H)
    assert not np.any((H @ K) & 1)
    assert not np.any((H @ D) & 1)
    assert not np.any((L @ D) & 1)
    assert not np.any((H @ Q) & 1)
    assert np.array_equal((L @ Q) & 1, logical_map)

    expected_kernel_dimension = n - rank_H
    if D.shape[1] + Q.shape[1] != expected_kernel_dimension:
        raise RuntimeError(
            "The decomposition did not span ker(H). "
            "Check whether L and H were supplied with the intended orientation."
        )

    return model


def particular_solution(model: QuotientModel, syndrome: np.ndarray) -> np.ndarray:
    syndrome = _gf2(syndrome).reshape(-1)
    if syndrome.size != model.H.shape[0]:
        raise ValueError("The syndrome length must equal H.shape[0].")

    x_s = (model.T @ syndrome) & 1
    if not np.array_equal((model.H @ x_s) & 1, syndrome):
        raise ValueError("The supplied syndrome is not in im(H).")
    return x_s


def quotient_graph_stats(model: QuotientModel) -> dict[str, Any]:
    """Basic graph diagnostics for the u-variable/factor graph defined by D."""
    D = model.D.astype(np.int64)
    factor_weights = D.sum(axis=1)
    variable_weights = D.sum(axis=0)

    overlap = D.T @ D
    four_cycles = 0
    for a in range(overlap.shape[0]):
        for b in range(a + 1, overlap.shape[1]):
            shared_factors = int(overlap[a, b])
            if shared_factors >= 2:
                four_cycles += comb(shared_factors, 2)

    return {
        "num_factors": int(D.shape[0]),
        "num_degeneracy_variables": int(D.shape[1]),
        "num_edges": int(D.sum()),
        "max_factor_degree": int(factor_weights.max(initial=0)),
        "max_variable_degree": int(variable_weights.max(initial=0)),
        "mean_factor_degree": float(factor_weights.mean()) if D.shape[0] else 0.0,
        "mean_variable_degree": float(variable_weights.mean()) if D.shape[1] else 0.0,
        "four_cycle_count": int(four_cycles),
    }


@dataclass(frozen=True)
class BPPartitionResult:
    log_partition: float
    iterations: int
    converged: bool
    max_message_change: float
    u_marginals: np.ndarray


def quotient_bp_log_partition(
    D: np.ndarray,
    offset: np.ndarray,
    error_probabilities: np.ndarray,
    *,
    max_iter: int = 100,
    damping: float = 0.2,
    tolerance: float = 1e-10,
) -> BPPartitionResult:
    """Approximate log sum_u P(offset ⊕ D u) using sum-product BP.

    Factor i has potential

        ψ_i(u_N(i)) = P_i(offset_i ⊕ XOR_{a in N(i)} u_a).

    The returned log partition is the Bethe approximation at the final
    message state.
    """
    D = _gf2(D)
    offset = _gf2(offset).reshape(-1)
    p = np.asarray(error_probabilities, dtype=float).reshape(-1)

    n, d = D.shape
    if offset.size != n or p.size != n:
        raise ValueError("offset and error_probabilities must have length D.shape[0].")
    if not 0.0 <= damping < 1.0:
        raise ValueError("damping must lie in [0, 1).")

    eps = 1e-15
    p = np.clip(p, eps, 1.0 - eps)

    factor_neighbors = [
        np.flatnonzero(D[i]).astype(int).tolist()
        for i in range(n)
    ]
    variable_neighbors = [
        np.flatnonzero(D[:, a]).astype(int).tolist()
        for a in range(d)
    ]

    # Messages are normalized two-component probability vectors.
    v_to_f = {
        (i, a): np.array([0.5, 0.5], dtype=float)
        for i in range(n)
        for a in factor_neighbors[i]
    }
    f_to_v = {
        (i, a): np.array([0.5, 0.5], dtype=float)
        for i in range(n)
        for a in factor_neighbors[i]
    }

    # ψ_i(parity) = 1/2 [1 + c_i (-1)^parity].
    c = (1.0 - 2.0 * p) * np.where(offset == 0, 1.0, -1.0)

    converged = False
    max_change = float("inf")

    for iteration in range(1, max_iter + 1):
        max_change = 0.0

        new_f_to_v: dict[tuple[int, int], np.ndarray] = {}
        for i, neighbors in enumerate(factor_neighbors):
            if not neighbors:
                continue

            magnetization = {
                a: float(v_to_f[(i, a)][0] - v_to_f[(i, a)][1])
                for a in neighbors
            }

            for a in neighbors:
                product_other = 1.0
                for b in neighbors:
                    if b != a:
                        product_other *= magnetization[b]

                q = float(np.clip(c[i] * product_other, -1.0 + eps, 1.0 - eps))
                proposed = np.array(
                    [(1.0 + q) / 2.0, (1.0 - q) / 2.0],
                    dtype=float,
                )

                old = f_to_v[(i, a)]
                message = (1.0 - damping) * proposed + damping * old
                message /= message.sum()

                new_f_to_v[(i, a)] = message
                max_change = max(
                    max_change,
                    float(np.max(np.abs(message - old))),
                )

        f_to_v = new_f_to_v

        new_v_to_f: dict[tuple[int, int], np.ndarray] = {}
        for a, neighbors in enumerate(variable_neighbors):
            for i in neighbors:
                proposed = np.ones(2, dtype=float)
                for other_factor in neighbors:
                    if other_factor != i:
                        proposed *= f_to_v[(other_factor, a)]

                total = proposed.sum()
                if total <= 0.0:
                    proposed[:] = 0.5
                else:
                    proposed /= total

                old = v_to_f[(i, a)]
                message = (1.0 - damping) * proposed + damping * old
                message /= message.sum()

                new_v_to_f[(i, a)] = message
                max_change = max(
                    max_change,
                    float(np.max(np.abs(message - old))),
                )

        v_to_f = new_v_to_f

        if max_change < tolerance:
            converged = True
            break

    # Bethe log-partition:
    # log Z_B = Σ_f log Z_f + Σ_v log Z_v - Σ_edges log Z_edge.
    log_Z = 0.0

    for i, neighbors in enumerate(factor_neighbors):
        if not neighbors:
            Z_factor = (1.0 - p[i]) if offset[i] == 0 else p[i]
        else:
            product_magnetization = 1.0
            for a in neighbors:
                message = v_to_f[(i, a)]
                product_magnetization *= float(message[0] - message[1])
            Z_factor = 0.5 * (1.0 + c[i] * product_magnetization)

        if Z_factor <= 0.0:
            log_Z = float("-inf")
            break
        log_Z += log(Z_factor)

    u_marginals = np.empty((d, 2), dtype=float)
    if np.isfinite(log_Z):
        for a, neighbors in enumerate(variable_neighbors):
            belief = np.ones(2, dtype=float)
            for i in neighbors:
                belief *= f_to_v[(i, a)]

            Z_variable = float(belief.sum())
            if Z_variable <= 0.0:
                log_Z = float("-inf")
                break

            log_Z += log(Z_variable)
            u_marginals[a] = belief / Z_variable

    if np.isfinite(log_Z):
        for i, neighbors in enumerate(factor_neighbors):
            for a in neighbors:
                Z_edge = float(np.dot(v_to_f[(i, a)], f_to_v[(i, a)]))
                if Z_edge <= 0.0:
                    log_Z = float("-inf")
                    break
                log_Z -= log(Z_edge)

    return BPPartitionResult(
        log_partition=float(log_Z),
        iterations=iteration,
        converged=converged,
        max_message_change=float(max_change),
        u_marginals=u_marginals,
    )


def exact_log_partition(
    D: np.ndarray,
    offset: np.ndarray,
    error_probabilities: np.ndarray,
    *,
    max_dimension: int = 24,
) -> float:
    """Exact log sum_u P(offset ⊕ D u), for small degeneracy dimension."""
    D = _gf2(D)
    offset = _gf2(offset).reshape(-1)
    p = np.asarray(error_probabilities, dtype=float).reshape(-1)

    n, d = D.shape
    if d > max_dimension:
        raise ValueError(
            f"Exact enumeration requires 2^{d} states; "
            f"increase max_dimension explicitly if intended."
        )

    total = 0.0
    for bits in product((0, 1), repeat=d):
        u = np.fromiter(bits, dtype=np.uint8, count=d)
        x = offset ^ ((D @ u) & 1)
        total += float(np.prod(np.where(x == 1, p, 1.0 - p)))

    return float(np.log(total))



@dataclass(frozen=True)
class HBPTraceResult:
    """Result of the first-stage BP on the original H Tanner graph.

    belief_history[t, i] is the posterior LLR of error variable i after
    iteration t + 1. A negative LLR gives hard decision 1.
    """

    posterior_llrs: np.ndarray
    belief_history: np.ndarray
    hard_decision_history: np.ndarray
    syndrome_consistent_history: np.ndarray
    residual_syndrome_weights: np.ndarray
    iterations: int
    converged: bool
    max_message_change: float


def parity_check_bp_trace(
    H: np.ndarray,
    syndrome: np.ndarray,
    error_probabilities: np.ndarray,
    *,
    max_iter: int = 50,
    damping: float = 0.2,
    tolerance: float = 1e-10,
    method: str = "sum_product",
    min_sum_scaling: float = 1.0,
    llr_clip: float = 50.0,
    stop_if_syndrome_satisfied: bool = False,
) -> HBPTraceResult:
    """Run BP on H and record the hard decision from every iteration.

    Parameters
    ----------
    H
        Detector-by-error parity-check matrix.
    syndrome
        Target detector syndrome.
    error_probabilities
        Independent probability of each elementary error mechanism.
    method
        ``"sum_product"`` or ``"min_sum"``.  D-BP remains sum-product,
        because its purpose is to estimate a coset partition function.

    Notes
    -----
    Damping uses

        message <- (1 - damping) * proposed + damping * old.
    """
    H = _gf2(H)
    syndrome = _gf2(syndrome).reshape(-1)
    p = np.asarray(error_probabilities, dtype=float).reshape(-1)

    m, n = H.shape
    if syndrome.size != m:
        raise ValueError("syndrome must have length H.shape[0].")
    if p.size != n:
        raise ValueError("error_probabilities must have length H.shape[1].")
    if max_iter < 1:
        raise ValueError("max_iter must be at least 1.")
    if not 0.0 <= damping < 1.0:
        raise ValueError("damping must lie in [0, 1).")
    if method not in {"sum_product", "min_sum"}:
        raise ValueError("method must be 'sum_product' or 'min_sum'.")
    if min_sum_scaling <= 0.0:
        raise ValueError("min_sum_scaling must be positive.")

    eps = 1e-15
    p = np.clip(p, eps, 1.0 - eps)
    prior_llr = np.log((1.0 - p) / p)
    prior_llr = np.clip(prior_llr, -llr_clip, llr_clip)

    edge_checks, edge_variables = np.nonzero(H)
    edge_checks = edge_checks.astype(int)
    edge_variables = edge_variables.astype(int)
    num_edges = edge_checks.size

    check_edges = [np.flatnonzero(edge_checks == a) for a in range(m)]
    variable_edges = [np.flatnonzero(edge_variables == i) for i in range(n)]

    # LLR messages. Variable-to-check starts from the channel prior.
    v_to_c = prior_llr[edge_variables].astype(float, copy=True)
    c_to_v = np.zeros(num_edges, dtype=float)

    belief_history: list[np.ndarray] = []
    hard_history: list[np.ndarray] = []
    valid_history: list[bool] = []
    residual_weight_history: list[int] = []

    converged = False
    max_change = float("inf")

    for iteration in range(1, max_iter + 1):
        proposed_c_to_v = np.zeros_like(c_to_v)

        for a, eids in enumerate(check_edges):
            degree = eids.size
            if degree == 0:
                if syndrome[a]:
                    raise ValueError(
                        "The syndrome is inconsistent: a zero-degree check is triggered."
                    )
                continue

            incoming = v_to_c[eids]
            syndrome_sign = -1.0 if syndrome[a] else 1.0

            if method == "sum_product":
                tanh_values = np.tanh(np.clip(incoming / 2.0, -20.0, 20.0))
                prefix = np.ones(degree + 1, dtype=float)
                suffix = np.ones(degree + 1, dtype=float)
                for j in range(degree):
                    prefix[j + 1] = prefix[j] * tanh_values[j]
                for j in range(degree - 1, -1, -1):
                    suffix[j] = suffix[j + 1] * tanh_values[j]

                for j, edge_id in enumerate(eids):
                    product_other = syndrome_sign * prefix[j] * suffix[j + 1]
                    product_other = float(
                        np.clip(product_other, -1.0 + eps, 1.0 - eps)
                    )
                    proposed_c_to_v[edge_id] = np.clip(
                        2.0 * np.arctanh(product_other),
                        -llr_clip,
                        llr_clip,
                    )
            else:
                signs = np.where(incoming < 0.0, -1.0, 1.0)
                magnitudes = np.abs(incoming)
                total_sign = syndrome_sign * float(np.prod(signs))

                if degree == 1:
                    min1 = llr_clip
                    min2 = llr_clip
                    min1_index = 0
                else:
                    order = np.argsort(magnitudes)
                    min1_index = int(order[0])
                    min1 = float(magnitudes[min1_index])
                    min2 = float(magnitudes[int(order[1])])

                for j, edge_id in enumerate(eids):
                    magnitude = min2 if j == min1_index else min1
                    outgoing_sign = total_sign * signs[j]
                    proposed_c_to_v[edge_id] = np.clip(
                        min_sum_scaling * outgoing_sign * magnitude,
                        -llr_clip,
                        llr_clip,
                    )

        new_c_to_v = (1.0 - damping) * proposed_c_to_v + damping * c_to_v
        new_c_to_v = np.clip(new_c_to_v, -llr_clip, llr_clip)

        posterior_llr = prior_llr.copy()
        for i, eids in enumerate(variable_edges):
            if eids.size:
                posterior_llr[i] += float(np.sum(new_c_to_v[eids]))
        posterior_llr = np.clip(posterior_llr, -llr_clip, llr_clip)

        proposed_v_to_c = np.empty_like(v_to_c)
        for i, eids in enumerate(variable_edges):
            if not eids.size:
                continue
            total = posterior_llr[i]
            proposed_v_to_c[eids] = total - new_c_to_v[eids]

        new_v_to_c = (1.0 - damping) * proposed_v_to_c + damping * v_to_c
        new_v_to_c = np.clip(new_v_to_c, -llr_clip, llr_clip)

        change_c = float(np.max(np.abs(new_c_to_v - c_to_v), initial=0.0))
        change_v = float(np.max(np.abs(new_v_to_c - v_to_c), initial=0.0))
        max_change = max(change_c, change_v)

        c_to_v = new_c_to_v
        v_to_c = new_v_to_c

        hard = (posterior_llr < 0.0).astype(np.uint8)
        residual = syndrome ^ ((H @ hard) & 1)
        syndrome_consistent = not np.any(residual)

        belief_history.append(posterior_llr.copy())
        hard_history.append(hard)
        valid_history.append(syndrome_consistent)
        residual_weight_history.append(int(residual.sum()))

        if max_change < tolerance:
            converged = True
            break
        if stop_if_syndrome_satisfied and syndrome_consistent:
            break

    return HBPTraceResult(
        posterior_llrs=belief_history[-1].copy(),
        belief_history=np.asarray(belief_history, dtype=float),
        hard_decision_history=np.asarray(hard_history, dtype=np.uint8),
        syndrome_consistent_history=np.asarray(valid_history, dtype=bool),
        residual_syndrome_weights=np.asarray(residual_weight_history, dtype=int),
        iterations=iteration,
        converged=converged,
        max_message_change=float(max_change),
    )


def logical_neutral_repair_map(model: QuotientModel) -> np.ndarray:
    """Return T0 with H T0 H = H and L T0 = 0, when LQ = I.

    Repairing a hard decision by T0 times its residual syndrome therefore
    makes it syndrome-valid without changing its observable label.
    """
    if not model.canonical_logical_coordinates:
        raise ValueError(
            "A logical-neutral repair map requires canonical logical coordinates "
            "with LQ = I."
        )

    LT = (model.L @ model.T) & 1
    T0 = model.T ^ ((model.Q @ LT) & 1)

    assert np.array_equal((model.H @ T0 @ model.H) & 1, model.H)
    assert not np.any((model.L @ T0) & 1)
    return T0


def _label_key(label: np.ndarray) -> tuple[int, ...]:
    return tuple(int(x) for x in _gf2(label).reshape(-1))


def _logical_coordinate_for_label(
    model: QuotientModel,
    x_s: np.ndarray,
    observable_label: np.ndarray,
) -> np.ndarray:
    """Convert an observable label into the Q-coordinate relative to x_s."""
    observable_label = _gf2(observable_label).reshape(-1)
    base_label = (model.L @ x_s) & 1
    delta = observable_label ^ base_label

    if model.canonical_logical_coordinates:
        return delta

    # Generic one-solution GF(2) solve for logical_map @ coordinate = delta.
    A = model.logical_map
    aug = np.concatenate([A, delta[:, None]], axis=1)
    R, pivots = gf2_rref(aug)
    num_variables = A.shape[1]

    for row in range(R.shape[0]):
        if not np.any(R[row, :num_variables]) and R[row, num_variables]:
            raise ValueError("The observable label is unreachable for this syndrome.")

    coordinate = np.zeros(num_variables, dtype=np.uint8)
    for row, pivot_col in enumerate(pivots):
        if pivot_col < num_variables:
            coordinate[pivot_col] = R[row, num_variables]
    return coordinate


def extract_trajectory_cosets(
    model: QuotientModel,
    syndrome: np.ndarray,
    h_bp: HBPTraceResult,
    *,
    burn_in: int = 0,
    stride: int = 1,
    only_syndrome_consistent: bool = False,
    repair_hard_decisions: bool = True,
    max_candidates: int | None = None,
) -> list[dict[str, Any]]:
    """Collect unique logical labels visited by H-BP hard decisions.

    Unique labels are preserved in first-appearance order.  When
    ``repair_hard_decisions`` is true, each trajectory hard decision is made
    syndrome-valid before its metadata representative is stored.  For the
    usual full-rank observable case, the repair is logical-neutral, so the
    label is exactly the raw label L @ hard_decision.
    """
    syndrome = _gf2(syndrome).reshape(-1)
    if stride < 1:
        raise ValueError("stride must be at least 1.")
    if burn_in < 0:
        raise ValueError("burn_in must be nonnegative.")

    x_s = particular_solution(model, syndrome)
    T0: np.ndarray | None = None
    if repair_hard_decisions and model.canonical_logical_coordinates:
        T0 = logical_neutral_repair_map(model)

    candidates_by_key: dict[tuple[int, ...], dict[str, Any]] = {}

    indices = range(burn_in, h_bp.iterations, stride)
    for index in indices:
        hard = h_bp.hard_decision_history[index].copy()
        raw_valid = bool(h_bp.syndrome_consistent_history[index])
        residual = syndrome ^ ((model.H @ hard) & 1)

        if only_syndrome_consistent and not raw_valid:
            continue

        repaired = hard
        if repair_hard_decisions and np.any(residual):
            if T0 is not None:
                repaired = hard ^ ((T0 @ residual) & 1)
            else:
                # Fallback for noncanonical logical maps.  This repair may
                # change the label, but always produces a legal representative.
                repaired = hard ^ ((model.T @ residual) & 1)

        if repair_hard_decisions:
            if not np.array_equal((model.H @ repaired) & 1, syndrome):
                raise RuntimeError("Internal error: trajectory repair failed.")
            label_source = repaired
        else:
            # The raw hard decision may violate the syndrome.  Its logical
            # label is still used for screening; a canonical syndrome-valid
            # representative is constructed below for D-BP.
            label_source = hard

        label = (model.L @ label_source) & 1
        key = _label_key(label)

        try:
            coordinate = _logical_coordinate_for_label(model, x_s, label)
        except ValueError:
            # This can only occur for a noncanonical/incomplete logical map.
            continue

        canonical_representative = x_s ^ ((model.Q @ coordinate) & 1)
        if not np.array_equal((model.H @ canonical_representative) & 1, syndrome):
            raise RuntimeError("Internal error: canonical representative has wrong syndrome.")
        if not np.array_equal((model.L @ canonical_representative) & 1, label):
            raise RuntimeError("Internal error: canonical representative has wrong label.")

        if key not in candidates_by_key:
            candidates_by_key[key] = {
                "observable_label": label.copy(),
                "logical_coordinate": coordinate.copy(),
                "representative": canonical_representative.copy(),
                "first_iteration": index + 1,
                "last_iteration": index + 1,
                "visit_count": 1,
                "best_residual_syndrome_weight": int(
                    h_bp.residual_syndrome_weights[index]
                ),
                "ever_syndrome_consistent": raw_valid,
                "trajectory_iterations": [index + 1],
                "example_repaired_hard_decision": repaired.copy(),
            }
        else:
            item = candidates_by_key[key]
            item["last_iteration"] = index + 1
            item["visit_count"] += 1
            item["best_residual_syndrome_weight"] = min(
                item["best_residual_syndrome_weight"],
                int(h_bp.residual_syndrome_weights[index]),
            )
            item["ever_syndrome_consistent"] = bool(
                item["ever_syndrome_consistent"] or raw_valid
            )
            item["trajectory_iterations"].append(index + 1)

    candidates = list(candidates_by_key.values())

    # A restrictive filter can remove everything.  Fall back to the final
    # trajectory state so the second stage is always well-defined.
    if not candidates:
        hard = h_bp.hard_decision_history[-1]
        residual = syndrome ^ ((model.H @ hard) & 1)
        repaired = hard
        if repair_hard_decisions and np.any(residual):
            if T0 is not None:
                repaired = hard ^ ((T0 @ residual) & 1)
            else:
                repaired = hard ^ ((model.T @ residual) & 1)
        label_source = repaired if repair_hard_decisions else hard
        label = (model.L @ label_source) & 1
        coordinate = _logical_coordinate_for_label(model, x_s, label)
        representative = x_s ^ ((model.Q @ coordinate) & 1)
        candidates = [{
            "observable_label": label.copy(),
            "logical_coordinate": coordinate.copy(),
            "representative": representative.copy(),
            "first_iteration": h_bp.iterations,
            "last_iteration": h_bp.iterations,
            "visit_count": 1,
            "best_residual_syndrome_weight": int(
                h_bp.residual_syndrome_weights[-1]
            ),
            "ever_syndrome_consistent": bool(h_bp.syndrome_consistent_history[-1]),
            "trajectory_iterations": [h_bp.iterations],
            "example_repaired_hard_decision": repaired.copy(),
        }]

    if max_candidates is not None:
        if max_candidates < 1:
            raise ValueError("max_candidates must be at least 1 when supplied.")
        # Prefer candidates that came closest to satisfying the syndrome,
        # then those seen later, then those visited more often.
        candidates = sorted(
            candidates,
            key=lambda item: (
                item["best_residual_syndrome_weight"],
                -item["last_iteration"],
                -item["visit_count"],
            ),
        )[:max_candidates]

    return candidates

def decode_all_logical_cosets_bp(
    model: QuotientModel,
    syndrome: np.ndarray,
    error_probabilities: np.ndarray,
    *,
    max_iter: int = 100,
    damping: float = 0.2,
    tolerance: float = 1e-10,
) -> dict[str, Any]:
    """Baseline: approximate DQMLD by comparing every logical hypothesis.

    For each λ, BP approximates

        Z_λ(s) = sum_u P(x_s ⊕ D u ⊕ Q λ).

    Since L usually has only a few rows, enumerating 2^logical_rank hypotheses
    is cheap. For two independent observable bits, this means four BP runs.
    """
    x_s = particular_solution(model, syndrome)
    p = np.asarray(error_probabilities, dtype=float).reshape(-1)

    if p.size != model.num_error_mechanisms:
        raise ValueError("error_probabilities must have length H.shape[1].")

    hypotheses: list[dict[str, Any]] = []

    for bits in product((0, 1), repeat=model.logical_rank):
        logical_coordinate = np.fromiter(
            bits,
            dtype=np.uint8,
            count=model.logical_rank,
        )

        offset = x_s ^ ((model.Q @ logical_coordinate) & 1)
        bp_result = quotient_bp_log_partition(
            model.D,
            offset,
            p,
            max_iter=max_iter,
            damping=damping,
            tolerance=tolerance,
        )

        observable_label = (model.L @ offset) & 1
        hypotheses.append(
            {
                "logical_coordinate": logical_coordinate,
                "observable_label": observable_label,
                "representative": offset,
                "log_partition": bp_result.log_partition,
                "bp": bp_result,
            }
        )

    log_values = np.array(
        [item["log_partition"] for item in hypotheses],
        dtype=float,
    )
    best_index = int(np.argmax(log_values))
    best = hypotheses[best_index]

    # Normalize the approximate coset masses for convenient diagnostics.
    finite = np.isfinite(log_values)
    posterior = np.zeros_like(log_values)
    if finite.any():
        shifted = log_values[finite] - np.max(log_values[finite])
        weights = np.exp(shifted)
        posterior[finite] = weights / weights.sum()

    for item, probability in zip(hypotheses, posterior):
        item["approximate_posterior"] = float(probability)

    correction = best["representative"].copy()

    # The representative automatically satisfies the supplied syndrome.
    assert np.array_equal((model.H @ correction) & 1, _gf2(syndrome).reshape(-1))

    return {
        "correction": correction,
        "chosen_observable_label": best["observable_label"].copy(),
        "chosen_logical_coordinate": best["logical_coordinate"].copy(),
        "chosen_log_partition": float(best["log_partition"]),
        "particular_solution": x_s,
        "hypotheses": hypotheses,
    }


def decode_logical_coset_bp(
    model: QuotientModel,
    syndrome: np.ndarray,
    error_probabilities: np.ndarray,
    *,
    h_max_iter: int = 50,
    h_damping: float = 0.2,
    h_tolerance: float = 1e-10,
    h_method: str = "sum_product",
    h_min_sum_scaling: float = 1.0,
    h_llr_clip: float = 50.0,
    h_stop_if_syndrome_satisfied: bool = False,
    trajectory_burn_in: int = 0,
    trajectory_stride: int = 1,
    trajectory_only_syndrome_consistent: bool = False,
    trajectory_repair_hard_decisions: bool = True,
    max_candidate_cosets: int | None = None,
    d_max_iter: int = 100,
    d_damping: float = 0.2,
    d_tolerance: float = 1e-10,
) -> dict[str, Any]:
    """H-BP trajectory screening followed by D-BP coset reranking.

    Stage 1
    -------
    Run BP once on the original H Tanner graph.  Convert the hard decision at
    every retained iteration into an observable label and keep the unique
    labels visited by the trajectory.

    Stage 2
    -------
    For each shortlisted label ell, construct a representative x_(s, ell) and
    run sum-product BP on D to approximate

        Z_ell(s) = sum_u P(x_(s, ell) xor D u).

    The returned posterior is normalized only over the shortlisted labels; it
    is not a calibrated posterior over logical labels that were never visited.
    """
    syndrome = _gf2(syndrome).reshape(-1)
    p = np.asarray(error_probabilities, dtype=float).reshape(-1)
    if p.size != model.num_error_mechanisms:
        raise ValueError("error_probabilities must have length H.shape[1].")

    h_bp = parity_check_bp_trace(
        model.H,
        syndrome,
        p,
        max_iter=h_max_iter,
        damping=h_damping,
        tolerance=h_tolerance,
        method=h_method,
        min_sum_scaling=h_min_sum_scaling,
        llr_clip=h_llr_clip,
        stop_if_syndrome_satisfied=h_stop_if_syndrome_satisfied,
    )

    candidates = extract_trajectory_cosets(
        model,
        syndrome,
        h_bp,
        burn_in=trajectory_burn_in,
        stride=trajectory_stride,
        only_syndrome_consistent=trajectory_only_syndrome_consistent,
        repair_hard_decisions=trajectory_repair_hard_decisions,
        max_candidates=max_candidate_cosets,
    )

    hypotheses: list[dict[str, Any]] = []
    for candidate in candidates:
        offset = candidate["representative"]
        bp_result = quotient_bp_log_partition(
            model.D,
            offset,
            p,
            max_iter=d_max_iter,
            damping=d_damping,
            tolerance=d_tolerance,
        )

        hypothesis = dict(candidate)
        hypothesis["log_partition"] = float(bp_result.log_partition)
        hypothesis["bp"] = bp_result
        hypotheses.append(hypothesis)

    log_values = np.asarray(
        [item["log_partition"] for item in hypotheses],
        dtype=float,
    )
    best_index = int(np.argmax(log_values))
    best = hypotheses[best_index]

    finite = np.isfinite(log_values)
    shortlist_posterior = np.zeros_like(log_values)
    if finite.any():
        shifted = log_values[finite] - np.max(log_values[finite])
        weights = np.exp(shifted)
        shortlist_posterior[finite] = weights / weights.sum()

    for item, probability in zip(hypotheses, shortlist_posterior):
        item["shortlist_posterior"] = float(probability)

    correction = best["representative"].copy()
    assert np.array_equal((model.H @ correction) & 1, syndrome)
    assert np.array_equal(
        (model.L @ correction) & 1,
        best["observable_label"],
    )

    trajectory_labels = np.asarray(
        [
            (model.L @ hard) & 1
            for hard in h_bp.hard_decision_history
        ],
        dtype=np.uint8,
    )

    return {
        "correction": correction,
        "chosen_observable_label": best["observable_label"].copy(),
        "chosen_logical_coordinate": best["logical_coordinate"].copy(),
        "chosen_log_partition": float(best["log_partition"]),
        "particular_solution": particular_solution(model, syndrome),
        "h_bp": h_bp,
        "trajectory_observable_labels_raw": trajectory_labels,
        "candidate_cosets": candidates,
        "candidate_count": len(candidates),
        "hypotheses": hypotheses,
        "screening": {
            "num_h_bp_iterations": h_bp.iterations,
            "num_unique_candidate_cosets": len(candidates),
            "used_all_logical_cosets": False,
            "posterior_is_shortlist_normalized": True,
        },
    }


def decode_reliability_frozen_quotient_bp(
    model: QuotientModel,
    syndrome: np.ndarray,
    error_probabilities: np.ndarray,
    *,
    reliability_threshold: float | None = None,
    max_unfrozen_variables: int | None = None,
    h_max_iter: int = 50,
    h_damping: float = 0.2,
    h_tolerance: float = 1e-10,
    h_method: str = "sum_product",
    h_min_sum_scaling: float = 1.0,
    h_llr_clip: float = 50.0,
    h_stop_if_syndrome_satisfied: bool = False,
    d_max_iter: int = 100,
    d_damping: float = 0.2,
    d_tolerance: float = 1e-10,
) -> dict[str, Any]:
    """Flip low-reliability logical variables and rerank unique cosets on D.

    Let ``x`` be the final H-BP hard decision, ``U`` the low-reliability
    variables, and

        V = {i in U : L[:, i] != 0}.

    The first stage directly enumerates every flip pattern supported on ``V``:

        delta_V in {0, 1}^|V|,
        logical_label(delta_V) = L (x xor delta).

    Duplicate labels are discarded before the unchanged D-BP stage.  This
    screening deliberately does not require ``H (x xor delta) = syndrome``;
    the original syndrome is used only when constructing each coset's canonical
    representative for D-BP.

    Exactly one of ``reliability_threshold`` and ``max_unfrozen_variables``
    must be supplied.  Reliability is the absolute value of the final belief.
    """
    syndrome = _gf2(syndrome).reshape(-1)
    p = np.asarray(error_probabilities, dtype=float).reshape(-1)
    if syndrome.size != model.H.shape[0]:
        raise ValueError("syndrome must have length H.shape[0].")
    if p.size != model.num_error_mechanisms:
        raise ValueError("error_probabilities must have length H.shape[1].")
    if isinstance(d_max_iter, (bool, np.bool_)) or not isinstance(
        d_max_iter,
        (int, np.integer),
    ):
        raise TypeError("d_max_iter must be an integer.")
    d_max_iter = int(d_max_iter)
    if d_max_iter < 1:
        raise ValueError("d_max_iter must be at least 1.")
    if (reliability_threshold is None) == (max_unfrozen_variables is None):
        raise ValueError(
            "Specify exactly one of reliability_threshold and "
            "max_unfrozen_variables."
        )

    if reliability_threshold is not None:
        if isinstance(reliability_threshold, (bool, np.bool_)) or not isinstance(
            reliability_threshold,
            (int, float, np.integer, np.floating),
        ):
            raise TypeError("reliability_threshold must be a real scalar.")
        reliability_threshold = float(reliability_threshold)
        if not np.isfinite(reliability_threshold) or reliability_threshold < 0.0:
            raise ValueError(
                "reliability_threshold must be finite and non-negative."
            )
    else:
        if isinstance(max_unfrozen_variables, (bool, np.bool_)) or not isinstance(
            max_unfrozen_variables,
            (int, np.integer),
        ):
            raise TypeError("max_unfrozen_variables must be an integer.")
        max_unfrozen_variables = int(max_unfrozen_variables)
        if not 0 <= max_unfrozen_variables <= model.num_error_mechanisms:
            raise ValueError(
                "max_unfrozen_variables must lie between 0 and H.shape[1]."
            )

    h_bp = parity_check_bp_trace(
        model.H,
        syndrome,
        p,
        max_iter=h_max_iter,
        damping=h_damping,
        tolerance=h_tolerance,
        method=h_method,
        min_sum_scaling=h_min_sum_scaling,
        llr_clip=h_llr_clip,
        stop_if_syndrome_satisfied=h_stop_if_syndrome_satisfied,
    )
    final_hard = h_bp.hard_decision_history[-1].copy()
    final_belief = h_bp.posterior_llrs.copy()
    reliability = np.abs(final_belief)
    if not np.all(np.isfinite(reliability)):
        raise ValueError("The final H-BP reliability contains a non-finite value.")

    variable_indices = np.arange(model.num_error_mechanisms, dtype=int)
    reliability_order = np.lexsort((variable_indices, reliability))
    if reliability_threshold is not None:
        unfrozen_variables = reliability_order[
            reliability[reliability_order] <= reliability_threshold
        ]
    else:
        unfrozen_variables = reliability_order[:max_unfrozen_variables]
    frozen_variables = np.setdiff1d(
        variable_indices,
        unfrozen_variables,
        assume_unique=True,
    )
    logical_support = np.any(model.L != 0, axis=0)
    logical_unfrozen_variables = unfrozen_variables[
        logical_support[unfrozen_variables]
    ]

    # Retained for diagnostics only; it does not filter the direct flip search.
    residual = syndrome ^ ((model.H @ final_hard) & 1)

    x_s = particular_solution(model, syndrome)
    candidates: list[dict[str, Any]] = []
    hypotheses: list[dict[str, Any]] = []
    observed_labels: set[tuple[int, ...]] = set()
    num_flip_variables = logical_unfrozen_variables.size
    for bits in product((0, 1), repeat=num_flip_variables):
        flip_pattern = np.fromiter(
            bits,
            dtype=np.uint8,
            count=num_flip_variables,
        )
        flip = np.zeros(
            model.num_error_mechanisms,
            dtype=np.uint8,
        )
        flip[logical_unfrozen_variables] = flip_pattern
        flipped_hard_decision = final_hard ^ flip

        observable_label = (model.L @ flipped_hard_decision) & 1
        label_key = _label_key(observable_label)
        if label_key in observed_labels:
            continue
        observed_labels.add(label_key)
        logical_coordinate = _logical_coordinate_for_label(
            model,
            x_s,
            observable_label,
        )
        representative = x_s ^ ((model.Q @ logical_coordinate) & 1)
        assert np.array_equal((model.H @ representative) & 1, syndrome)
        assert np.array_equal(
            (model.L @ representative) & 1,
            observable_label,
        )

        candidate = {
            "observable_label": observable_label.copy(),
            "logical_coordinate": logical_coordinate.copy(),
            "representative": representative.copy(),
            "flip_pattern": flip_pattern.copy(),
            "flip": flip.copy(),
            "flipped_hard_decision": flipped_hard_decision.copy(),
        }
        candidates.append(candidate)

        bp_result = quotient_bp_log_partition(
            model.D,
            representative,
            p,
            max_iter=d_max_iter,
            damping=d_damping,
            tolerance=d_tolerance,
        )
        hypothesis = dict(candidate)
        hypothesis["log_partition"] = float(bp_result.log_partition)
        hypothesis["bp"] = bp_result
        hypotheses.append(hypothesis)

    log_values = np.asarray(
        [item["log_partition"] for item in hypotheses],
        dtype=float,
    )
    best_index = int(np.argmax(log_values))
    best = hypotheses[best_index]

    finite = np.isfinite(log_values)
    shortlist_posterior = np.zeros_like(log_values)
    if finite.any():
        shifted = log_values[finite] - np.max(log_values[finite])
        weights = np.exp(shifted)
        shortlist_posterior[finite] = weights / weights.sum()
    for item, probability in zip(hypotheses, shortlist_posterior):
        item["shortlist_posterior"] = float(probability)

    correction = best["representative"].copy()
    assert np.array_equal((model.H @ correction) & 1, syndrome)
    assert np.array_equal(
        (model.L @ correction) & 1,
        best["observable_label"],
    )

    return {
        "correction": correction,
        "chosen_observable_label": best["observable_label"].copy(),
        "chosen_logical_coordinate": best["logical_coordinate"].copy(),
        "chosen_log_partition": float(best["log_partition"]),
        "particular_solution": x_s,
        "h_bp": h_bp,
        "final_hard_decision": final_hard,
        "final_belief": final_belief,
        "final_reliability": reliability,
        "residual_syndrome": residual,
        "frozen_variables": frozen_variables,
        "unfrozen_variables": unfrozen_variables,
        "logical_unfrozen_variables": logical_unfrozen_variables,
        "num_flip_patterns": 1 << num_flip_variables,
        "candidate_cosets": candidates,
        "candidate_count": len(candidates),
        "hypotheses": hypotheses,
        "screening": {
            "method": "reliability_frozen_logical_support_flips",
            "num_h_bp_iterations": h_bp.iterations,
            "num_frozen_variables": frozen_variables.size,
            "num_unfrozen_variables": unfrozen_variables.size,
            "num_logical_unfrozen_variables": logical_unfrozen_variables.size,
            "num_flip_patterns": 1 << num_flip_variables,
            "num_unique_candidate_cosets": len(candidates),
            "used_all_logical_cosets": len(candidates) == 2 ** model.logical_rank,
            "posterior_is_shortlist_normalized": True,
        },
    }


# Explicit descriptive alias for the main two-stage decoder.
decode_trajectory_screened_quotient_bp = decode_logical_coset_bp
