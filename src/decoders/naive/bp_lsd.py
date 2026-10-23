from itertools import combinations
from typing import Literal

import numpy as np
import stim

from ..base import Decoder
from .bp import BP


class BPLSD(Decoder):
    """BP decoder with a lightweight localized-statistics post-processing step."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        lsd_order: int = 2,
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)

        self.max_iter = max_iter
        self.bp_method = bp_method
        self.lsd_order = lsd_order
        self.bp = None

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse_pcm(self.code)
        self.bp = BP(self.code, self.error_rate, self.max_iter, self.bp_method)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).flatten() % 2
        llr = self.bp.min_sum(syndrome) if self.bp_method == "minimum_sum" else self.bp.prod_sum(syndrome)
        correction = (llr < 0).astype(np.uint8)

        if np.array_equal((self.H @ correction) % 2, syndrome):
            return correction

        return self._localized_statistics_postprocess(llr, syndrome)

    def _localized_statistics_postprocess(self, llr: np.ndarray, syndrome: np.ndarray) -> np.ndarray:
        current = (llr < 0).astype(np.uint8)
        residual = ((self.H @ current) % 2) ^ syndrome

        if not np.any(residual) or self.lsd_order <= 0:
            return current

        clusters = self._unsatisfied_check_clusters(residual)
        for cluster_checks in clusters:
            local_delta = self._solve_local_cluster(cluster_checks, residual, llr)
            if local_delta is None:
                continue
            current ^= local_delta
            residual = ((self.H @ current) % 2) ^ syndrome
            if not np.any(residual):
                return current

        global_delta = self._solve_global_residual(residual, llr)
        if global_delta is not None:
            current ^= global_delta

        return current

    def _unsatisfied_check_clusters(self, residual: np.ndarray) -> list[np.ndarray]:
        unsatisfied = set(np.flatnonzero(residual))
        clusters: list[np.ndarray] = []

        while unsatisfied:
            seed = unsatisfied.pop()
            queue = [seed]
            cluster_checks = {seed}
            cluster_vars = set()

            while queue:
                check = queue.pop()
                for var in self.check_conn[check]:
                    if var in cluster_vars:
                        continue
                    cluster_vars.add(var)
                    for neighbor_check in self.var_conn[var]:
                        if neighbor_check in unsatisfied:
                            unsatisfied.remove(neighbor_check)
                            cluster_checks.add(neighbor_check)
                            queue.append(neighbor_check)

            clusters.append(np.array(sorted(cluster_checks), dtype=int))

        return clusters

    def _solve_local_cluster(
        self,
        cluster_checks: np.ndarray,
        residual: np.ndarray,
        llr: np.ndarray,
    ) -> np.ndarray | None:
        candidate_vars = sorted(
            {
                var
                for check in cluster_checks
                for var in self.check_conn[int(check)]
            }
        )
        if not candidate_vars:
            return None

        ranked_vars = self._rank_candidates(candidate_vars, cluster_checks, llr)
        max_candidates = min(
            len(ranked_vars),
            max(self.lsd_order, len(cluster_checks)) + 2 * self.lsd_order,
        )
        ranked_vars = ranked_vars[:max_candidates]

        local_target = np.zeros(self.m, dtype=np.uint8)
        local_target[cluster_checks] = residual[cluster_checks]

        return self._search_exact_delta(ranked_vars, llr, local_target)

    def _solve_global_residual(self, residual: np.ndarray, llr: np.ndarray) -> np.ndarray | None:
        unsatisfied_checks = np.flatnonzero(residual)
        if unsatisfied_checks.size == 0:
            return np.zeros(self.n, dtype=np.uint8)

        candidate_vars = sorted(
            {
                var
                for check in unsatisfied_checks
                for var in self.check_conn[int(check)]
            }
        )
        ranked_vars = self._rank_candidates(candidate_vars, unsatisfied_checks, llr)
        max_candidates = min(
            len(ranked_vars),
            max(self.lsd_order, unsatisfied_checks.size) + 3 * self.lsd_order,
        )
        ranked_vars = ranked_vars[:max_candidates]

        return self._search_exact_delta(ranked_vars, llr, residual)

    def _rank_candidates(
        self,
        candidate_vars: list[int],
        active_checks: np.ndarray,
        llr: np.ndarray,
    ) -> list[int]:
        active_set = set(map(int, active_checks))

        def key(var: int) -> tuple[float, int, int]:
            support = sum(1 for check in self.var_conn[var] if check in active_set)
            return (float(abs(llr[var])), -support, var)

        return sorted(candidate_vars, key=key)

    def _search_exact_delta(
        self,
        candidate_vars: list[int],
        llr: np.ndarray,
        target: np.ndarray,
    ) -> np.ndarray | None:
        if not np.any(target):
            return np.zeros(self.n, dtype=np.uint8)

        best_delta = None
        best_cost = None
        max_weight = min(self.lsd_order, len(candidate_vars))

        for weight in range(1, max_weight + 1):
            for subset in combinations(candidate_vars, weight):
                delta = np.zeros(self.n, dtype=np.uint8)
                delta[list(subset)] = 1

                if not np.array_equal((self.H @ delta) % 2, target):
                    continue

                cost = float(np.sum(np.abs(llr[list(subset)])))
                if best_cost is None or cost < best_cost:
                    best_delta = delta
                    best_cost = cost

            if best_delta is not None:
                break

        return best_delta
