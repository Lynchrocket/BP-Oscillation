from itertools import combinations
from typing import Literal

import numpy as np
import stim

from ..base import Decoder
from .bp import BP


class BPOSD(Decoder):
    """BP decoder with a lightweight ordered-statistics post-processing step."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        osd_order: int = 2,
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.max_iter = max_iter
        self.bp_method = bp_method
        self.osd_order = osd_order
        self.bp = None

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)
        self.bp = BP(self.code, self.error_rate, self.max_iter, self.bp_method)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).flatten() % 2

        llr = self.bp.min_sum(syndrome) if self.bp_method == "minimum_sum" else self.bp.prod_sum(syndrome)
        correction = (llr < 0).astype(np.uint8)

        if np.array_equal((self.H @ correction) % 2, syndrome):
            return correction

        return self._ordered_statistics_postprocess(llr, syndrome)

    def _ordered_statistics_postprocess(self, llr: np.ndarray, syndrome: np.ndarray) -> np.ndarray:
        m, n = self.H.shape

        reliability_order = np.argsort(np.abs(llr))
        H_ordered = self.H[:, reliability_order]

        # OSD-0
        J = self._select_info_set(H_ordered)
        J_bar = np.setdiff1d(np.arange(n), J)
        
        Lambda_J = H_ordered[:, J]
        Lambda_J_bar = H_ordered[:, J_bar]

        Lambda_J_inv = self.gf2_inv(Lambda_J)
        e0 = (Lambda_J_inv @ syndrome) % 2

        best_e = np.zeros(n)
        best_e[J] = e0
        best_weight = np.sum(e0)

        # OSD-w
        for w in range(1, self.osd_order+1):
            for flip_indices in combinations(range(len(J_bar)), w):
                e_J_bar = np.zeros(len(J_bar))
                e_J_bar[list(flip_indices)] = 1

                e_J = (e0 + (Lambda_J_inv @ (Lambda_J_bar @ e_J_bar % 2) % 2)) % 2
                e = np.zeros(n)
                e[J] = e_J
                e[J_bar] = e_J_bar

                weight = np.sum(e)
                if weight < best_weight:
                    best_e = e
                    best_weight = weight

        return best_e

    def _select_info_set(self, H):
        """Select the information set using Gaussian elimination."""
        H = H.copy() % 2
        m, n = H.shape

        pivot_cols = []
        row = 0

        for col in range(n):
            if np.all(H[row:, col] == 0):
                continue
            pivot = np.argmax(H[row:, col] == 1)

            # swap
            if pivot != row:
                H[[row, pivot]] = H[[pivot, row]]

            for r in range(m):
                if r != row and H[r, col] == 1:
                    H[r] ^= H[row]

            pivot_cols.append(col)
            row += 1

            if row == m:
                break

        return np.array(pivot_cols)
    
    def gf2_inv(self, A):
        A = A.copy()
        n = A.shape[0]
        I = np.eye(n, dtype=int)

        for i in range(n):
            if A[i, i] == 0:
                for j in range(i+1, n):
                    if A[j, i] == 1:
                        A[[i, j]] = A[[j, i]]
                        I[[i, j]] = I[[j, i]]
                        break

            for j in range(n):
                if j != i and A[j, i] == 1:
                    A[j] ^= A[i]
                    I[j] ^= I[i]

        return I
    