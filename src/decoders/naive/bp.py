from typing import Literal

import numpy as np
import stim

from ..base import Decoder


class BP(Decoder):
    """Plain belief-propagation decoder on the static CSS Tanner graph."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.max_iter = max_iter
        self.bp_method = bp_method

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None, need_parse: bool = True) -> None:
        if need_parse:
            self._parse(self.code, dem)

        p = float(self.error_rate if self.error_rate is not None else 1e-3)
        p = min(max(p, 1e-12), 1.0 - 1e-12)
        self.prior = float(np.log((1.0 - p) / p))

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if not hasattr(self, "H"):
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).flatten() % 2
        llr = self.min_sum(syndrome) if self.bp_method == "minimum_sum" else self.prod_sum(syndrome)
        return (llr < 0).astype(np.uint8)

    def min_sum(self, syndrome: np.ndarray) -> np.ndarray:
        msg_v2c = np.zeros((self.n, self.m), dtype=float)
        msg_c2v = np.zeros((self.m, self.n), dtype=float)

        for j in range(self.n):
            for i in self.var_conn[j]:
                msg_v2c[j, i] = self.prior

        for _ in range(self.max_iter):
            for i in range(self.m):
                syndrome_sign = -1.0 if int(syndrome[i]) else 1.0
                for j in self.check_conn[i]:
                    incoming = [msg_v2c[k, i] for k in self.check_conn[i] if k != j]
                    if not incoming:
                        msg_c2v[i, j] = 0.0
                        continue

                    msg_c2v[i, j] = syndrome_sign * np.prod(np.sign(incoming)) * np.min(np.abs(incoming))

            for j in range(self.n):
                for i in self.var_conn[j]:
                    incoming = [msg_c2v[k, j] for k in self.var_conn[j] if k != i]
                    msg_v2c[j, i] = self.prior + float(sum(incoming))

            llr = np.zeros(self.n, dtype=float)
            for j in range(self.n):
                llr[j] = self.prior + float(sum(msg_c2v[i, j] for i in self.var_conn[j]))
            correction = (llr < 0).astype(np.uint8)
            if np.array_equal((self.H @ correction) % 2, syndrome):
                return llr

        llr = np.zeros(self.n, dtype=float)
        for j in range(self.n):
            llr[j] = self.prior + float(sum(msg_c2v[i, j] for i in self.var_conn[j]))
        return llr

    def prod_sum(self, syndrome: np.ndarray) -> np.ndarray:
        msg_v2c = np.zeros((self.n, self.m), dtype=float)
        msg_c2v = np.zeros((self.m, self.n), dtype=float)

        for j in range(self.n):
            for i in self.var_conn[j]:
                msg_v2c[j, i] = self.prior

        for _ in range(self.max_iter):
            for i in range(self.m):
                syndrome_sign = -1.0 if int(syndrome[i]) else 1.0
                for j in self.check_conn[i]:
                    incoming = [msg_v2c[k, i] for k in self.check_conn[i] if k != j]
                    if not incoming:
                        msg_c2v[i, j] = 0.0
                        continue

                    tanh_prod = np.prod(np.tanh(np.array(incoming, dtype=float) / 2.0))
                    tanh_prod = float(np.clip(tanh_prod, -1.0 + 1e-12, 1.0 - 1e-12))
                    msg_c2v[i, j] = 2.0 * syndrome_sign * np.arctanh(tanh_prod)

            for j in range(self.n):
                for i in self.var_conn[j]:
                    incoming = [msg_c2v[k, j] for k in self.var_conn[j] if k != i]
                    msg_v2c[j, i] = self.prior + float(sum(incoming))

            llr = np.zeros(self.n, dtype=float)
            for j in range(self.n):
                llr[j] = self.prior + float(sum(msg_c2v[i, j] for i in self.var_conn[j]))
            correction = (llr < 0).astype(np.uint8)
            if np.array_equal((self.H @ correction) % 2, syndrome):
                return llr

        llr = np.zeros(self.n, dtype=float)
        for j in range(self.n):
            llr[j] = self.prior + float(sum(msg_c2v[i, j] for i in self.var_conn[j]))
        return llr
