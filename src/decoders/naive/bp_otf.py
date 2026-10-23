"""
BP+OTF Decoder - Belief Propagation with Ordered Tanner forest

An almost-linear time decoding algorithm for quantum LDPC codes under circuit-level noise (Antonio et al.)
"""

from typing import Literal

import numpy as np
import stim

from ..base import Decoder
from .bp import BP


class BPOTF(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 0,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)

        self.max_iter = max_iter
        self.bp_method = bp_method
        self.bp = None

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse_pcm(self.code)
        self.bp = BP(self.code, self.error_rate, self.max_iter, self.bp_method)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        llr = self.bp.min_sum(syndrome) if self.bp_method == "minimum_sum" else self.bp.prod_sum(syndrome)
        correction = (llr < 0).astype(np.uint8)
        if np.array_equal((self.H @ correction) % 2, syndrome):
            return correction

        h_otf = self._build_otf(self.H, llr)
        self.bp.update_params(H=h_otf)
        self.bp._init_backend(need_parse=False)
        llr_otf = self.bp.prod_sum(syndrome)
        return (llr_otf < 0).astype(np.uint8)

    def _build_otf(self, h: np.ndarray, llr: np.ndarray) -> np.ndarray:
        m, n = h.shape
        order = np.argsort(-np.abs(llr))
        uf = UnionFind(n + m)

        h_otf = np.zeros_like(h)
        for v in order:
            checks = np.where(h[:, v])[0]
            for c in checks:
                if uf.union(int(v), int(c + n)):
                    h_otf[c, v] = 1

        return h_otf


class UnionFind:
    def __init__(self, n: int):
        self.p = list(range(n))

    def find(self, x: int) -> int:
        if self.p[x] != x:
            self.p[x] = self.find(self.p[x])
        return self.p[x]

    def union(self, a: int, b: int) -> bool:
        pa, pb = self.find(a), self.find(b)
        if pa == pb:
            return False
        self.p[pb] = pa
        return True
