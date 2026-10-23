from typing import Literal

import numpy as np
import stim

from ..base import Decoder
from .cython_core.bp_core import min_sum_trace, product_sum_trace
from .cython_core.osd_core import ordered_statistics_postprocess


class BPOSD(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 60,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        osd_order: int = 10,
        osd_method: str = "osd_cs",
        need_tracking: bool = True,
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.max_iter = max_iter
        self.bp_method = bp_method
        self.osd_order = osd_order
        assert osd_method == "osd_cs", "Only the benchmark osd_cs search is supported"
        self.osd_method = osd_method
        self.ms_scaling_factor = ms_scaling_factor
        self.num_observables = None
        self.priors = None
        self.need_tracking = need_tracking
        
        self.belief_history = None
        self._graph_ready = False
        
        self.is_bp_converged = False
        self.is_osd_satisfied = False

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)

        if self.priors is None:
            p = float(self.error_rate if self.error_rate is not None else 1e-3)
            p = min(max(p, 1e-12), 1.0 - 1e-12)
            self.priors = np.full(self.n, p, dtype=np.float64)
        self._prepare_cython_graph()

    def _prepare_cython_graph(self):
        self.H = np.ascontiguousarray(self.H, dtype=np.uint8) % 2

        # edges
        num_edges = np.sum(self.H)
        msg = np.zeros(2 * num_edges, dtype=np.float64) # c2v and v2c are stored interchangebly
        head = np.full(self.m + self.n, -1, dtype=np.int32) # head[u] -> 1st edge starting from node u
        head_nxt = np.full(2 * num_edges, -1, dtype=np.int32) # head_nxt[edge] -> next edge from the same head node
        to = np.full(2 * num_edges, -1, dtype=np.int32) # to[edge] -> destination node of the edge
        def _add(u, v, cnt):
            to[cnt] = v
            head_nxt[cnt] = head[u]
            head[u] = cnt
        cnt = 0
        for i, conn in enumerate(self.check_conn):
            for j in conn:
                _add(i, j + self.m, cnt) # check to variable edge
                cnt += 1
                _add(j + self.m, i, cnt) # variable to check edge
                cnt += 1
        
        self._channel_llr = np.ascontiguousarray(np.log((1.0 - self.priors) / self.priors))
        self._msg = np.ascontiguousarray(msg)
        self._head = np.ascontiguousarray(head)
        self._head_nxt = np.ascontiguousarray(head_nxt)
        self._to = np.ascontiguousarray(to)
        self._graph_ready = True

    def decode_to_correction(self, syndrome: np.ndarray, force_osd: bool = False) -> np.ndarray:
        if not self._graph_ready:
            self._init_backend()

        syndrome = np.ascontiguousarray(np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2)
        if syndrome.size != self.m:
            raise ValueError(f"syndrome has length {syndrome.size}, expected {self.m}")

        if self.bp_method == "minimum_sum":
            llr, self.belief_history = min_sum_trace(
                syndrome,
                self._msg,
                self._head,
                self._head_nxt,
                self._to,
                self.m,
                self.n,
                self._channel_llr,
                self.max_iter,
                self.ms_scaling_factor,
                self.need_tracking,
            )
        elif self.bp_method == "product_sum":
            llr, self.belief_history = product_sum_trace(
                syndrome,
                self._msg,
                self._head,
                self._head_nxt,
                self._to,
                self.m,
                self.n,
                self._channel_llr,
                self.max_iter,
                self.need_tracking,
            )
        else:
            raise ValueError(f"unsupported bp_method: {self.bp_method!r}")

        self.bp_correction = (llr < 0).astype(np.uint8)
        if np.array_equal((self.H @ self.bp_correction) % 2, syndrome):
            self.is_bp_converged = True
            if not force_osd:
                return self.bp_correction
        
        self.osd_correction = ordered_statistics_postprocess(
            np.ascontiguousarray(llr), syndrome, self.H, self.osd_order,
            np.ascontiguousarray(self.priors, dtype=np.float64), self.osd_method,
        )
        if np.array_equal((self.H @ self.osd_correction) % 2, syndrome):
            self.is_osd_satisfied = True

        return self.osd_correction
