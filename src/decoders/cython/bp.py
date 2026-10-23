from typing import Literal
from dataclasses import dataclass

import numpy as np
import stim

from ..base import Decoder
from .cython_core.bp_core import min_sum_trace, product_sum_trace, research_trace


BP_TRACE_VERSION = "cython_research_bp_v1"


@dataclass(frozen=True)
class BPTrace:
    """Row-major edge histories, including the initial state at index zero."""

    c2v: np.ndarray
    v2c: np.ndarray
    llr: np.ndarray
    hard: np.ndarray
    valid: np.ndarray
    runtime_seconds: np.ndarray
    first_valid: int | None


class BP(Decoder):
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
        self.num_observables = None
        self.priors = None
        self.belief_history = np.empty((0, 0), dtype=np.float64)
        self._graph_ready = False

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None) -> None:
        self._parse(self.code, dem)

        if self.priors is None:
            p = float(self.error_rate if self.error_rate is not None else 1e-3)
            p = min(max(p, 1e-12), 1.0 - 1e-12)
            self.priors = np.asarray([p])
        self.ms_scaling_factor = self.config.get("ms_scaling_factor", 1.0)
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
        self._trace_graph = None

    @classmethod
    def from_pcm(cls, pcm, priors, **kwargs):
        """Prepare a decoder for an explicit PCM and independent error priors."""
        decoder = cls(**kwargs)
        decoder.H = np.array(pcm, dtype=np.uint8, order="C", copy=True) % 2
        assert decoder.H.ndim == 2
        decoder.m, decoder.n = decoder.H.shape
        decoder.priors = np.array(priors, dtype=np.float64, copy=True)
        assert decoder.priors.shape == (decoder.n,)
        assert np.all(np.isfinite(decoder.priors))
        assert np.all((decoder.priors > 0) & (decoder.priors < 1))
        decoder.check_conn = [np.flatnonzero(row) for row in decoder.H]
        decoder.ms_scaling_factor = decoder.config.get("ms_scaling_factor", 1.0)
        decoder._prepare_cython_graph()
        return decoder

    def trace(
        self,
        syndrome,
        iterations=None,
        *,
        update_rule="min_sum",
        scaling=1.0,
        damping=0.0,
        offset=0.0,
        schedule="parallel",
        rng=None,
        prior_llr=None,
        initial_message_sigma=0.0,
        stop_early=False,
        record_messages=True,
        include_initial_valid=True,
        accumulation="sequential",
        serial_accumulation="sequential",
        combine_scaling_damping=False,
    ) -> BPTrace:
        """Inspect BP without changing normal decode/OSD stopping semantics.

        Serial schedules update checks, unlike ldpc's variable serial schedule.
        Accumulation modes preserve the corresponding experiment's reduction
        order: sequential prior-plus-messages, sum-then-prior, or NumPy sum.
        Runtime includes per-shot native initialization and history allocation;
        it excludes graph preparation and Python argument preparation.
        """
        if not self._graph_ready:
            self._init_backend()
        syndrome = np.ascontiguousarray(np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2)
        iterations = self.max_iter if iterations is None else iterations
        assert isinstance(iterations, (int, np.integer)) and iterations >= 1
        assert syndrome.shape == (self.m,)
        assert update_rule in {"min_sum", "sum_product"}
        assert schedule in {"parallel", "serial", "random_serial_fixed", "random_serial"}
        assert np.isfinite(offset) and offset >= 0
        assert not combine_scaling_damping or offset == 0
        assert np.isfinite(initial_message_sigma) and initial_message_sigma >= 0
        modes = {"sequential": 0, "sum_then_prior": 1, "numpy": 2}
        assert accumulation in modes and serial_accumulation in modes

        if self._trace_graph is None:
            checks, variables = np.nonzero(self.H)
            check_ptr = np.concatenate(([0], np.cumsum(np.bincount(checks, minlength=self.m)))).astype(np.int32)
            variable_ptr = np.concatenate(([0], np.cumsum(np.bincount(variables, minlength=self.n)))).astype(np.int32)
            variable_edges = np.argsort(variables, kind="stable").astype(np.int32)
            self._trace_graph = (check_ptr, variables.astype(np.int32), variable_ptr, variable_edges)
        check_ptr, variables, variable_ptr, variable_edges = self._trace_graph
        edges = variables.size
        prior = np.ascontiguousarray(self._channel_llr if prior_llr is None else prior_llr, dtype=np.float64)
        if prior.size == 1:
            prior = np.full(self.n, prior[0], dtype=np.float64)
        assert prior.shape == (self.n,) and np.all(np.isfinite(prior))
        scaling = np.array(np.broadcast_to(scaling, (edges,)), dtype=np.float64, order="C", copy=True)
        damping = np.array(np.broadcast_to(damping, (edges,)), dtype=np.float64, order="C", copy=True)
        assert np.all(np.isfinite(scaling)) and np.all(scaling > 0)
        assert np.all(np.isfinite(damping)) and np.all((damping >= 0) & (damping < 1))

        rng = np.random.default_rng(0) if rng is None else rng
        initial_v2c = np.ascontiguousarray(prior[variables])
        if initial_message_sigma:
            initial_v2c += rng.normal(0, initial_message_sigma, edges)
        order = np.arange(self.m, dtype=np.int32)
        if schedule == "random_serial_fixed":
            order = np.asarray(rng.permutation(self.m), dtype=np.int32)
        if schedule == "random_serial":
            orders = np.empty((iterations, self.m), dtype=np.int32)
            for iteration in range(iterations):
                # The reference reshuffles a fresh canonical order each round.
                order = np.arange(self.m, dtype=np.int32)
                rng.shuffle(order)
                orders[iteration] = order
        else:
            orders = order.reshape(1, -1)
        result = research_trace(
            syndrome, prior, check_ptr, variables, variable_ptr, variable_edges,
            iterations, scaling, damping, offset, update_rule == "sum_product",
            orders, schedule != "parallel", initial_v2c, modes[accumulation],
            modes[serial_accumulation], include_initial_valid, stop_early,
            record_messages, combine_scaling_damping,
        )
        return BPTrace(*result)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
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
            )
        else:
            raise ValueError(f"Unsupported BP method: {self.bp_method}")
        return (llr < 0).astype(np.uint8)
