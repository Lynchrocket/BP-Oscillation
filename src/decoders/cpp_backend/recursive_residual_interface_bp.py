from typing import Literal

import numpy as np
import stim
from ldpc_modified import BpDecoder, RecursiveResidualInterfaceRefiner

from .recursive_frontier_bp import RecursiveFrontierBP, _dense_binary


class RecursiveResidualInterfaceBP(RecursiveFrontierBP):
    """BP followed by recursive region-interface inference and residual compression."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["product_sum", "minimum_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        schedule: Literal["parallel", "serial", "serial_relative"] = "parallel",
        omp_thread_count: int = 1,
        max_interface_states: int = 512,
        hierarchy_mode: Literal["recursive", "generic", "shuffled"] = "recursive",
        compression_mode: Literal["residual", "mass"] = "residual",
        hierarchy_seed: int = 20260903,
        **kwargs,
    ):
        assert not kwargs, f"Unknown RecursiveResidualInterfaceBP parameters: {sorted(kwargs)}"
        assert compression_mode in {"residual", "mass"}
        super().__init__(
            code=code,
            error_rate=error_rate,
            max_iter=max_iter,
            bp_method=bp_method,
            ms_scaling_factor=ms_scaling_factor,
            schedule=schedule,
            omp_thread_count=omp_thread_count,
            max_frontier_states=max_interface_states,
            ordering_mode=hierarchy_mode,
            ordering_seed=hierarchy_seed,
        )
        self.max_interface_states = int(max_interface_states)
        self.hierarchy_mode = hierarchy_mode
        self.compression_mode = compression_mode
        self.hierarchy_seed = int(hierarchy_seed)

    def update_params(self, **params):
        allowed = {
            "error_rate", "max_iter", "bp_method", "ms_scaling_factor", "schedule",
            "omp_thread_count", "max_interface_states", "hierarchy_mode",
            "compression_mode", "hierarchy_seed",
        }
        unknown = set(params) - allowed
        assert not unknown, f"Unknown RecursiveResidualInterfaceBP parameters: {sorted(unknown)}"
        if "max_interface_states" in params:
            assert isinstance(params["max_interface_states"], (int, np.integer))
            assert 0 < params["max_interface_states"] <= np.iinfo(np.int32).max
        if "hierarchy_mode" in params:
            assert params["hierarchy_mode"] in {"recursive", "generic", "shuffled"}
        if "compression_mode" in params:
            assert params["compression_mode"] in {"residual", "mass"}
        if "schedule" in params:
            assert params["schedule"] in {"parallel", "serial", "serial_relative"}
        if "omp_thread_count" in params:
            assert isinstance(params["omp_thread_count"], (int, np.integer))
            assert params["omp_thread_count"] > 0
        if "hierarchy_seed" in params:
            assert isinstance(params["hierarchy_seed"], (int, np.integer))
        for key, value in params.items():
            setattr(self, key, value)
        self.max_frontier_states = self.max_interface_states
        self.ordering_mode = self.hierarchy_mode
        self.ordering_seed = self.hierarchy_seed
        if self.bp is not None:
            self._init_backend(self._dem)

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._dem = dem
        self.max_frontier_states = self.max_interface_states
        self.ordering_mode = self.hierarchy_mode
        self.ordering_seed = self.hierarchy_seed
        if dem is None and self._pcm_override is not None:
            self.mode = "static"
            self.H = self._pcm_override.copy()
            self.m, self.n = self.H.shape
            self.var_conn = [np.flatnonzero(self.H[:, j]) for j in range(self.n)]
            self.check_conn = [np.flatnonzero(self.H[i, :]) for i in range(self.m)]
        else:
            self._parse(self.code, dem)
        priors = (
            np.asarray(self.priors, dtype=np.float64)
            if self.mode == "dem"
            else np.full(self.n, self.error_rate, dtype=np.float64)
        )
        self.bp = BpDecoder(
            pcm=self.H,
            error_channel=priors,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            schedule=self.schedule,
            omp_thread_count=self.omp_thread_count,
            input_vector_type="syndrome",
        )
        with np.errstate(divide="ignore"):
            self._initial_llr = np.log1p(-priors) - np.log(priors)
        self.variable_order, _ = self._variable_order(dem)
        self.bag_ends = None
        self.refiner = RecursiveResidualInterfaceRefiner(
            _dense_binary(self.H),
            self._logical_matrix(),
            priors,
            self.variable_order,
            self.max_interface_states,
            self.compression_mode == "residual",
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()
        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        assert syndrome.size == self.m
        main = np.asarray(self.bp.decode(syndrome), dtype=np.uint8)
        self.last_main_correction = main.copy()
        correction = self.refiner.refine(syndrome, main, self.bp.log_prob_ratios)
        self.last_stats = dict(self.refiner.stats)
        self.last_stats["main_bp_converged"] = bool(self.bp.converge)
        return correction
