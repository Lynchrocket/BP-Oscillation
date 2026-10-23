from typing import Literal

import numpy as np
import stim
from ldpc_modified import BpDecoder, CorrelationAdaptiveRegionRefiner

from .recursive_loop_core_bp import RecursiveLoopCoreBP, _dense_binary


_INT32_MIN = -(2**31)
_INT32_MAX = 2**31 - 1


def _is_integer(value) -> bool:
    return isinstance(value, (int, np.integer)) and not isinstance(
        value, (bool, np.bool_)
    )


def _is_finite_number(value) -> bool:
    return (
        isinstance(value, (int, float, np.integer, np.floating))
        and not isinstance(value, (bool, np.bool_))
        and bool(np.isfinite(value))
    )


class CorrelationAdaptiveBP(RecursiveLoopCoreBP):
    """BP followed by factorization-gap-guided exact regional inference."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["product_sum", "minimum_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        schedule: Literal["parallel", "serial", "serial_relative"] = "parallel",
        omp_thread_count: int = 1,
        max_region_states: int = 4096,
        max_probe_regions: int = 32,
        region_mode: Literal["generic", "recursive", "shuffled"] = "recursive",
        score_mode: Literal[
            "full_gap", "correlation_only", "marginal_only", "structural", "random"
        ] = "full_gap",
        region_seed: int = 20260903,
        refine_zero_syndrome: bool = False,
        **kwargs,
    ):
        assert not kwargs, f"Unknown CorrelationAdaptiveBP parameters: {sorted(kwargs)}"
        assert _is_finite_number(error_rate) and 0 < error_rate < 1
        assert _is_integer(max_iter) and max_iter > 0
        assert bp_method in {"product_sum", "minimum_sum"}
        assert _is_finite_number(ms_scaling_factor) and ms_scaling_factor > 0
        assert schedule in {"parallel", "serial", "serial_relative"}
        assert _is_integer(omp_thread_count) and omp_thread_count > 0
        assert _is_integer(max_region_states)
        assert 2 <= max_region_states <= 2**62
        assert _is_integer(max_probe_regions)
        assert 0 < max_probe_regions <= _INT32_MAX
        assert region_mode in {"generic", "recursive", "shuffled"}
        assert score_mode in {
            "full_gap", "correlation_only", "marginal_only", "structural", "random"
        }
        assert _is_integer(region_seed)
        assert _INT32_MIN <= region_seed <= _INT32_MAX
        assert isinstance(refine_zero_syndrome, (bool, np.bool_))
        super().__init__(
            code=code,
            error_rate=error_rate,
            max_iter=max_iter,
            bp_method=bp_method,
            ms_scaling_factor=ms_scaling_factor,
            schedule=schedule,
            omp_thread_count=omp_thread_count,
            max_factor_states=max_region_states,
            cluster_mode="generic",
            sum_logical_mass=True,
        )
        self.max_region_states = int(max_region_states)
        self.max_probe_regions = int(max_probe_regions)
        self.region_mode = region_mode
        self.score_mode = score_mode
        self.region_seed = int(region_seed)
        self.refine_zero_syndrome = bool(refine_zero_syndrome)

    def update_params(self, **params):
        allowed = {
            "error_rate", "max_iter", "bp_method", "ms_scaling_factor", "schedule",
            "omp_thread_count", "max_region_states", "max_probe_regions", "region_mode",
            "score_mode", "region_seed", "refine_zero_syndrome",
        }
        unknown = set(params) - allowed
        assert not unknown, f"Unknown CorrelationAdaptiveBP parameters: {sorted(unknown)}"
        if "error_rate" in params:
            assert _is_finite_number(params["error_rate"])
            assert 0 < params["error_rate"] < 1
        if "max_iter" in params:
            assert _is_integer(params["max_iter"]) and params["max_iter"] > 0
        if "bp_method" in params:
            assert params["bp_method"] in {"product_sum", "minimum_sum"}
        if "ms_scaling_factor" in params:
            assert _is_finite_number(params["ms_scaling_factor"])
            assert params["ms_scaling_factor"] > 0
        if "max_region_states" in params:
            assert _is_integer(params["max_region_states"])
            assert 2 <= params["max_region_states"] <= 2**62
        if "max_probe_regions" in params:
            assert _is_integer(params["max_probe_regions"])
            assert 0 < params["max_probe_regions"] <= _INT32_MAX
        if "region_mode" in params:
            assert params["region_mode"] in {"generic", "recursive", "shuffled"}
        if "score_mode" in params:
            assert params["score_mode"] in {
                "full_gap", "correlation_only", "marginal_only", "structural", "random"
            }
        if "region_seed" in params:
            assert _is_integer(params["region_seed"])
            assert _INT32_MIN <= params["region_seed"] <= _INT32_MAX
        if "refine_zero_syndrome" in params:
            assert isinstance(params["refine_zero_syndrome"], (bool, np.bool_))
        if "schedule" in params:
            assert params["schedule"] in {"parallel", "serial", "serial_relative"}
        if "omp_thread_count" in params:
            assert _is_integer(params["omp_thread_count"])
            assert params["omp_thread_count"] > 0
        for key, value in params.items():
            setattr(self, key, value)
        if self.bp is not None:
            self._init_backend(self._dem)

    def _logical_matrix(self) -> np.ndarray:
        if self.mode == "dem":
            return _dense_binary(self.L)
        n = int(self.code.n)
        lz = _dense_binary(self.code.L_Z)
        if self.H.shape[1] == n:
            assert lz.shape[1] == n
            return lz
        lx = _dense_binary(self.code.L_X)
        assert self.H.shape[1] == 2 * n
        return np.block([
            [lz, np.zeros((lz.shape[0], n), dtype=np.uint8)],
            [np.zeros((lx.shape[0], n), dtype=np.uint8), lx],
        ])

    def _logical_from_correction(self, correction: np.ndarray) -> np.ndarray:
        if self.mode == "dem" or self.H.shape[1] == 2 * self.code.n:
            return super()._logical_from_correction(correction)
        correction = np.asarray(correction, dtype=np.uint8).reshape(-1) % 2
        assert correction.size == self.code.n
        return (_dense_binary(self.code.L_Z) @ correction) % 2

    def _static_recursive_metadata(self):
        times = np.zeros(self.m, dtype=np.int32)
        groups = np.arange(self.m, dtype=np.int32)
        layers = np.full(self.m, -1, dtype=np.int32)
        if self.code is None or not hasattr(self.code, "rows"):
            return times, groups, layers
        variable_layer = np.full(self.n, -1, dtype=np.int32)
        base_n = int(self.code.n)
        assert self.n in {base_n, 2 * base_n}
        for layer, variables in enumerate(self.code.rows):
            for variable in variables:
                assert 0 <= variable < base_n
                assert variable_layer[variable] < 0
                variable_layer[variable] = layer
                if self.n == 2 * base_n:
                    variable_layer[variable + base_n] = layer
        assert np.all(variable_layer >= 0), "code.rows must partition all base variables"
        for check in range(self.m):
            support = np.flatnonzero(self.H[check])
            if support.size and np.all(variable_layer[support] >= 0):
                layers[check] = int(np.max(variable_layer[support]))
        return times, groups, layers

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._dem = dem
        self._parse(self.code, dem)
        if dem is None and self.region_mode == "recursive":
            assert self.code is not None and hasattr(self.code, "rows"), (
                "region_mode='recursive' requires code.rows metadata"
            )
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
        if dem is None:
            times, groups, layers = self._static_recursive_metadata()
        else:
            times, groups, layers = self._check_metadata(dem)
        region_mode = {"generic": 0, "recursive": 1, "shuffled": 2}[self.region_mode]
        score_mode = {
            "full_gap": 0,
            "correlation_only": 1,
            "marginal_only": 2,
            "structural": 3,
            "random": 4,
        }[self.score_mode]
        self.refiner = CorrelationAdaptiveRegionRefiner(
            _dense_binary(self.H),
            self._logical_matrix(),
            priors,
            times,
            groups,
            layers,
            self.max_region_states,
            self.max_probe_regions,
            region_mode,
            score_mode,
            self.region_seed,
            self.refine_zero_syndrome,
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()
        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        assert syndrome.size == self.m
        main = np.asarray(self.bp.decode(syndrome), dtype=np.uint8)
        self.last_main_correction = main.copy()
        posterior_llr = self.bp.log_prob_ratios if np.any(syndrome) else self._initial_llr
        correction = self.refiner.refine(syndrome, main, posterior_llr)
        assert np.array_equal((_dense_binary(self.H) @ correction) % 2, syndrome)
        self.last_stats = dict(self.refiner.stats)
        self.last_stats["main_bp_converged"] = bool(self.bp.converge)
        return correction

    def converge(self) -> bool:
        return bool(self.last_stats.get("syndrome_valid", False))
