from typing import Literal

import numpy as np
import scipy.sparse as sp
import stim
from ldpc_modified import BpDecoder, RecursiveLoopCoreRefiner

from ..base import Decoder


def _dense_binary(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    return np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) % 2)


class RecursiveLoopCoreBP(Decoder):
    """BP followed by structure-guided, bounded exact local contraction."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["product_sum", "minimum_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        schedule: Literal["parallel", "serial", "serial_relative"] = "parallel",
        omp_thread_count: int = 1,
        max_factor_states: int = 65_536,
        cluster_mode: Literal["tree", "generic", "recursive"] = "recursive",
        sum_logical_mass: bool = True,
        **kwargs,
    ):
        assert not kwargs, f"Unknown RecursiveLoopCoreBP parameters: {sorted(kwargs)}"
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        assert isinstance(max_factor_states, (int, np.integer)) and max_factor_states > 0
        assert cluster_mode in {"tree", "generic", "recursive"}
        assert schedule in {"parallel", "serial", "serial_relative"}
        assert isinstance(omp_thread_count, (int, np.integer)) and omp_thread_count > 0
        self.max_iter = max_iter
        self.bp_method = bp_method
        self.ms_scaling_factor = ms_scaling_factor
        self.schedule = schedule
        self.omp_thread_count = int(omp_thread_count)
        self.max_factor_states = int(max_factor_states)
        self.cluster_mode = cluster_mode
        self.sum_logical_mass = bool(sum_logical_mass)
        self.bp = None
        self.refiner = None
        self.last_main_correction = None
        self.last_stats = {}
        self._initial_llr = None
        self._dem = None

    def update_params(self, **params):
        allowed = {
            "error_rate", "max_iter", "bp_method", "ms_scaling_factor", "schedule",
            "omp_thread_count", "max_factor_states", "cluster_mode", "sum_logical_mass",
        }
        unknown = set(params) - allowed
        assert not unknown, f"Unknown RecursiveLoopCoreBP parameters: {sorted(unknown)}"
        if "max_factor_states" in params:
            assert isinstance(params["max_factor_states"], (int, np.integer)) and params["max_factor_states"] > 0
        if "cluster_mode" in params:
            assert params["cluster_mode"] in {"tree", "generic", "recursive"}
        if "schedule" in params:
            assert params["schedule"] in {"parallel", "serial", "serial_relative"}
        if "omp_thread_count" in params:
            assert isinstance(params["omp_thread_count"], (int, np.integer)) and params["omp_thread_count"] > 0
        super().update_params(**params)
        if self.bp is not None:
            self._init_backend(self._dem)

    def _logical_matrix(self) -> np.ndarray:
        if self.mode == "dem":
            return _dense_binary(self.L)
        n = self.code.n
        lz = _dense_binary(self.code.L_Z)
        lx = _dense_binary(self.code.L_X)
        if self.H.shape[1] == n:
            compatible = [matrix for matrix in (lz, lx) if matrix.shape[1] == n and matrix.shape[0]]
            return np.vstack(compatible) if compatible else np.zeros((0, n), dtype=np.uint8)
        assert self.H.shape[1] == 2 * n
        zeros_z = np.zeros((lz.shape[0], n), dtype=np.uint8)
        zeros_x = np.zeros((lx.shape[0], n), dtype=np.uint8)
        return np.block([[lz, zeros_z], [zeros_x, lx]])

    def _check_metadata(self, dem: stim.DetectorErrorModel | None):
        times = np.full(self.m, -1, dtype=np.int32)
        groups = np.full(self.m, -1, dtype=np.int32)
        layers = np.full(self.m, -1, dtype=np.int32)
        if dem is None:
            return times, groups, layers

        coordinates = dem.get_detector_coordinates()
        assert len(coordinates) == self.m
        coordinate_rows = {
            check: tuple(float(x) for x in raw_coord)
            for check, raw_coord in coordinates.items()
        }
        base_h = None
        base_layers = None
        coordinate_to_group = {}
        if self.code is not None and hasattr(self.code, "H"):
            base_h = _dense_binary(self.code.H)
            base_layers = np.full(base_h.shape[0], -1, dtype=np.int32)
            if hasattr(self.code, "rows"):
                variable_layer = np.full(base_h.shape[1], -1, dtype=np.int32)
                for layer, variables in enumerate(self.code.rows):
                    variable_layer[np.asarray(variables, dtype=np.int64)] = layer
                for check in range(base_h.shape[0]):
                    support = np.flatnonzero(base_h[check])
                    if support.size:
                        base_layers[check] = int(np.max(variable_layer[support]))
            if hasattr(self.code, "stabilizer_coords"):
                base_coordinates = [tuple(float(x) for x in row) for row in self.code.stabilizer_coords()]
                if len(base_coordinates) == base_h.shape[0] and len(set(base_coordinates)) == len(base_coordinates):
                    coordinate_to_group = {coord: row for row, coord in enumerate(base_coordinates)}

        fallback_groups = {}
        if not coordinate_to_group:
            for coord in coordinate_rows.values():
                if len(coord) >= 2 and coord[:-1] not in fallback_groups:
                    fallback_groups[coord[:-1]] = len(fallback_groups)

        for check, coord in coordinate_rows.items():
            if coord:
                time = int(round(coord[-1]))
                if abs(coord[-1] - time) < 1e-9:
                    times[check] = time
            group = coordinate_to_group.get(coord[:-1], -1)
            if group < 0 and not coordinate_to_group and len(coord) == 2 and base_h is not None:
                candidate = int(round(coord[0]))
                if abs(coord[0] - candidate) < 1e-9 and 0 <= candidate < base_h.shape[0]:
                    group = candidate
            if group < 0:
                group = fallback_groups.get(coord[:-1], -1)
            if group >= 0:
                groups[check] = group
                if base_layers is not None and 0 <= group < base_layers.size:
                    layers[check] = base_layers[group]
        return times, groups, layers

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._dem = dem
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
        times, groups, layers = self._check_metadata(dem)
        mode = {"tree": 0, "generic": 1, "recursive": 2}[self.cluster_mode]
        self.refiner = RecursiveLoopCoreRefiner(
            _dense_binary(self.H),
            self._logical_matrix(),
            priors,
            times,
            groups,
            layers,
            self.max_factor_states,
            mode,
            self.sum_logical_mass,
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
