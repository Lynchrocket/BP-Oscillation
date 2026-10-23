from typing import Literal

import numpy as np
import scipy.sparse as sp
import stim
from ldpc_modified import BpDecoder, RecursiveFrontierRefiner

from ..base import Decoder


def _dense_binary(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    return np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) % 2)


class RecursiveFrontierBP(Decoder):
    """BP followed by bounded inference over recursive space-time separators."""

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["product_sum", "minimum_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        schedule: Literal["parallel", "serial", "serial_relative"] = "parallel",
        omp_thread_count: int = 1,
        max_frontier_states: int = 4096,
        ordering_mode: Literal["recursive", "generic", "shuffled"] = "recursive",
        ordering_seed: int = 20260902,
        **kwargs,
    ):
        assert not kwargs, f"Unknown RecursiveFrontierBP parameters: {sorted(kwargs)}"
        assert code is not None and hasattr(code, "rows"), (
            "RecursiveFrontierBP requires a structured code with recursive rows"
        )
        super().__init__(code=code, error_rate=error_rate)
        assert isinstance(max_frontier_states, (int, np.integer)) and 0 < max_frontier_states <= np.iinfo(np.int32).max
        assert ordering_mode in {"recursive", "generic", "shuffled"}
        assert schedule in {"parallel", "serial", "serial_relative"}
        assert isinstance(omp_thread_count, (int, np.integer)) and omp_thread_count > 0
        assert isinstance(ordering_seed, (int, np.integer))
        self.max_iter = max_iter
        self.bp_method = bp_method
        self.ms_scaling_factor = ms_scaling_factor
        self.schedule = schedule
        self.omp_thread_count = int(omp_thread_count)
        self.max_frontier_states = int(max_frontier_states)
        self.ordering_mode = ordering_mode
        self.ordering_seed = int(ordering_seed)
        self.bp = None
        self.refiner = None
        self.last_main_correction = None
        self.last_stats = {}
        self.variable_order = None
        self.bag_ends = None
        self._initial_llr = None
        self._dem = None
        self._pcm_override = None

    def update_params(self, **params):
        allowed = {
            "error_rate", "max_iter", "bp_method", "ms_scaling_factor", "schedule",
            "omp_thread_count", "max_frontier_states", "ordering_mode", "ordering_seed",
        }
        unknown = set(params) - allowed
        assert not unknown, f"Unknown RecursiveFrontierBP parameters: {sorted(unknown)}"
        if "max_frontier_states" in params:
            assert isinstance(params["max_frontier_states"], (int, np.integer))
            assert 0 < params["max_frontier_states"] <= np.iinfo(np.int32).max
        if "ordering_mode" in params:
            assert params["ordering_mode"] in {"recursive", "generic", "shuffled"}
        if "schedule" in params:
            assert params["schedule"] in {"parallel", "serial", "serial_relative"}
        if "omp_thread_count" in params:
            assert isinstance(params["omp_thread_count"], (int, np.integer))
            assert params["omp_thread_count"] > 0
        if "ordering_seed" in params:
            assert isinstance(params["ordering_seed"], (int, np.integer))
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
        return np.block([
            [lz, np.zeros((lz.shape[0], n), dtype=np.uint8)],
            [np.zeros((lx.shape[0], n), dtype=np.uint8), lx],
        ])

    def _logical_from_correction(self, correction: np.ndarray) -> np.ndarray:
        if self.mode == "dem" or self.H.shape[1] == 2 * self.code.n:
            return super()._logical_from_correction(correction)
        correction = np.asarray(correction, dtype=np.uint8).reshape(-1) % 2
        assert correction.size == self.code.n
        return (self._logical_matrix() @ correction) % 2

    def _base_structure(self):
        assert self.code is not None and hasattr(self.code, "rows")
        assert hasattr(self.code, "H"), "RecursiveFrontierBP requires the base code.H"
        base_h = _dense_binary(self.code.H)
        variable_layer = np.full(base_h.shape[1], -1, dtype=np.int32)
        variable_position = np.full(base_h.shape[1], -1, dtype=np.int32)
        for layer, variables in enumerate(self.code.rows):
            for position, variable in enumerate(variables):
                assert 0 <= variable < base_h.shape[1]
                assert variable_layer[variable] < 0
                variable_layer[variable] = layer
                variable_position[variable] = position
        assert np.all(variable_layer >= 0), "code.rows must partition all base variables"

        check_layer = np.full(base_h.shape[0], -1, dtype=np.int32)
        check_position = np.full(base_h.shape[0], -1, dtype=np.int32)
        next_position = {}
        for check in range(base_h.shape[0]):
            support = np.flatnonzero(base_h[check])
            assert support.size > 0
            layer = int(np.max(variable_layer[support]))
            check_layer[check] = layer
            check_position[check] = next_position.get(layer, 0)
            next_position[layer] = int(check_position[check]) + 1
        return base_h, variable_layer, variable_position, check_layer, check_position

    def _detector_structure(self, dem: stim.DetectorErrorModel | None):
        base_h, variable_layer, variable_position, base_layers, base_positions = self._base_structure()
        if dem is None:
            assert self.H.shape == base_h.shape
            return (
                np.zeros(self.m, dtype=np.int32),
                base_layers,
                base_positions,
                variable_layer,
                variable_position,
            )

        times = np.full(self.m, -1, dtype=np.int32)
        layers = np.full(self.m, -1, dtype=np.int32)
        positions = np.full(self.m, -1, dtype=np.int32)
        coordinates = dem.get_detector_coordinates()
        assert len(coordinates) == self.m
        base_coordinate_to_group = {}
        if hasattr(self.code, "stabilizer_coords"):
            base_coordinates = [tuple(float(x) for x in row) for row in self.code.stabilizer_coords()]
            assert len(base_coordinates) == base_h.shape[0]
            assert len(set(base_coordinates)) == len(base_coordinates)
            base_coordinate_to_group = {coord: index for index, coord in enumerate(base_coordinates)}

        for check in range(self.m):
            coord = tuple(float(x) for x in coordinates[check])
            assert len(coord) >= 2, "detectors require a spatial coordinate and time"
            time = int(round(coord[-1]))
            assert abs(coord[-1] - time) < 1e-9
            times[check] = time
            if base_coordinate_to_group:
                group = base_coordinate_to_group.get(coord[:-1], -1)
            else:
                assert len(coord) == 2
                group = int(round(coord[0]))
                assert abs(coord[0] - group) < 1e-9
            assert 0 <= group < base_h.shape[0]
            layers[check] = base_layers[group]
            positions[check] = base_positions[group]
        return times, layers, positions, variable_layer, variable_position

    @staticmethod
    def _order_and_bag_ends(keys):
        ordered = sorted(keys)
        order = np.asarray([key[-1] for key in ordered], dtype=np.int32)
        bag_ends = []
        for position, key in enumerate(ordered):
            if position + 1 == len(ordered) or ordered[position + 1][:-1] != key[:-1]:
                bag_ends.append(position)
        return order, np.asarray(bag_ends, dtype=np.int32)

    def _recursive_order(self, dem: stim.DetectorErrorModel | None):
        times, layers, positions, variable_layers, variable_positions = self._detector_structure(dem)
        if dem is None:
            keys = [
                (int(variable_layers[v]), int(variable_positions[v]), v)
                for v in range(self.n)
            ]
            return self._order_and_bag_ends(keys)

        keys = []
        logical_only_key = (
            int(np.max(times, initial=-1)) + 1,
            int(np.max(layers, initial=-1)) + 1,
            int(np.max(positions, initial=-1)) + 1,
            0,
            0,
            0,
            0,
        )
        for variable, checks in enumerate(self.var_conn):
            checks = np.asarray(checks, dtype=np.int64)
            if checks.size == 0:
                keys.append((*logical_only_key, variable))
                continue
            variable_times = times[checks]
            variable_layers = layers[checks]
            variable_positions = positions[checks]
            keys.append((
                int(np.max(variable_times)),
                int(np.max(variable_layers)),
                int(np.max(variable_positions)),
                int(np.ptp(variable_times)),
                int(np.ptp(variable_layers)),
                int(np.ptp(variable_positions)),
                int(checks.size),
                variable,
            ))
        return self._order_and_bag_ends(keys)

    def _generic_order(self) -> np.ndarray:
        remaining = np.ones(self.n, dtype=bool)
        remaining_degree = np.asarray([len(checks) for checks in self.check_conn], dtype=np.int32)
        opened = np.zeros(self.m, dtype=bool)
        active_count = 0
        order = []
        for _ in range(self.n):
            best_key = None
            best_variable = None
            for variable in np.flatnonzero(remaining):
                checks = self.var_conn[variable]
                new_count = int(np.count_nonzero(~opened[checks]))
                close_count = int(np.count_nonzero(remaining_degree[checks] == 1))
                key = (active_count + new_count - close_count, -close_count, new_count, len(checks), int(variable))
                if best_key is None or key < best_key:
                    best_key = key
                    best_variable = int(variable)
            order.append(best_variable)
            remaining[best_variable] = False
            for check in self.var_conn[best_variable]:
                if not opened[check]:
                    opened[check] = True
                    active_count += 1
                remaining_degree[check] -= 1
                if remaining_degree[check] == 0:
                    active_count -= 1
        return np.asarray(order, dtype=np.int32)

    def _variable_order(self, dem: stim.DetectorErrorModel | None):
        if self.ordering_mode == "generic":
            order = self._generic_order()
            return order, np.arange(self.n, dtype=np.int32)
        recursive, bag_ends = self._recursive_order(dem)
        if self.ordering_mode == "shuffled":
            recursive = recursive.copy()
            np.random.default_rng(self.ordering_seed).shuffle(recursive)
            bag_ends = np.arange(self.n, dtype=np.int32)
        return recursive, bag_ends

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._dem = dem
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
        self.variable_order, self.bag_ends = self._variable_order(dem)
        self.refiner = RecursiveFrontierRefiner(
            _dense_binary(self.H),
            self._logical_matrix(),
            priors,
            self.variable_order,
            self.bag_ends,
            self.max_frontier_states,
        )

    def compile_decoder_for_pcm(self, *, pcm: np.ndarray):
        self._pcm_override = _dense_binary(pcm)
        self._dem = None
        self.bp = None
        self.refiner = None
        self.mode = "static"
        return super().compile_decoder_for_pcm(pcm=self._pcm_override)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()
        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        assert syndrome.size == self.m
        main = np.asarray(self.bp.decode(syndrome), dtype=np.uint8)
        self.last_main_correction = main.copy()
        posterior_llr = self.bp.log_prob_ratios if np.any(syndrome) else self._initial_llr
        correction = self.refiner.refine(syndrome, main, posterior_llr)
        self.last_stats = dict(self.refiner.stats)
        self.last_stats["main_bp_converged"] = bool(self.bp.converge)
        return correction

    def converge(self) -> bool:
        return bool(self.last_stats.get("syndrome_valid", False))
