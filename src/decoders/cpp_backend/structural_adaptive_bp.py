from typing import Literal

import numpy as np
import scipy.sparse as sp
import stim
from ldpc_modified import StructuralAdaptiveBpDecoder

from ..base import Decoder


def _dense_binary(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.uint8) % 2


class StructuralAdaptiveBP(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal["product_sum", "minimum_sum"] = "minimum_sum",
        adaptive_max_views: int = 4,
        view_thread_count: int = 0,
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        assert isinstance(adaptive_max_views, (int, np.integer))
        assert isinstance(view_thread_count, (int, np.integer))
        assert adaptive_max_views >= 0 and view_thread_count >= 0
        self.max_iter = max_iter
        self.bp_method = bp_method
        self.adaptive_max_views = adaptive_max_views
        self.view_thread_count = view_thread_count
        self.bp = None

    def update_params(self, **params):
        if "adaptive_max_views" in params:
            assert isinstance(params["adaptive_max_views"], (int, np.integer))
            assert params["adaptive_max_views"] >= 0
        if "view_thread_count" in params:
            assert isinstance(params["view_thread_count"], (int, np.integer))
            assert params["view_thread_count"] >= 0
        super().update_params(**params)
        if self.bp is None:
            return
        if "max_iter" in params:
            self.bp.max_iter = self.max_iter
        if "bp_method" in params:
            self.bp.bp_method = self.bp_method
        if "error_rate" in params and self.mode == "static":
            self.bp.error_rate = self.error_rate
        if "adaptive_max_views" in params:
            self.bp.adaptive_max_views = self.adaptive_max_views
        if "view_thread_count" in params:
            self.bp.view_thread_count = self.view_thread_count

    def _base_structure(self):
        if self.code is None or not hasattr(self.code, "H"):
            return None, None
        base_h = _dense_binary(self.code.H)
        check_layers = np.full(base_h.shape[0], -1, dtype=np.int32)
        if hasattr(self.code, "rows"):
            qubit_layer = np.full(base_h.shape[1], -1, dtype=np.int32)
            for layer, variables in enumerate(self.code.rows):
                qubit_layer[np.asarray(variables, dtype=np.int64)] = layer
            for check in range(base_h.shape[0]):
                support = np.flatnonzero(base_h[check])
                if support.size:
                    check_layers[check] = int(np.max(qubit_layer[support]))
        return base_h, check_layers

    def _structure_metadata(self, dem: stim.DetectorErrorModel | None):
        check_count = self.H.shape[0]
        neighbors = [set() for _ in range(check_count)]
        times = np.full(check_count, -1, dtype=np.int32)
        groups = np.full(check_count, -1, dtype=np.int32)
        layers = np.full(check_count, -1, dtype=np.int32)
        base_h, base_layers = self._base_structure()

        if dem is not None:
            coordinates = dem.get_detector_coordinates()
            assert len(coordinates) == check_count
            coordinate_rows = [tuple(float(x) for x in coordinates[row]) for row in range(check_count)]
            coordinate_to_group = {}
            if base_h is not None and hasattr(self.code, "stabilizer_coords"):
                base_coordinates = [tuple(float(x) for x in coord) for coord in self.code.stabilizer_coords()]
                if len(base_coordinates) == base_h.shape[0] and len(set(base_coordinates)) == len(base_coordinates):
                    coordinate_to_group = {coord: row for row, coord in enumerate(base_coordinates)}

            spatial_keys = []
            for row, coord in enumerate(coordinate_rows):
                group = coordinate_to_group.get(coord[:-1], -1)
                if group < 0 and not coordinate_to_group and len(coord) == 2 and base_h is not None:
                    candidate = int(round(coord[0]))
                    if abs(coord[0] - candidate) < 1e-9 and 0 <= candidate < base_h.shape[0]:
                        group = candidate
                if group >= 0:
                    spatial_keys.append((group,))
                    groups[row] = group
                    layers[row] = base_layers[group]
                    time = int(round(coord[-1]))
                    if abs(coord[-1] - time) < 1e-9:
                        times[row] = time
                else:
                    spatial_keys.append(coord)

            temporal_buckets = {}
            for row, key in enumerate(spatial_keys):
                if times[row] >= 0:
                    temporal_buckets.setdefault(key, []).append(row)
            for rows in temporal_buckets.values():
                rows.sort(key=lambda row: (times[row], row))
                for left, right in zip(rows, rows[1:]):
                    if times[left] != times[right]:
                        neighbors[left].add(right)
                        neighbors[right].add(left)

            if base_h is not None:
                wide_base_h = base_h.astype(np.int32, copy=False)
                base_overlap = (wide_base_h @ wide_base_h.T) > 0
                time_buckets = {}
                for row, time in enumerate(times):
                    if time >= 0 and groups[row] >= 0:
                        time_buckets.setdefault(int(time), []).append(row)
                for rows in time_buckets.values():
                    for i, left in enumerate(rows):
                        for right in rows[i + 1 :]:
                            if base_overlap[groups[left], groups[right]]:
                                neighbors[left].add(right)
                                neighbors[right].add(left)

        # Structural edges augment the complete Tanner check graph.
        dense_h = _dense_binary(self.H)
        wide_h = dense_h.astype(np.int32, copy=False)
        check_overlap = (wide_h @ wide_h.T) > 0
        for left in range(check_count):
            for right in np.flatnonzero(check_overlap[left, left + 1 :]) + left + 1:
                neighbors[left].add(int(right))
                neighbors[int(right)].add(left)

        has_time = bool(np.any(times >= 0))
        has_code_groups = bool(np.any(groups >= 0))
        mode = 3 if has_time and has_code_groups else 2 if dem is not None else 1
        return (
            [sorted(row) for row in neighbors],
            times.tolist(),
            groups.tolist(),
            layers.tolist(),
            mode,
        )

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)
        self.bp = StructuralAdaptiveBpDecoder(
            pcm=self.H,
            error_rate=self.error_rate if self.mode == "static" else None,
            error_channel=list(self.priors) if self.mode == "dem" else None,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            adaptive_max_views=self.adaptive_max_views,
            view_thread_count=self.view_thread_count,
            input_vector_type="syndrome",
        )
        self.bp.set_structure_metadata(*self._structure_metadata(dem))

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()
        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        return self.bp.decode(syndrome).astype(np.uint8)

    def converge(self) -> bool:
        if self.bp is None:
            self._init_backend()
        return self.bp.converge

    def auxiliary_converge(self) -> bool:
        if self.bp is None:
            self._init_backend()
        return self.bp.auxiliary_converge
