"""
Cython-backed search-based decoder.

This decoder follows the local cython convention H = diag(H_X, H_Z).
"""

import numpy as np
import stim
from typing import Literal

from ..base import Decoder

from .cython_core.tesseract_core import search_sector_core


class Tesseract(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        beam: int = 10,
        pq_limit: int = 1000000,
        det_penalty: float = 0.0,
        no_revisit_dets: bool = False,
        low_confidence_action: Literal["assert", "return_zero"] = "assert",
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        assert beam >= 0
        assert pq_limit >= 1
        assert low_confidence_action in ("assert", "return_zero")

        self.beam = beam
        self.pq_limit = pq_limit
        self.detector_penalty = det_penalty
        self.no_revisit_dets = no_revisit_dets
        self.low_confidence_action = low_confidence_action
        self.num_observables = None
        self.search_paths = {"x": [], "z": []}
        self.search_costs = {"x": np.inf, "z": np.inf}
        self.low_confidence = {"x": False, "z": False}
        self._graph_ready = False

        if self.code is not None:
            self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None) -> None:
        assert dem is None, "Cython Tesseract currently supports code-level decoding only"
        assert self.code is not None, "code is required to initialize Cython Tesseract"

        self.H_X = np.ascontiguousarray(self.code.H_X, dtype=np.uint8) % 2
        self.H_Z = np.ascontiguousarray(self.code.H_Z, dtype=np.uint8) % 2
        self.c_x, n_x = self.H_X.shape
        self.c_z, n_z = self.H_Z.shape
        assert n_x == n_z
        p = float(self.error_rate if self.error_rate is not None else 1e-3)
        self._search_error_rate = min(max(p, 1e-12), 1.0 - 1e-12)
        self.n = n_x
        self.m = self.c_x + self.c_z
        self.H = np.block(
            [
                [self.H_Z, np.zeros((self.c_x, self.n), dtype=np.uint8)],
                [np.zeros((self.c_z, self.n), dtype=np.uint8), self.H_X],
            ]
        ).astype(np.uint8)
        self._graph_ready = True

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if not self._graph_ready:
            self._init_backend()

        syndrome = np.ascontiguousarray(np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2)
        assert syndrome.size == self.m, f"syndrome has length {syndrome.size}, expected {self.m}"

        syndrome_x = np.ascontiguousarray(syndrome[: self.c_x], dtype=np.uint8)
        syndrome_z = np.ascontiguousarray(syndrome[self.c_x :], dtype=np.uint8)

        # With the local cython convention, the first syndrome block belongs to
        # H_X and the second block belongs to H_Z, matching H = diag(H_X, H_Z).
        correction_x, cost_x, path_x, low_x = search_sector_core(
            self.H_X,
            syndrome_x,
            self._search_error_rate,
            int(self.beam),
            int(self.pq_limit),
            float(self.detector_penalty),
            bool(self.no_revisit_dets),
        )
        correction_z, cost_z, path_z, low_z = search_sector_core(
            self.H_Z,
            syndrome_z,
            self._search_error_rate,
            int(self.beam),
            int(self.pq_limit),
            float(self.detector_penalty),
            bool(self.no_revisit_dets),
        )

        self.search_paths = {"x": path_x, "z": path_z}
        self.search_costs = {"x": cost_x, "z": cost_z}
        self.low_confidence = {"x": low_x, "z": low_z}
        correction = np.concatenate(
            [
                np.asarray(correction_x, dtype=np.uint8),
                np.asarray(correction_z, dtype=np.uint8),
            ]
        )
        if low_x or low_z:
            if self.low_confidence_action == "return_zero":
                return np.zeros(2 * self.n, dtype=np.uint8)
            raise AssertionError("Tesseract search ended with low confidence")

        assert np.array_equal((self.H @ correction) % 2, syndrome), "Tesseract search returned an inconsistent correction"
        return correction
