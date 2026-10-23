from typing import Literal

import numpy as np
import stim
from ldpc_modified import MbbpLdDecoder

from ..base import Decoder


class MBBPLD(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        error_channel=None,
        max_iter: int = 100,
        bp_method: Literal["product_sum", "minimum_sum"] = "minimum_sum",
        ms_scaling_factor: float = 0.875,
        schedule: Literal["serial", "parallel"] = "serial",
        check_permutations=None,
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.error_channel = error_channel
        if error_channel is not None:
            self.error_rate = None
        self.max_iter = max_iter
        self.bp_method = bp_method
        self.ms_scaling_factor = ms_scaling_factor
        self.schedule = schedule
        self.check_permutations = check_permutations
        self.bp = None

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)
        self.bp = MbbpLdDecoder(
            pcm=self.H,
            error_rate=(
                self.error_rate
                if self.mode == "static" and self.error_channel is None
                else None
            ),
            error_channel=(
                list(self.priors)
                if self.mode == "dem"
                else self.error_channel
            ),
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            schedule=self.schedule,
            check_permutations=self.check_permutations,
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp is None:
            self._init_backend()
        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        return self.bp.decode(syndrome)

    def converge(self) -> bool:
        if self.bp is None:
            self._init_backend()
        return self.bp.converge
