from typing import Literal
import numpy as np
import stim
from ldpc import BpLsdDecoder

from ..base import Decoder


class BPLSD(Decoder):

    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 60,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        lsd_order: int = 10,
        lsd_method: str = "lsd_cs",
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.max_iter          = max_iter
        self.bp_method         = bp_method
        self.ms_scaling_factor = ms_scaling_factor
        self.lsd_order         = lsd_order
        self.lsd_method        = lsd_method
        self.bp_lsd = None

        # if self.code is not None:
        #     self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)
            
        self.bp_lsd = BpLsdDecoder(
            pcm=self.H,
            error_rate=self.error_rate,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            lsd_method=self.lsd_method,
            lsd_order=self.lsd_order
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.bp_lsd is None:
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        return self.bp_lsd.decode(syndrome).astype(np.uint8)


