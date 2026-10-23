from typing import Literal
import numpy as np
import stim
from ldpc_modified import AdaptiveBpOsdDecoder

from ..base import Decoder


class AdaptiveBPOSD(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 60,
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        osd_order: int = 10,
        osd_method: str = "osd_cs",
        schedule: Literal["parallel", "serial"] = "parallel",
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.max_iter          = max_iter
        self.bp_method         = bp_method
        self.osd_order         = osd_order
        self.osd_method        = osd_method
        self.ms_scaling_factor = ms_scaling_factor
        self.schedule = schedule

        self.adaptive_bp_osd = None
        self.correction = None

        # if self.code is not None:
        #     self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)

        self.adaptive_bp_osd = AdaptiveBpOsdDecoder(
            pcm=self.H,
            error_rate=self.error_rate if self.mode == "static" else None,
            error_channel=list(self.priors) if self.mode == "dem" else None,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            osd_method=self.osd_method,
            osd_order=self.osd_order,
            schedule=self.schedule,
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.adaptive_bp_osd is None:
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        self.correction = self.adaptive_bp_osd.decode(syndrome).astype(np.uint8)
        return self.correction

    def converge(self) -> bool:
        if self.adaptive_bp_osd is None:
            self._init_backend()

        return self.adaptive_bp_osd.converge

    def channel_probs(self) -> np.ndarray:
        if self.adaptive_bp_osd is None:
            self._init_backend()

        return self.adaptive_bp_osd.channel_probs
    
    def update_channel_probs(self, channel_probs: np.ndarray):
        if self.adaptive_bp_osd is None:
            self._init_backend()

        self.adaptive_bp_osd.update_channel_probs(channel_probs)

    def correction(self) -> np.ndarray:
        if self.adaptive_bp_osd is None:
            self._init_backend()

        if self.correction is None:
            raise ValueError("No correction available. Please run decode_to_correction first.")
        
        return self.correction