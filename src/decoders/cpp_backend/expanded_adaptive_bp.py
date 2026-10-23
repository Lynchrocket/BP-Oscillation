from typing import Literal
import numpy as np
import stim
from ldpc_modified import ExpandedAdaptiveBpDecoder

from ..base import Decoder


class ExpandedAdaptiveBP(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal['product_sum', 'minimum_sum'] = "minimum_sum",
        adaptive_max_views: int = 4,
        view_thread_count: int = 0,
        **kwargs,
    ):

        super().__init__(code=code, error_rate=error_rate, **kwargs)
        assert isinstance(adaptive_max_views, (int, np.integer))
        assert isinstance(view_thread_count, (int, np.integer))
        assert adaptive_max_views >= 0 and view_thread_count >= 0
        self.max_iter   = max_iter
        self.bp_method  = bp_method
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
        if self.bp is not None:
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

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)

        self.bp = ExpandedAdaptiveBpDecoder(
            pcm=self.H,
            error_rate=self.error_rate if self.mode == "static" else None,
            error_channel=list(self.priors) if self.mode == "dem" else None,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            adaptive_max_views=self.adaptive_max_views,
            view_thread_count=self.view_thread_count,
            input_vector_type="syndrome",
        )

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
