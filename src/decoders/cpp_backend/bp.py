from typing import Literal
import numpy as np
import stim
from ldpc import BpDecoder

from ..base import Decoder


class BP(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.001,
        max_iter: int = 30,
        bp_method: Literal['product_sum', 'minimum_sum'] = "minimum_sum",
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.max_iter   = max_iter
        self.bp_method  = bp_method

        self.bp = None

        # if self.code is not None:
        #     self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)

        self.bp = BpDecoder(
            pcm=self.H,
            error_rate=self.error_rate,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
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