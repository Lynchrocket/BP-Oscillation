"""
BP+OTF Decoder - Belief Propagation with Ordered Tanner forest

An almost-linear time decoding algorithm for quantum LDPC codes under circuit-level noise (Antonio et al.)
"""

import numpy as np
import scipy.sparse as sp
from typing import Literal, Tuple
import BPOTF as OTF
import stim

from ..base import Decoder


class BPOTF(Decoder):
    def __init__(
        self, 
        code=None,
        error_rate=0.001,
        max_iters_3s: Tuple[int, int, int] = (30, 100, 100),
        bp_method: Literal["minimum_sum", "product_sum"] = "minimum_sum", 
        **kwargs
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)

        self.max_iters_3s = max_iters_3s
        self.bp_method = bp_method
        self.otf = None

        # if self.code is not None:
        #     self._init_backend()

    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse(self.code, dem)

        H = self.H.toarray() if sp.issparse(self.H) else self.H
        self.otf = OTF.OBPOTF(
            np.array(H).astype(np.uint8),
            p=self.error_rate,
            po_ext_bp_iters=np.array(self.max_iters_3s).astype(np.int32),
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.otf is None:
            self._init_backend()

        syndrome = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        return self.otf.decode(syndrome.astype(np.uint8)).astype(np.uint8)