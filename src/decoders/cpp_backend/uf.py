import numpy as np
import scipy.sparse as sp
from typing import Literal
from ldpc import UnionFindDecoder

from ..base import Decoder


class UnionFind(Decoder):
    def __init__(
        self, 
        code=None, 
        error_rate=0.01,
        uf_method = "union_find", # wrong !! look at issue https://github.com/quantumgizmos/ldpc/issues/88
        **kwargs
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)

        
        self.uf_method = uf_method
        self.uf = None

        # if self.code is not None:
        #     self._init_backend()

    def _init_backend(self):
        self._parse_pcm()
        
        self.uf = UnionFindDecoder(
            self.H,
            uf_method=self.uf_method
        )

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if self.uf is None:
            self._init_backend()
        return self.uf.decode(syndrome)
