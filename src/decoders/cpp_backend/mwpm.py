import numpy as np
import pymatching
import stim

from ..base import Decoder


class MWPM(Decoder):
    def __init__(
        self, 
        code=None, 
        error_rate: float = 0.1, 
        **kwargs
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)

        self.matcher_x          = None
        self.matcher_z          = None
        self.matching_available = False
        self.matching_error     = None

        # if self.code is not None:
        #     self._init_backend()
        
    def _init_backend(self, dem: stim.DetectorErrorModel = None):
        self._parse_pcm(self.code)

        weight = -np.log(self.error_rate / (1 - self.error_rate))
        try:
            if hasattr(self.code, "H_X") and hasattr(self.code, "H_Z"):
                # Matcher for X errors (using H_Z)
                self.matcher_x = pymatching.Matching(self.code.H_Z, weights=weight)
                # Matcher for Z errors (using H_X)
                self.matcher_z = pymatching.Matching(self.code.H_X, weights=weight)
                self.matching_available = True
            else:
                self.matching_available = False
                self.matching_error = "Code must have H_X and H_Z matrices for MWPM decoding."
        except Exception as e:
            self.matching_available = False
            self.matching_error = str(e)
            
    
    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        if not self.matching_available:
            if self.matching_error is None:
                self._init_backend()
            else:
                raise RuntimeError(
                    f"MWPM is unavailable for this code. Reason: {self.matching_error}"
                )
        
        correction_x = self.matcher_x.decode(syndrome[:self.c_z])
        correction_z = self.matcher_z.decode(syndrome[self.c_z:])
        
        return np.concatenate([correction_x, correction_z]).astype(np.uint8)
    