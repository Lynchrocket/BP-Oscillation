from abc import ABC, abstractmethod
import numpy as np
import scipy.sparse as sp
import stim
import sinter
import numba as nb


class Decoder(ABC, sinter.Decoder):
    """Abstract base class providing a unified interface for all decoders.

    Parameters
    ----------
    code : CSS-like object
        Must expose H_X, H_Z, L_X, L_Z as numpy arrays or sparse matrices,
        along with integer attributes n, k, c_x, c_z.
    error_rate : float
        Physical error rate used to initialise channel log-likelihood ratios.
    """

    def __init__(self, code=None, error_rate: float | None = None, **kwargs):
        self.code       = code
        self.error_rate = error_rate
        self.config     = dict(kwargs)
        self.mode       = None

    def update_params(self, **params):
        """Inject updated hyperparameters at runtime."""
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.config[key] = value

    @abstractmethod
    def _init_backend(self, dem: stim.DetectorErrorModel = None) -> None:
        """Initialise the internal decoding backend."""

    @abstractmethod
    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        """Decode a single syndrome vector and return a corresponding correction."""

    def _parse_pcm(self, code):
        """Parse the parity check matrix from the code wrapper."""
        if hasattr(code, 'H_X') and hasattr(code, 'H_Z'):
            hx = code.H_X
            hz = code.H_Z
            mx, nx = hx.shape
            mz, nz = hz.shape
            H = np.block([
                [hz, np.zeros((mz, nx))],
                [np.zeros((mx, nz)), hx]
            ])
        elif hasattr(code, 'H'):
            H = code.H
        else:
            raise ValueError("Code wrapper missing H matrix.")
        
        self.mode = "static"
        self.H = H.astype(np.uint8) % 2
        self.L_X = code.L_X
        self.L_Z = code.L_Z

    def _parse_dem(self, dem: stim.DetectorErrorModel) -> np.ndarray:
        """Sparse fault matrix from DEM."""
        from beliefmatching import detector_error_model_to_check_matrices
        self.mode = "dem"
        demMatrices = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
        self.H = demMatrices.check_matrix
        self.L = demMatrices.observables_matrix
        self.num_errors = dem.num_errors
        self.num_detectors = dem.num_detectors
        self.num_observables = dem.num_observables
        self.priors = demMatrices.priors

    def _parse(self, code=None, dem=None):
        assert code is not None or dem is not None, "Either code or DEM must be provided."

        if dem is not None:
            self._parse_dem(dem)
        else:
            self._parse_pcm(code)
        
        self.m, self.n = self.H.shape
        self.H = self.H.toarray() if sp.issparse(self.H) else self.H
        self.var_conn = [np.flatnonzero(self.H[:, j]) for j in range(self.n)]
        self.check_conn = [np.flatnonzero(self.H[i, :]) for i in range(self.m)]

    def _logical_from_correction(self, correction: np.ndarray) -> np.ndarray:
        """Map a correction to logical observable predictions."""
        correction = np.asarray(correction, dtype=np.uint8).reshape(-1) % 2

        if self.mode == "dem":
            all_obs = (self.L @ correction) % 2
            all_obs = all_obs[:self.num_observables]
        else:
            # correction[:n]   corresponds to X-type faults constrained by H_Z
            # correction[n:2n] corresponds to Z-type faults constrained by H_X
            n = self.code.n
            logical_x = (self.code.L_Z @ correction[:n]) % 2
            logical_z = (self.code.L_X @ correction[n:]) % 2
            all_obs = np.concatenate([logical_x, logical_z]).astype(np.uint8)

        return all_obs
    
    def decode(self, syndrome: np.ndarray) -> np.ndarray:
        """Decode a single syndrome vector and return observable predictions."""
        correction = self.decode_to_correction(syndrome)
        return self._logical_from_correction(correction)
    
    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> sinter.CompiledDecoder:
        """Return a compiled decoder for the given Detector Error Model."""
        return _SinterCompiledWrapper(self, dem)

    def compile_decoder_for_pcm(self, *, pcm: np.ndarray):
        self.H = pcm
        self.m, self.n = self.H.shape
        self.var_conn = [np.flatnonzero(self.H[:, j]) for j in range(self.n)]
        self.check_conn = [np.flatnonzero(self.H[i, :]) for i in range(self.m)]
        return self

class _SinterCompiledWrapper(sinter.CompiledDecoder):
    """Sinter CompiledDecoder that delegates to a Python Decoder instance.

    Handles the bit-packing/unpacking required by the sinter interface and
    routes each unpacked shot to the decoder's decode() method.
    This wrapper is used by all decoders that do not provide a native
    sinter batch interface.
    """

    def __init__(self, decoder: Decoder, dem: stim.DetectorErrorModel):
        self.decoder = decoder
        self.num_errors = dem.num_errors
        self.num_detectors = dem.num_detectors
        self.num_observables = dem.num_observables
        # Initialise the backend once here, in the sinter worker process.
        self.decoder._init_backend(dem=dem)

    def decode_shots_bit_packed(
        self, *, bit_packed_detection_event_data: np.ndarray
    ) -> np.ndarray:
        """Decode a batch of bit-packed shots and return bit-packed predictions.

        Parameters
        ----------
        bit_packed_detection_event_data : np.ndarray
            Shape (num_shots, ceil(num_detectors / 8)), dtype uint8.

        Returns
        -------
        np.ndarray
            Shape (num_shots, ceil(num_observables / 8)), dtype uint8.
        """
        shots = np.unpackbits(
            bit_packed_detection_event_data,
            axis=1,
            count=self.num_detectors,
            bitorder="little",
        )
        num_shots = shots.shape[0]
        predictions = np.zeros((num_shots, self.num_observables), dtype=np.uint8)
        for i in range(num_shots):
            predictions[i, :] = self.decoder.decode(shots[i, :])
        return np.packbits(predictions, axis=1, bitorder="little")
    
    def decode_with_correction_shots_bit_packed(
        self, *, bit_packed_detection_event_data: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        shots = np.unpackbits(
            bit_packed_detection_event_data,
            axis=1,
            count=self.num_detectors,
            bitorder="little",
        )
        num_shots = shots.shape[0]
        corrections = np.zeros((num_shots, self.num_errors), dtype=np.uint8)
        predictions = np.zeros((num_shots, self.num_observables), dtype=np.uint8)
        for i in range(num_shots):
            corrections[i, :] = self.decoder.decode_to_correction(shots[i, :])
            predictions[i, :] = self.decoder._logical_from_correction(corrections[i, :])
        return np.packbits(corrections, axis=1, bitorder="little"), np.packbits(predictions, axis=1, bitorder="little")
    