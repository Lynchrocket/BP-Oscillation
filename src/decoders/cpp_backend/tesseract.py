import numpy as np
import sinter
import stim

from tesseract_decoder import TesseractSinterDecoder

from ..base import Decoder


class Tesseract(Decoder):
    def __init__(
        self,
        code=None,
        error_rate: float = 0.1,
        beam: int = 5,
        pq_limit: int = 500,
        num_det_orders: int = 1,
        **kwargs,
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        self.beam           = beam
        self.pq_limit       = pq_limit
        self.num_det_orders = num_det_orders
        self.seed           = kwargs.get("seed", 0)
        self.verbose        = kwargs.get("verbose", False)
        self.create_visualization = kwargs.get("create_visualization", False)

    def _init_backend(self, dem: stim.DetectorErrorModel = None) -> None:
        """No-op"""

    def compile_decoder_for_dem(
        self, *, dem: stim.DetectorErrorModel
    ) -> sinter.CompiledDecoder:
        sinter_decoder = TesseractSinterDecoder(
            det_beam=self.beam,
            pqlimit=self.pq_limit,
            beam_climbing=True,
            no_revisit_dets=True,
            merge_errors=True,
            num_det_orders=self.num_det_orders,
            verbose=self.verbose,
            create_visualization=self.create_visualization,
            seed=self.seed,
        )
        return sinter_decoder.compile_decoder_for_dem(dem=dem)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        raise NotImplementedError(
            "Use sinter.collect() with this decoder instead."
        )