# failure_regions, force_osd_failures, low_weight_failures, quotient, trapping_sets
"""Ordinary decoder construction helper; callers supply decoder types explicitly."""
from decoders.cpp_backend import BPOSD, Tesseract

def build_decoder_pair(code, p: float, *, bposd_type=BPOSD, tesseract_type=Tesseract):
    bposd = bposd_type(
        code=code,
        error_rate=p,
        max_iter=60,
        bp_method='minimum_sum',
        osd_order=10,
        osd_method='osd_cs',
    )

    tess = tesseract_type(
        code=code,
        error_rate=p,
        beam=8,
        pq_limit=200000,
        det_penalty=0.0,
    )

    return bposd, tess

def build_bp_serial(code, p: float):
    bposd = BPOSD(
        code=code,
        error_rate=p,
        max_iter=60,
        bp_method='minimum_sum',
        osd_order=10,
        osd_method='osd_cs',
        schedule='serial',
    )

    tess = Tesseract(
        code=code,
        error_rate=p,
        beam=8,
        pq_limit=200000,
        det_penalty=0.0,
    )

    return bposd, tess

def build_force_osd(code, p: float):
    bposd = BPOSD(
        code=code,
        error_rate=p,
        max_iter=60,
        bp_method='minimum_sum',
        osd_order=10,
        osd_method='osd_cs',
        force_osd=True,
    )

    tess = Tesseract(
        code=code,
        error_rate=p,
        beam=8,
        pq_limit=200000,
        det_penalty=0.0,
    )

    return bposd, tess