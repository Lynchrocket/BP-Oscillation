# failure_regions, force_osd_failures, low_weight_failures
import numpy as np
import pandas as pd
import stim
import networkx as nx
from collections import Counter
from typing import Literal, Callable
from itertools import combinations

from noises_ca import CircuitBuilder, ErrorRate
from decoders.cpp_backend import BPOSD, Tesseract

from utils.graph_struct import error_llr_weight

from .decoder_setup import build_decoder_pair, build_bp_serial, build_force_osd

def packed_rows_equal(a: np.ndarray, b: np.ndarray, n_bits: int) -> np.ndarray:
    """Row-wise equality of little-endian packed bit arrays, ignoring padding bits."""
    if n_bits == 0:
        return np.ones(a.shape[0], dtype=bool)

    if a.shape != b.shape:
        raise ValueError(f"Packed arrays have different shapes: {a.shape} vs {b.shape}")

    n_full = n_bits // 8
    rem = n_bits % 8

    ok = np.ones(a.shape[0], dtype=bool)

    if n_full:
        ok &= np.all(a[:, :n_full] == b[:, :n_full], axis=1)

    if rem:
        mask = (1 << rem) - 1
        ok &= ((a[:, n_full] ^ b[:, n_full]) & mask) == 0

    return ok

_POPCOUNT_U8 = np.array([int(i).bit_count() for i in range(256)], dtype=np.uint8)

def packed_weight(x: np.ndarray, n_bits: int) -> np.ndarray:
    """Hamming weight of little-endian packed rows, ignoring padding bits."""
    if n_bits == 0:
        return np.zeros(x.shape[0], dtype=np.int64)

    y = x.copy()
    rem = n_bits % 8
    if rem:
        y[:, -1] &= (1 << rem) - 1

    return _POPCOUNT_U8[y].sum(axis=1).astype(np.int64)

def packed_to_str(row: np.ndarray, n_bits: int) -> str:
    """Convert one little-endian packed row to a '0101...' string."""
    if n_bits == 0:
        return ""

    bits = np.unpackbits(np.asarray(row, dtype=np.uint8), bitorder="little")[:n_bits]

    return "".join("1" if b else "0" for b in bits)

def collect_samples(
    code_name: str,
    ca_code,
    error_rate: ErrorRate,
    noise_model: str,
    d_rounds: int = None,
    total_shots: int = 1_000_000_000,
    batch_shots: int = 100_000,
    sample_times: int = 1000,
    num_samples: int = 1000,
    kind: Literal['both_success', 'both_fail', 'bposd_fail_tess_success', 'bposd_success_tess_fail'] = 'bposd_fail_tess_success',
) -> tuple[list[dict], stim.DetectorErrorModel]:
    """
    Collect cases

    Depending on the noise model and error rate, it can be BPOSD-fail/Tesseract-success cases, or the opposite or both success/fail cases.
    """
    builder = CircuitBuilder(ca_code)
    circuit = builder.get_circuit(noise_model=noise_model, error_rate=error_rate, d_rounds=d_rounds)
    dem = circuit.detector_error_model()

    bposd, tess = build_decoder_pair(ca_code, error_rate.p)

    compiled_b = bposd.compile_decoder_for_dem(dem=dem)
    compiled_t = tess.compile_decoder_for_dem(dem=dem)

    sampler = circuit.compile_detector_sampler()
    num_det = dem.num_detectors
    num_obs = dem.num_observables

    recs: list[dict] = []
    sampled_shot = 0
    for batch_idx in range(sample_times):
        cur_shots = min(batch_shots, total_shots - sampled_shot)
        if cur_shots <= 0:
            break
        dets, obs = sampler.sample(shots=cur_shots, separate_observables=True, bit_packed=True)

        _, unique_idx = np.unique(dets, axis=0, return_index=True)
        if unique_idx.size == 0:
            sampled_shot += cur_shots
            continue

        dets = dets[unique_idx]
        obs = obs[unique_idx]

        pred_b = compiled_b.decode_shots_bit_packed(bit_packed_detection_event_data=dets)
        pred_t = compiled_t.decode_shots_bit_packed(bit_packed_detection_event_data=dets)

        b_succ = packed_rows_equal(pred_b, obs, num_obs)
        t_succ = packed_rows_equal(pred_t, obs, num_obs)

        match kind:
            case 'both_success':
                mask = b_succ & t_succ
            case 'both_fail':
                mask = (~b_succ) & (~t_succ)
            case 'bposd_fail_tess_success':
                mask = (~b_succ) & t_succ
            case 'bposd_success_tess_fail':
                mask = b_succ & (~t_succ)
            case _:
                raise ValueError(f'Unsupported kind: {kind}')
            
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            sampled_shot += cur_shots
            continue
        need = num_samples - len(recs)
        if idx.size > need:
            idx = idx[:need]
        
        dets_sel = dets[idx]
        obs_sel = obs[idx]
        pred_b_sel = pred_b[idx]
        pred_t_sel = pred_t[idx]

        weights = packed_weight(dets_sel, num_det)

        recs.extend([
            {
                'code': code_name,
                'p': error_rate.p,
                'p_meas': error_rate.p_meas,
                'p_cx': error_rate.p_cx,
                'p_prep': error_rate.p_prep,
                'num_detectors': int(num_det),
                'num_observables': int(num_obs),
                'detector_weight': int(w),
                'detector_events': packed_to_str(d, num_det),
                'true_observables': packed_to_str(o, num_obs),
                'bposd_pred_observables': packed_to_str(pb, num_obs),
                'tesseract_pred_observables': packed_to_str(pt, num_obs),
            }
            for w, d, o, pb, pt in zip(weights, dets_sel, obs_sel, pred_b_sel, pred_t_sel)
        ])

        sampled_shot += cur_shots
        
        if len(recs) >= num_samples:
            break

    return recs, dem


def count_samples_union(
    code_name: str,
    ca_code,
    error_rate: ErrorRate,
    noise_model: str,
    d_rounds: int = None,
    total_shots: int = 1_000_000_000,
    batch_shots: int = 100_000,
    sample_times: int = 1000,
    num_samples: int = 10_000,
) -> tuple[int, int, int, int, stim.DetectorErrorModel]:
    """
    Count samples in all different kinds
    """
    builder = CircuitBuilder(ca_code)
    circuit = builder.get_circuit(noise_model=noise_model, error_rate=error_rate, d_rounds=d_rounds)
    dem = circuit.detector_error_model()

    bposd, tess = build_decoder_pair(ca_code, error_rate.p)

    compiled_b = bposd.compile_decoder_for_dem(dem=dem)
    compiled_t = tess.compile_decoder_for_dem(dem=dem)

    sampler = circuit.compile_detector_sampler()
    num_det = dem.num_detectors
    num_obs = dem.num_observables

    bs_cnt, bf_cnt, bfts_cnt, bstf_cnt = 0, 0, 0, 0
    sampled_shot = 0
    for batch_idx in range(sample_times):
        cur_shots = min(batch_shots, total_shots - sampled_shot)
        if cur_shots <= 0:
            break
        dets, obs = sampler.sample(shots=cur_shots, separate_observables=True, bit_packed=True)

        _, unique_idx = np.unique(dets, axis=0, return_index=True)
        if unique_idx.size == 0:
            sampled_shot += cur_shots
            continue

        dets = dets[unique_idx]
        obs = obs[unique_idx]

        pred_b = compiled_b.decode_shots_bit_packed(bit_packed_detection_event_data=dets)
        pred_t = compiled_t.decode_shots_bit_packed(bit_packed_detection_event_data=dets)

        b_succ = packed_rows_equal(pred_b, obs, num_obs)
        t_succ = packed_rows_equal(pred_t, obs, num_obs)

        bs_cnt += np.sum(b_succ & t_succ)
        bf_cnt += np.sum((~b_succ) & (~t_succ))
        bfts_cnt += np.sum((~b_succ) & t_succ)
        bstf_cnt += np.sum(b_succ & (~t_succ))

        total_fails = bf_cnt + bfts_cnt + bstf_cnt
        if total_fails >= num_samples:
            break

    return bs_cnt, bf_cnt, bfts_cnt, bstf_cnt, dem


def collect_samples_with_errors(
    code_name: str,
    ca_code,
    error_rate: ErrorRate,
    noise_model: str,
    d_rounds: int = None,
    total_shots: int = 1_000_000_000,
    batch_shots: int = 100_000,
    sample_times: int = 1000,
    num_samples: int = 10_000,
    build_decoders: Literal['normal', 'serial_bp', 'force_osd'] = 'normal',
    kind: Literal['both_success', 'both_fail', 'bposd_fail_tess_success', 'bposd_success_tess_fail'] = 'bposd_fail_tess_success',
) -> tuple[list[dict], stim.DetectorErrorModel]:
    """
    Collect cases with errors

    Parameters
    ----------
    code_name: str
    ca_code: CA code object
    error_rate: ErrorRate
    noise_model: str
    d_rounds: int, optional
    total_shots: int, default 1_000_000_000
        Total shot budget, not per-iteration shots
    batch_shots: int, default 100_000
        Batch size per sampler.sample call
    sample_times: int, default 1000
        Maximum number of sampling batches
    num_samples: int, default 10_000
        Number of samples to collect, should be less than total_shots
    kind: str, default 'bposd_fail_tess_success'
    """
    builder = CircuitBuilder(ca_code)
    circuit = builder.get_circuit(noise_model=noise_model, error_rate=error_rate, d_rounds=d_rounds)
    dem = circuit.detector_error_model()

    if build_decoders == 'normal':
        bposd, tess = build_decoder_pair(ca_code, error_rate.p)
    elif build_decoders == 'serial_bp':
        bposd, tess = build_bp_serial(ca_code, error_rate.p)
    elif build_decoders == 'force_osd':
        bposd, tess = build_force_osd(ca_code, error_rate.p)
    else:
        raise ValueError(f'Unsupported build_decoders: {build_decoders}')

    compiled_b = bposd.compile_decoder_for_dem(dem=dem)
    compiled_t = tess.compile_decoder_for_dem(dem=dem)

    sampler = dem.compile_sampler()
    num_det = dem.num_detectors
    num_obs = dem.num_observables
    num_err = dem.num_errors

    recs: list[dict] = []
    sampled_shot = 0
    for batch_idx in range(sample_times):
        cur_shots = min(batch_shots, total_shots - sampled_shot)
        if cur_shots <= 0:
            break
        dets, obs, errs = sampler.sample(shots=cur_shots, return_errors=True, bit_packed=True)

        _, unique_idx = np.unique(errs, axis=0, return_index=True)
        if unique_idx.size == 0:
            sampled_shot += cur_shots
            continue

        dets = dets[unique_idx]
        obs = obs[unique_idx]
        errs = errs[unique_idx]

        pred_b = compiled_b.decode_shots_bit_packed(bit_packed_detection_event_data=dets)
        pred_t = compiled_t.decode_shots_bit_packed(bit_packed_detection_event_data=dets)

        b_succ = packed_rows_equal(pred_b, obs, num_obs)
        t_succ = packed_rows_equal(pred_t, obs, num_obs)

        match kind:
            case 'both_success':
                mask = b_succ & t_succ
            case 'both_fail':
                mask = (~b_succ) & (~t_succ)
            case 'bposd_fail_tess_success':
                mask = (~b_succ) & t_succ
            case 'bposd_success_tess_fail':
                mask = b_succ & (~t_succ)
            case _:
                raise ValueError(f'Unsupported kind: {kind}')
            
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            sampled_shot += cur_shots
            continue
        need = num_samples - len(recs)
        if idx.size > need:
            idx = idx[:need]
        
        dets_sel = dets[idx]
        obs_sel = obs[idx]
        errs_sel = errs[idx]
        pred_b_sel = pred_b[idx]
        pred_t_sel = pred_t[idx]

        det_weights = packed_weight(dets_sel, num_det)
        err_weights = packed_weight(errs_sel, num_err)

        recs.extend([
            {
                'code': code_name,
                'p': error_rate.p,
                'p_meas': error_rate.p_meas,
                'p_cx': error_rate.p_cx,
                'p_prep': error_rate.p_prep,
                'num_detectors': int(num_det),
                'num_observables': int(num_obs),
                'num_errors': int(num_err),
                'detector_weight': int(dw),
                'detector_events': packed_to_str(d, num_det),
                'error_weight': int(ew),
                'error_llr_weight': error_llr_weight(dem, packed_to_str(e, num_err)),
                'error_events': packed_to_str(e, num_err),
                'true_observables': packed_to_str(o, num_obs),
                'bposd_pred_observables': packed_to_str(pb, num_obs),
                'tesseract_pred_observables': packed_to_str(pt, num_obs),
            }
            for dw, ew, d, o, e, pb, pt in zip(det_weights, err_weights, dets_sel, obs_sel, errs_sel, pred_b_sel, pred_t_sel)
        ])

        sampled_shot += cur_shots
        
        if len(recs) >= num_samples:
            break

    return recs, dem

def sample_on_error_vec(
    num_errors: int,
    weight: int,
    H: np.ndarray,
    L: np.ndarray,
):
    """
    Sample error vectors of given length and weight.

    Parameters
    ----------
    length: int
        Length of the error vectors
    weight: int
        Hamming weight of the error vectors to sample
    H: np.ndarray
        Parity check matrix of the code, shape (num_checks, num_errors)
    L: np.ndarray
        Logical operator matrix of the code, shape (num_observables, num_errors)

    Returns
    ----------
    errs: np.ndarray
        Array of sampled error vectors, shape (num_cases, num_errors)
    dets: np.ndarray
        Array of corresponding detector events, shape (num_cases, num_detectors)
    obs: np.ndarray
        Array of corresponding observable events, shape (num_cases, num_observables)
    """
    errs = []
    for comb in combinations(range(num_errors), weight):
        err = np.zeros(num_errors, dtype=np.uint8)
        err[list(comb)] = 1
        errs.append(err)

    errs = np.stack(errs)
    dets = (errs @ H.T) % 2
    obs = (errs @ L.T) % 2
    
    return errs, dets, obs
