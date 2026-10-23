# bp_failure_dynamics, failure_regions, force_osd_failures, gap_mechanisms, instability_guided_bp, low_weight_failures, osd, quotient, trapping_sets
from pathlib import Path

import pandas as pd
import numpy as np

def id2str(id_list, length=None):
    """
    [2, 4, 7] => '00101001'
    """
    if length is None:
        length = max(id_list) + 1
    s = ['0'] * length
    for i in id_list:
        s[i] = '1'
    return ''.join(s)

def str2id(string, offset=0):
    """
    '00101001' => [2, 4, 7]
    """
    return [i + offset for i, ch in enumerate(string) if ch == '1']

def id2bits(id_list, length=None):
    """
    [2, 4, 7] => [0, 0, 1, 0, 1, 0, 0, 1]
    """
    if length is None:
        length = max(id_list) + 1
    bits = np.zeros(length, dtype=np.uint8)
    bits[list(id_list)] = 1
    return bits

def bits2id(bits, offset=0):
    """
    [0, 0, 1, 0, 1, 0, 0, 1] => [2, 4, 7]
    """
    return [i + offset for i in range(len(bits)) if bits[i] == 1]

def str2bits(string):
    """
    '00101001' => [0, 0, 1, 0, 1, 0, 0, 1]
    """
    return [int(ch) for ch in string]

def bits2str(x: np.ndarray) -> str:
    """
    [0, 0, 1, 0, 1, 0, 0, 1] => '00101001'
    ord('0')=48
    """
    # return ''.join(str(v) for v in x.tolist())
    return (x + 48).tobytes().decode('ascii')

def xor_strs(a: str, b: str) -> str:
    """
    XOR two binary strings of equal length.
    """
    assert len(a) == len(b), "Strings must be of equal length"
    return ''.join('1' if x != y else '0' for x, y in zip(a, b))

def and_strs(a: str, b: str) -> str:
    """
    AND two binary strings of equal length.
    """
    assert len(a) == len(b), "Strings must be of equal length"
    return ''.join('1' if x == '1' and y == '1' else '0' for x, y in zip(a, b))

def or_strs(a: str, b: str) -> str:
    """
    OR two binary strings of equal length.
    """
    assert len(a) == len(b), "Strings must be of equal length"
    return ''.join('1' if x == '1' or y == '1' else '0' for x, y in zip(a, b))

def pack_bits_little(arr2d_uint8: np.ndarray) -> np.ndarray:
    return np.packbits(arr2d_uint8.astype(np.uint8), axis=1, bitorder='little')

def unpack_bits_little(arr2d_packed: np.ndarray, n_bits: int) -> np.ndarray:
    return np.unpackbits(arr2d_packed, axis=1, count=n_bits, bitorder='little').astype(np.uint8)

def p2llr(p: np.ndarray) -> np.ndarray:
    """
    Convert probability to log-likelihood ratio (LLR).
    LLR = log(p / (1 - p))
    """
    eps = 1e-12  # To avoid division by zero
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))

def llr2p(llr: np.ndarray) -> np.ndarray:
    """
    Convert log-likelihood ratio (LLR) to probability.
    p = 1 / (1 + exp(-LLR))
    """
    return 1 / (1 + np.exp(-llr))

def bit_array(text):
    return np.fromiter((ord(x) - 48 for x in text), dtype=np.uint8)

def log(level, payload, *, CONFIG):
    if CONFIG['print_level'] >= level:
        print(payload)

def bit_text(bits):
    return ''.join((str(int(x)) for x in np.asarray(bits).reshape(-1)))

def aggregate_path(path, *, data_paths):
    path = Path(path)
    target = data_paths.stats() / path.relative_to(data_paths.trace())
    target.parent.mkdir(parents=True, exist_ok=True)
    return target

def bit_array_variant_2(text):
    return np.fromiter((ord(value) - 48 for value in str(text)), dtype=np.uint8)

def bit_text_variant_2(bits):
    return ''.join((str(int(value)) for value in np.asarray(bits).reshape(-1)))

def save_csv(frame, path):
    frame.to_csv(path, index=False)
    return frame

def bits(values):
    return ''.join((str(int(value)) for value in np.asarray(values, dtype=np.uint8).reshape(-1)))

def bit_array_variant_3(value, length):
    text = str(value).strip().zfill(length)
    assert len(text) == length and set(text) <= {'0', '1'}
    return np.fromiter((int(char) for char in text), dtype=np.uint8, count=length)

def stage_cached(paths, *, CONFIG):
    present = [path.is_file() for path in paths]
    assert all(present) or not any(present), f'partial cache: {[str(p) for p, ok in zip(paths, present) if ok]}'
    return CONFIG['reuse_existing_data'] and all(present)

def read_csv(path):
    assert path.is_file(), path
    return pd.read_csv(path, dtype={'shot_id': str, 'syndrome': str, 'true_logical': str, 'error': str}, low_memory=False)

def bits_variant_2(values):
    return ''.join((str(int(x)) for x in np.asarray(values, dtype=np.uint8).reshape(-1)))


def check_algorithm_cache(data_paths, config, versions):
    """Reject old algorithm results without changing or deleting their files."""
    import json
    marker = data_paths.output('config') / 'algorithm_versions.json'
    has_results = any(folder.exists() and any(folder.iterdir())
                      for folder in (data_paths.output('stats'), data_paths.output('trace')))
    if has_results:
        assert marker.is_file(), 'Existing results predate the algorithm fixes; select a fresh run path before rerunning'
        assert json.loads(marker.read_text()) == versions, 'Algorithm cache mismatch; select a fresh run path'
    config['algorithm_versions'] = dict(versions)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(versions, indent=2) + '\n')
