# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t
import builtins as py_builtins
import numpy as np


cdef vector[vector[uint8_t]] _binary_matrix(object matrix):
    cdef object dense = np.asarray(matrix, dtype=np.uint8)
    if dense.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    dense = np.ascontiguousarray(dense % 2)
    cdef vector[vector[uint8_t]] result
    cdef vector[uint8_t] row
    cdef Py_ssize_t i, j
    for i in range(dense.shape[0]):
        row.clear()
        for j in range(dense.shape[1]):
            row.push_back(<uint8_t>dense[i, j])
        result.push_back(row)
    return result


cdef vector[double] _double_vector(object values):
    cdef vector[double] result
    cdef object array = np.asarray(values, dtype=np.float64).reshape(-1)
    cdef Py_ssize_t i
    for i in range(array.size):
        result.push_back(<double>array[i])
    return result


cdef vector[int] _int_vector(object values):
    cdef vector[int] result
    cdef object array = np.asarray(values, dtype=np.int32).reshape(-1)
    cdef Py_ssize_t i
    for i in range(array.size):
        result.push_back(<int>array[i])
    return result


cdef vector[uint8_t] _bit_vector(object values):
    cdef vector[uint8_t] result
    cdef object array = np.asarray(values, dtype=np.uint8).reshape(-1) % 2
    cdef Py_ssize_t i
    for i in range(array.size):
        result.push_back(<uint8_t>array[i])
    return result


cdef class CorrelationAdaptiveRegionRefiner:
    def __cinit__(
        self,
        pcm,
        logicals,
        priors,
        check_times,
        check_groups,
        check_layers,
        max_region_states=4096,
        max_probe_regions=32,
        region_mode=1,
        score_mode=0,
        region_seed=20260903,
        refine_zero_syndrome=False,
    ):
        cdef vector[vector[uint8_t]] h = _binary_matrix(pcm)
        cdef vector[vector[uint8_t]] l = _binary_matrix(logicals)
        cdef vector[double] p = _double_vector(priors)
        cdef vector[int] times = _int_vector(check_times)
        cdef vector[int] groups = _int_vector(check_groups)
        cdef vector[int] layers = _int_vector(check_layers)
        if isinstance(max_region_states, (py_builtins.bool, np.bool_)) or not isinstance(
            max_region_states, (int, np.integer)
        ):
            raise ValueError("max_region_states must be an integer")
        if max_region_states < 2:
            raise ValueError("max_region_states must be at least two")
        if isinstance(max_probe_regions, (py_builtins.bool, np.bool_)) or not isinstance(
            max_probe_regions, (int, np.integer)
        ):
            raise ValueError("max_probe_regions must be an integer")
        if max_probe_regions < 1 or max_probe_regions > 2147483647:
            raise ValueError("max_probe_regions must fit a positive signed 32-bit integer")
        if isinstance(region_mode, (py_builtins.bool, np.bool_)) or not isinstance(
            region_mode, (int, np.integer)
        ):
            raise ValueError("region_mode must be an integer")
        if region_mode not in (0, 1, 2):
            raise ValueError("region_mode must be 0, 1, or 2")
        if isinstance(score_mode, (py_builtins.bool, np.bool_)) or not isinstance(
            score_mode, (int, np.integer)
        ):
            raise ValueError("score_mode must be an integer")
        if score_mode not in (0, 1, 2, 3, 4):
            raise ValueError("score_mode must be 0, 1, 2, 3, or 4")
        if isinstance(region_seed, (py_builtins.bool, np.bool_)) or not isinstance(
            region_seed, (int, np.integer)
        ):
            raise ValueError("region_seed must be an integer")
        if region_seed < -2147483648 or region_seed > 2147483647:
            raise ValueError("region_seed must fit a signed 32-bit integer")
        if not isinstance(refine_zero_syndrome, (py_builtins.bool, np.bool_)):
            raise ValueError("refine_zero_syndrome must be boolean")
        self._refiner = new CorrelationAdaptiveRegionRefinerCpp(
            h, l, p, times, groups, layers,
            <uint64_t>max_region_states,
            <int>max_probe_regions,
            <int>region_mode,
            <int>score_mode,
            <int>region_seed,
            <bool>refine_zero_syndrome,
        )

    def __dealloc__(self):
        del self._refiner

    def refine(self, syndrome, main_decoding, posterior_llr):
        cdef vector[uint8_t] s = _bit_vector(syndrome)
        cdef vector[uint8_t] main = _bit_vector(main_decoding)
        cdef vector[double] llr = _double_vector(posterior_llr)
        cdef vector[uint8_t] result = self._refiner.refine(s, main, llr)
        cdef Py_ssize_t i
        output = np.empty(result.size(), dtype=np.uint8)
        for i in range(result.size()):
            output[i] = result[i]
        return output

    @property
    def stats(self):
        return {
            "structural_atom_count": self._refiner.structural_atom_count,
            "eligible_atom_count": self._refiner.eligible_atom_count,
            "probe_region_count": self._refiner.probe_region_count,
            "candidate_evaluation_count": self._refiner.candidate_evaluation_count,
            "selected_atom_count": self._refiner.selected_atom_count,
            "selected_region_variables": list(self._refiner.selected_region_variables),
            "selected_region_size": self._refiner.selected_region_variables.size(),
            "probe_enumerated_state_count": self._refiner.probe_enumerated_state_count,
            "selected_region_state_count": self._refiner.selected_region_state_count,
            "selected_feasible_state_count": self._refiner.selected_feasible_state_count,
            "total_correlation_bits": self._refiner.total_correlation_bits,
            "marginal_kl_bits": self._refiner.marginal_kl_bits,
            "factorization_gap_bits": self._refiner.factorization_gap_bits,
            "selected_utility": self._refiner.selected_utility,
            "selected_logical_label": self._refiner.selected_logical_label,
            "logical_log_masses": list(self._refiner.logical_log_masses),
            "main_syndrome_valid": True if self._refiner.main_syndrome_valid else False,
            "used_feasible_seed": True if self._refiner.used_feasible_seed else False,
            "refined": True if self._refiner.refined else False,
            "syndrome_valid": True if self._refiner.syndrome_valid else False,
        }
