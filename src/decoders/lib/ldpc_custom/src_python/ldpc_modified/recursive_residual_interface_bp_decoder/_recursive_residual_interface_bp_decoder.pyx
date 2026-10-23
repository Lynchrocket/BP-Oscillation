# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t
import numpy as np


cdef vector[vector[uint8_t]] _binary_matrix(object matrix):
    cdef object dense = np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) % 2)
    if dense.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
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


cdef class RecursiveResidualInterfaceRefiner:
    def __cinit__(self, pcm, logicals, priors, variable_order,
                  max_interface_states=512, residual_compression=True):
        cdef vector[vector[uint8_t]] h = _binary_matrix(pcm)
        cdef vector[vector[uint8_t]] l = _binary_matrix(logicals)
        cdef vector[double] p = _double_vector(priors)
        cdef vector[int] order = _int_vector(variable_order)
        if max_interface_states < 1:
            raise ValueError("max_interface_states must be positive")
        self._refiner = new RecursiveResidualInterfaceRefinerCpp(
            h, l, p, order, <uint64_t>max_interface_states, residual_compression,
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
            "region_count": self._refiner.region_count,
            "compressed_region_count": self._refiner.compressed_region_count,
            "peak_boundary_check_count": self._refiner.peak_boundary_check_count,
            "peak_state_count": self._refiner.peak_state_count,
            "expanded_pair_count": self._refiner.expanded_pair_count,
            "merged_path_count": self._refiner.merged_path_count,
            "pruned_state_count": self._refiner.pruned_state_count,
            "closed_check_rejection_count": self._refiner.closed_check_rejection_count,
            "terminal_state_count": self._refiner.terminal_state_count,
            "selected_logical_label": self._refiner.selected_logical_label,
            "logical_log_masses": list(self._refiner.logical_log_masses),
            "residual_spectral_energy": self._refiner.residual_spectral_energy,
            "mean_effective_walsh_rank": self._refiner.mean_effective_walsh_rank,
            "truncated": True if self._refiner.truncated else False,
            "used_main_fallback": True if self._refiner.used_main_fallback else False,
            "refined": True if self._refiner.refined else False,
            "syndrome_valid": True if self._refiner.syndrome_valid else False,
        }
