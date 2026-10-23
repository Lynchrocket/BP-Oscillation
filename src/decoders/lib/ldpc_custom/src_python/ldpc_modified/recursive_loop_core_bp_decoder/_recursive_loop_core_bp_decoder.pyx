# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t
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


cdef class RecursiveLoopCoreRefiner:
    def __cinit__(
        self,
        pcm,
        logicals,
        priors,
        check_times,
        check_groups,
        check_layers,
        max_factor_states=65536,
        cluster_mode=2,
        sum_logical_mass=True,
    ):
        cdef vector[vector[uint8_t]] h = _binary_matrix(pcm)
        cdef vector[vector[uint8_t]] l = _binary_matrix(logicals)
        cdef vector[double] p = _double_vector(priors)
        cdef vector[int] times = _int_vector(check_times)
        cdef vector[int] groups = _int_vector(check_groups)
        cdef vector[int] layers = _int_vector(check_layers)
        if max_factor_states < 1:
            raise ValueError("max_factor_states must be positive")
        if cluster_mode not in (0, 1, 2):
            raise ValueError("cluster_mode must be 0, 1, or 2")
        self._refiner = new RecursiveLoopCoreRefinerCpp(
            h, l, p, times, groups, layers,
            <uint64_t>max_factor_states, <int>cluster_mode, <bool>sum_logical_mass,
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
            "tree_variable_count": self._refiner.tree_variable_count,
            "loop_core_variable_count": self._refiner.loop_core_variable_count,
            "loop_component_count": self._refiner.loop_component_count,
            "exact_component_count": self._refiner.exact_component_count,
            "truncated_component_count": self._refiner.truncated_component_count,
            "active_variable_count": self._refiner.active_variables.size(),
            "active_variables": list(self._refiner.active_variables),
            "fundamental_circuits": [list(row) for row in self._refiner.fundamental_circuits],
            "enumerated_state_count": self._refiner.enumerated_state_count,
            "feasible_state_count": self._refiner.feasible_state_count,
            "refined": True if self._refiner.refined else False,
            "syndrome_valid": True if self._refiner.syndrome_valid else False,
            "selected_logical_label": self._refiner.selected_logical_label,
            "logical_log_masses": list(self._refiner.logical_log_masses),
        }
