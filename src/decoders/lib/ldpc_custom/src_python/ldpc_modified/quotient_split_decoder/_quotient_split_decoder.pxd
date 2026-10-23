# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++

from libcpp cimport bool
from libcpp.vector cimport vector
cimport numpy as np

ctypedef np.uint8_t uint8_t


cdef extern from "quotient_split.hpp" namespace "ldpc_modified::quotient_split":
    cdef cppclass BpEntry "ldpc_modified::quotient_split::BpEntry":
        BpEntry() except +

    cdef cppclass BpSparse "ldpc_modified::quotient_split::BpSparse":
        int m
        int n
        BpSparse() except +
        BpSparse(int m, int n, int entry_count) except +
        BpEntry& insert_entry(int row, int column)

    cdef cppclass QuotientModelCpp "ldpc_modified::quotient_split::QuotientModel":
        QuotientModelCpp(BpSparse& h, BpSparse& logical) except +
        int check_count
        int bit_count
        int logical_count
        int rank_h
        int logical_rank
        int kernel_dimension
        int degeneracy_dimension
        vector[vector[int]] T
        vector[vector[int]] K
        vector[vector[int]] Q
        vector[vector[int]] D
        vector[vector[int]] logical_map

    cdef cppclass QuotientSplitLeafResultCpp "ldpc_modified::quotient_split::QuotientSplitLeafResult":
        vector[uint8_t] logical_label
        vector[uint8_t] correction
        double path_score
        double log_partition
        bool bp_converged
        int bp_iterations

    cdef cppclass QuotientSplitDecodeResultCpp "ldpc_modified::quotient_split::QuotientSplitDecodeResult":
        vector[uint8_t] best_logical_label
        vector[uint8_t] correction
        vector[QuotientSplitLeafResultCpp] leaves
        vector[int] surviving_prefix_counts

    cdef cppclass QuotientSplitDecoderCpp "ldpc_modified::quotient_split::QuotientSplitDecoder":
        QuotientSplitDecoderCpp(
            QuotientModelCpp& model,
            vector[double] probabilities,
            int list_width,
            vector[int] split_order,
            bool adaptive_split,
            int internal_iterations,
            int leaf_iterations,
            double damping,
            double tolerance,
            double probability_floor,
        ) except +
        vector[double] channel_probabilities
        int list_width
        vector[int] split_order
        bool adaptive_split
        int internal_iterations
        int leaf_iterations
        double damping
        double tolerance
        double probability_floor
        void set_channel_probabilities(
            const vector[double]& probabilities,
        ) except +
        QuotientSplitDecodeResultCpp decode(
            vector[uint8_t] syndrome,
        ) except +


cdef class QuotientSplitDecoder:
    cdef BpSparse* h_pcm
    cdef BpSparse* logical_pcm
    cdef QuotientModelCpp* model
    cdef QuotientSplitDecoderCpp* split_decoder
    cdef vector[uint8_t] _syndrome
    cdef vector[double] _error_channel
    cdef QuotientSplitDecodeResultCpp _result
    cdef int m
    cdef int n
    cdef int k
