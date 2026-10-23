#cython: language_level=3
from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
cimport numpy as np

ctypedef np.uint8_t uint8_t

cdef extern from "mbbp_ld.hpp" namespace "ldpc_modified::mbbp_ld":
    cdef cppclass BpEntry:
        bool at_end()

    cdef cppclass BpSparse:
        int m
        int n
        BpSparse(int m, int n, int entry_count) except +
        BpEntry& insert_entry(int row, int column)

    cdef cppclass MbbpLdDecoderCpp "ldpc_modified::mbbp_ld::MbbpLdDecoder":
        MbbpLdDecoderCpp(
            BpSparse& parity_check_matrix,
            const vector[double]& channel_probabilities,
            int maximum_iterations,
            int bp_method,
            int schedule,
            double min_sum_scaling_factor,
            const vector[vector[int]]& check_permutations) except +
        bool converge
        vector[uint8_t] decoding
        int candidate_count
        int unique_candidate_count
        int selected_frequency
        int selected_weight
        vector[uint8_t] decode(vector[uint8_t]& syndrome) except +
        int partition_count()
        int partition_subtree_count(int partition_index)
        int partition_subtree_size(int partition_index, int subtree_index)
        int partition_subtree_row(int partition_index, int subtree_index, int row_index)
        int basis_count()
        int basis_row_count(int basis_index)
        int basis_source_row(int basis_index, int row_index)
        int basis_candidate_bit(int basis_index, int bit_index)
        bool basis_converged_at(int basis_index)
        int basis_iterations_at(int basis_index)
        int basis_edge_count_at(int basis_index)
        int converged_candidate_bit(int candidate_index, int bit_index)
        int unique_candidate_frequency_at(int candidate_index)
        int unique_candidate_weight_at(int candidate_index)
        int unique_candidate_bit(int candidate_index, int bit_index)
        double channel_probability(int bit_index)

cdef class MbbpLdDecoder:
    cdef BpSparse *pcm
    cdef MbbpLdDecoderCpp *decoder
    cdef vector[uint8_t] _syndrome
    cdef int m
    cdef int n
