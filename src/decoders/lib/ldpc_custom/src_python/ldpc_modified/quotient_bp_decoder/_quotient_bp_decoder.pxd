# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++

from libcpp cimport bool
from libcpp.vector cimport vector
cimport numpy as np

ctypedef np.uint8_t uint8_t


cdef extern from "quotient_bp.hpp" namespace "ldpc_modified::quotient_bp":
    cdef const vector[int] NULL_INT_VECTOR

    cdef enum BpMethod:
        PRODUCT_SUM = 0
        MINIMUM_SUM = 1

    cdef enum BpInputType:
        SYNDROME = 0
        RECEIVED_VECTOR = 1
        AUTO = 2

    cdef enum BpSchedule:
        SERIAL = 0
        PARALLEL = 1
        SERIAL_RELATIVE = 2

    cdef cppclass BpEntry "ldpc_modified::quotient_bp::BpEntry":
        BpEntry() except +

    cdef cppclass BpSparse "ldpc_modified::quotient_bp::BpSparse":
        int m
        int n
        BpSparse() except +
        BpSparse(int m, int n, int entry_count) except +
        BpEntry& insert_entry(int row, int column)
        vector[vector[int]] nonzero_coordinates()
        int entry_count()
        vector[uint8_t] mulvec(vector[uint8_t]& input_vector)

    cdef cppclass BpDecoderCpp "ldpc_modified::quotient_bp::BpDecoder":
        BpDecoderCpp(
            BpSparse& parity_check_matrix,
            vector[double] channel_probabilities,
            int maximum_iterations,
            BpMethod bp_method,
            BpSchedule schedule,
            double min_sum_scaling_factor,
            int omp_threads,
            vector[int] serial_schedule,
            int random_schedule_seed,
            bool random_serial_schedule,
            BpInputType bp_input_type,
        ) except +
        vector[double] channel_probabilities
        int check_count
        int bit_count
        int maximum_iterations
        BpMethod bp_method
        BpSchedule schedule
        double ms_scaling_factor
        vector[uint8_t] decoding
        vector[double] log_prob_ratios
        vector[int] serial_schedule_order
        int iterations
        int omp_thread_count
        bool converge
        int random_schedule_seed
        bool random_serial_schedule
        vector[vector[uint8_t]] trajectory_logical_cosets
        vector[vector[uint8_t]] logical_cosets
        vector[uint8_t] decode(vector[uint8_t]& syndrome) except +
        void enable_logical_tracking(BpSparse& logical, bool collect_full) except +
        void set_random_schedule_seed(int seed)

    cdef cppclass QuotientModelCpp "ldpc_modified::quotient_bp::QuotientModel":
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
        BpSparse& d_matrix()
        vector[uint8_t] particular_solution(vector[uint8_t] syndrome) except +
        vector[uint8_t] representative(
            vector[uint8_t] syndrome,
            vector[uint8_t] observable_label,
        ) except +

    cdef cppclass DPartitionBpDecoderCpp "ldpc_modified::quotient_bp::DPartitionBpDecoder":
        DPartitionBpDecoderCpp(
            BpSparse& d,
            vector[double] probabilities,
            int maximum_iterations,
            double damping,
            double tolerance,
        ) except +
        vector[double] channel_probabilities
        int maximum_iterations
        double damping
        double tolerance
        int iterations
        bool converge
        double max_message_change
        double log_partition
        void set_channel_probabilities(const vector[double]& probabilities) except +
        double decode(vector[uint8_t]& offset) except +


cdef class QuotientBpDecoder:
    cdef BpSparse* h_pcm
    cdef BpSparse* logical_pcm
    cdef QuotientModelCpp* model
    cdef BpDecoderCpp* h_bpd
    cdef DPartitionBpDecoderCpp* d_bpd
    cdef vector[uint8_t] _syndrome
    cdef vector[double] _error_channel
    cdef vector[int] _serial_schedule_order
    cdef vector[uint8_t] _decoding
    cdef vector[uint8_t] _chosen_logical_coset
    cdef vector[vector[uint8_t]] _candidate_representatives
    cdef vector[double] _candidate_log_partitions
    cdef vector[int] _d_iterations
    cdef vector[uint8_t] _d_converged
    cdef int m
    cdef int n
    cdef int k
