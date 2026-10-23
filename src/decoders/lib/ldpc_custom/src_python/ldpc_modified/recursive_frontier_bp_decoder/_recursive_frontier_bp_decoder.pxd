from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t

cdef extern from "recursive_frontier_bp.hpp" namespace "ldpc_modified::recursive_frontier_bp":
    cdef cppclass RecursiveFrontierRefinerCpp "ldpc_modified::recursive_frontier_bp::RecursiveFrontierRefiner":
        RecursiveFrontierRefinerCpp(
            const vector[vector[uint8_t]]& pcm,
            const vector[vector[uint8_t]]& logicals,
            const vector[double]& priors,
            const vector[int]& variable_order,
            const vector[int]& bag_ends,
            uint64_t max_frontier_states,
        ) except +
        vector[uint8_t]& refine(
            const vector[uint8_t]& syndrome,
            const vector[uint8_t]& main_decoding,
            const vector[double]& posterior_llr,
        ) except +
        vector[uint8_t] decoding
        vector[double] logical_log_masses
        uint64_t expanded_transition_count
        uint64_t merged_path_count
        uint64_t pruned_state_count
        uint64_t closed_check_rejection_count
        uint64_t terminal_state_count
        uint64_t peak_state_count
        int peak_frontier_check_count
        int recursive_bag_count
        int peak_recursive_bag_size
        int selected_logical_label
        bool truncated
        bool used_main_fallback
        bool refined
        bool syndrome_valid

cdef class RecursiveFrontierRefiner:
    cdef RecursiveFrontierRefinerCpp* _refiner
