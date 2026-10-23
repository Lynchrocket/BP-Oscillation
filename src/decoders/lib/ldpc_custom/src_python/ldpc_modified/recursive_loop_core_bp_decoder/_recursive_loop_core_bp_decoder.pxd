from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t

cdef extern from "recursive_loop_core_bp.hpp" namespace "ldpc_modified::recursive_loop_core_bp":
    cdef cppclass RecursiveLoopCoreRefinerCpp "ldpc_modified::recursive_loop_core_bp::RecursiveLoopCoreRefiner":
        RecursiveLoopCoreRefinerCpp(
            const vector[vector[uint8_t]]& pcm,
            const vector[vector[uint8_t]]& logicals,
            const vector[double]& priors,
            const vector[int]& check_times,
            const vector[int]& check_groups,
            const vector[int]& check_layers,
            uint64_t max_factor_states,
            int cluster_mode,
            bool sum_logical_mass,
        ) except +
        vector[uint8_t]& refine(
            const vector[uint8_t]& syndrome,
            const vector[uint8_t]& main_decoding,
            const vector[double]& posterior_llr,
        ) except +
        vector[uint8_t] decoding
        vector[int] active_variables
        vector[vector[int]] fundamental_circuits
        vector[double] logical_log_masses
        int tree_variable_count
        int loop_core_variable_count
        int loop_component_count
        int exact_component_count
        int truncated_component_count
        uint64_t enumerated_state_count
        uint64_t feasible_state_count
        bool refined
        bool syndrome_valid
        int selected_logical_label

cdef class RecursiveLoopCoreRefiner:
    cdef RecursiveLoopCoreRefinerCpp* _refiner
