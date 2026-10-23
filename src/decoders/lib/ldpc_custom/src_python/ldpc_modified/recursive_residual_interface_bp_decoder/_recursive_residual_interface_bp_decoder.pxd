from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t

cdef extern from "recursive_residual_interface_bp.hpp" namespace "ldpc_modified::recursive_residual_interface_bp":
    cdef cppclass RecursiveResidualInterfaceRefinerCpp "ldpc_modified::recursive_residual_interface_bp::RecursiveResidualInterfaceRefiner":
        RecursiveResidualInterfaceRefinerCpp(
            const vector[vector[uint8_t]]& pcm,
            const vector[vector[uint8_t]]& logicals,
            const vector[double]& priors,
            const vector[int]& variable_order,
            uint64_t max_interface_states,
            bool residual_compression,
        ) except +
        vector[uint8_t]& refine(
            const vector[uint8_t]& syndrome,
            const vector[uint8_t]& main_decoding,
            const vector[double]& posterior_llr,
        ) except +
        vector[uint8_t] decoding
        vector[double] logical_log_masses
        uint64_t expanded_pair_count
        uint64_t merged_path_count
        uint64_t pruned_state_count
        uint64_t closed_check_rejection_count
        uint64_t terminal_state_count
        uint64_t peak_state_count
        int region_count
        int compressed_region_count
        int peak_boundary_check_count
        int selected_logical_label
        double residual_spectral_energy
        double mean_effective_walsh_rank
        bool truncated
        bool used_main_fallback
        bool refined
        bool syndrome_valid

cdef class RecursiveResidualInterfaceRefiner:
    cdef RecursiveResidualInterfaceRefinerCpp* _refiner
