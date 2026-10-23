from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t

cdef extern from "correlation_adaptive_bp.hpp" namespace "ldpc_modified::correlation_adaptive_bp":
    cdef cppclass CorrelationAdaptiveRegionRefinerCpp "ldpc_modified::correlation_adaptive_bp::CorrelationAdaptiveRegionRefiner":
        CorrelationAdaptiveRegionRefinerCpp(
            const vector[vector[uint8_t]]& pcm,
            const vector[vector[uint8_t]]& logicals,
            const vector[double]& priors,
            const vector[int]& check_times,
            const vector[int]& check_groups,
            const vector[int]& check_layers,
            uint64_t max_region_states,
            int max_probe_regions,
            int region_mode,
            int score_mode,
            int region_seed,
            bool refine_zero_syndrome,
        ) except +
        vector[uint8_t]& refine(
            const vector[uint8_t]& syndrome,
            const vector[uint8_t]& main_decoding,
            const vector[double]& posterior_llr,
        ) except +
        vector[uint8_t] decoding
        vector[int] selected_region_variables
        vector[double] logical_log_masses
        int structural_atom_count
        int eligible_atom_count
        int probe_region_count
        int candidate_evaluation_count
        int selected_atom_count
        uint64_t probe_enumerated_state_count
        uint64_t selected_region_state_count
        uint64_t selected_feasible_state_count
        double total_correlation_bits
        double marginal_kl_bits
        double factorization_gap_bits
        double selected_utility
        bool main_syndrome_valid
        bool used_feasible_seed
        bool refined
        bool syndrome_valid
        int selected_logical_label

cdef class CorrelationAdaptiveRegionRefiner:
    cdef CorrelationAdaptiveRegionRefinerCpp* _refiner
