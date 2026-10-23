#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
from libc.stdlib cimport malloc, calloc, free
from libcpp cimport bool
from libcpp.vector cimport vector
cimport numpy as np
ctypedef np.uint8_t uint8_t

cdef extern from "expanded_adaptive_bp.hpp" namespace "ldpc_modified::expanded_adaptive_bp":

    cdef const vector[int] NULL_INT_VECTOR
    int expanded_adaptive_bp_live_instance_count()

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

    cdef cppclass BpEntry "ldpc_modified::expanded_adaptive_bp::BpEntry":
        BpEntry() except +
        bool at_end()

    cdef cppclass BpSparse "ldpc_modified::expanded_adaptive_bp::BpSparse":
        int m
        int n
        BpSparse() except +
        BpSparse(int m, int n, int entry_count) except +
        BpEntry& insert_entry(int i, int j)
        BpEntry& get_entry(int i, int j)
        vector[uint8_t]& mulvec(vector[uint8_t]& input_vector, vector[uint8_t]& output_vector)
        vector[uint8_t] mulvec(vector[uint8_t]& input_vector)

        vector[vector[int]] nonzero_coordinates()
        int entry_count()

        int get_col_degree(int col)
        int get_row_degree(int row)

    cdef cppclass ExpandedAdaptiveBpDecoderCpp "ldpc_modified::expanded_adaptive_bp::ExpandedAdaptiveBpDecoder":
            ExpandedAdaptiveBpDecoderCpp(
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
                BpInputType bp_input_type) except +
            BpSparse& pcm
            vector[double] channel_probabilities
            int check_count
            int bit_count
            int maximum_iterations
            BpMethod bp_method
            BpSchedule schedule
            double ms_scaling_factor
            vector[uint8_t] decoding
            vector[uint8_t] candidate_syndrome
            vector[double] log_prob_ratios
            vector[double] initial_log_prob_ratios
            vector[int] serial_schedule_order
            vector[double] llr_max
            vector[double] llr_min
            vector[int] sign_flip_count
            vector[double] instability_scores
            vector[double] instability_flip_component
            vector[double] instability_movement_component
            vector[double] instability_message_conflict_component
            vector[double] instability_residual_component
            vector[double] instability_weak_belief_component
            vector[double] instability_spatial_component
            vector[double] instability_seed_scores
            vector[int] adaptive_region_variables
            vector[uint8_t] adaptive_final_decoding
            int adaptive_executed_view_count
            int adaptive_main_bp_run_count
            int adaptive_core_variable_count
            int adaptive_unstable_variable_count
            int adaptive_feasible_unstable_variable_count
            int adaptive_covered_unstable_variable_count
            int adaptive_candidate_count
            int adaptive_valid_candidate_count
            int adaptive_selected_view
            int adaptive_fallback_residual_weight
            double adaptive_best_candidate_cost
            bool adaptive_selected_candidate_valid
            int iterations
            int omp_thread_count
            bool converge
            int random_schedule_seed
            bool random_serial_schedule
            vector[uint8_t] decode(vector[uint8_t]& syndrome)
            void prepare_adaptive_redundant_views(vector[uint8_t]& syndrome, int adaptive_max_views) except +
            void build_prepared_adaptive_redundant_view(int view_index, vector[uint8_t]& syndrome, BpSparse& augmented_pcm, vector[uint8_t]& augmented_syndrome) except +
            void ensure_instability_scores()
            int adaptive_prepared_view_count()
            int adaptive_candidate_pair_count()
            int adaptive_selected_row_count(int view_index)
            int adaptive_weight_limit(int view_index)
            int adaptive_selected_row_weight(int view_index, int row_index)
            void reset_adaptive_multi_view_state()
            void consider_adaptive_candidate(vector[uint8_t]& syndrome, const vector[uint8_t]& candidate, const vector[uint8_t]& candidate_check_values, int source_view) except +
            void update_decoder_parameters(const vector[double]& channel_probabilities, int maximum_iterations, BpMethod bp_method, BpSchedule schedule, double min_sum_scaling_factor, int omp_threads, vector[int]& serial_schedule, int random_schedule_seed, bool random_serial_schedule, BpInputType bp_input_type) except +
            void refresh_initial_log_prob_ratios()
            void reset_zero_syndrome_state()
            void set_omp_thread_count(int count)
            BpInputType bp_input_type
            void set_random_schedule_seed(int seed)

    cdef cppclass AdaptiveViewDecoderPool "ldpc_modified::expanded_adaptive_bp::AdaptiveViewDecoderPool":
            AdaptiveViewDecoderPool() except +
            void run(ExpandedAdaptiveBpDecoderCpp& main_decoder, vector[uint8_t]& syndrome, int view_thread_count) except +
            void reset()
            int active_view_count()
            int last_worker_count()
            ExpandedAdaptiveBpDecoderCpp& decoder(int view_index) except +
            BpSparse& pcm(int view_index) except +

cdef class ExpandedAdaptiveBpDecoderBase:
    cdef BpSparse *pcm
    cdef int m, n
    cdef vector[uint8_t] _syndrome
    cdef vector[double] _error_channel
    cdef vector[int] _serial_schedule_order
    cdef bool MEMORY_ALLOCATED
    cdef ExpandedAdaptiveBpDecoderCpp *bpd
    cdef AdaptiveViewDecoderPool *view_pool
    cdef list _view_candidates
    cdef list _view_iterations
    cdef list _view_edge_counts
    cdef list _view_validity
    cdef list _view_convergence
    cdef void _reset_view_telemetry(self)
    cdef void _record_view_telemetry(self, ExpandedAdaptiveBpDecoderCpp *decoder, BpSparse *view_pcm) except *
    cdef str user_dtype
    # cdef int random_schedule_seed
    cdef int _adaptive_max_views
    cdef int _view_thread_count

cdef class ExpandedAdaptiveBpDecoder(ExpandedAdaptiveBpDecoderBase):
    cdef vector[uint8_t] _received_vector
    pass
