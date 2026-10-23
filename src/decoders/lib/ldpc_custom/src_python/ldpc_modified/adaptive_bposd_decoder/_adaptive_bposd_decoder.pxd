#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
from libc.stdlib cimport malloc, calloc, free
from libcpp cimport bool
from libcpp.vector cimport vector
cimport numpy as np
from ldpc_modified.adaptive_bposd_decoder._adaptive_bp_decoder cimport BpSparse, AdaptiveBpDecoderBase, AdaptiveBpDecoderCpp, NULL_INT_VECTOR, SYNDROME
ctypedef np.uint8_t uint8_t

cdef extern from "adaptive_osd.hpp" namespace "ldpc_modified::adaptive_osd" nogil:

    cdef enum OsdMethod:
        OSD_OFF = 0
        OSD_0 = 1
        EXHAUSTIVE = 2
        COMBINATION_SWEEP = 3

    cdef cppclass AdaptiveOsdDecoderCpp "ldpc_modified::adaptive_osd::AdaptiveOsdDecoder":
        AdaptiveOsdDecoderCpp(BpSparse& parity_check_matrix, OsdMethod osd_method, int osd_order, vector[double]& channel_probabilities) except +
        vector[uint8_t]& decode(vector[uint8_t]& syndrome, vector[double]& log_prob_ratios)
        vector[uint8_t]& adaptive_decode_from_llrs(vector[uint8_t]& syndrome, vector[double]& l0, vector[double]& l1, vector[uint8_t]& bp_l0_decoding, vector[uint8_t]& bp_l1_decoding, int auxiliary_candidate_pair_count, int auxiliary_selected_row_count, bool auxiliary_bp_converge) except +
        void reset_zero_syndrome_stats()
        vector[uint8_t] osd0_decoding
        vector[uint8_t] osdw_decoding
        vector[uint8_t] osd0_l0_decoding
        vector[uint8_t] osdw_l0_decoding
        vector[uint8_t] osd0_l1_decoding
        vector[uint8_t] osdw_l1_decoding
        vector[uint8_t] final_decoding
        OsdMethod osd_method
        int osd_setup()
        int osd_order
        int auxiliary_candidate_pair_count
        int auxiliary_selected_row_count
        int auxiliary_view_disagreement_weight
        int merged_candidate_count
        double best_candidate_cost
        bool auxiliary_view_disagreement
        bool auxiliary_bp_converge

cdef class AdaptiveBpOsdDecoder(AdaptiveBpDecoderBase):
    cdef AdaptiveOsdDecoderCpp* osdD
    cdef public bint force_osd
