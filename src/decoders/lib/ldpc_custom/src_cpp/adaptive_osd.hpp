#ifndef ADAPTIVE_OSD_H
#define ADAPTIVE_OSD_H

#include <vector>
#include <memory>
#include <iterator>
#include <cmath> 
#include <limits>
#include <algorithm>
#include <cassert>

#include "adaptive_bp.hpp"
#include "sort.hpp"
#include "util.hpp"
#include "gf2sparse_linalg.hpp"

namespace ldpc_modified::adaptive_osd{


enum OsdMethod{
    OSD_OFF,
    OSD_0,
    EXHAUSTIVE,
    COMBINATION_SWEEP
};


class AdaptiveOsdDecoder{
    public:
        OsdMethod osd_method;
        int osd_order;
        int k, bit_count, check_count;
        ldpc_modified::adaptive_bp::BpSparse& pcm;
        std::vector<double>& channel_probabilities;
        std::vector<uint8_t> osd0_decoding;
        std::vector<uint8_t> osdw_decoding;
        std::vector<uint8_t> osd0_l0_decoding;
        std::vector<uint8_t> osdw_l0_decoding;
        std::vector<uint8_t> osd0_l1_decoding;
        std::vector<uint8_t> osdw_l1_decoding;
        std::vector<uint8_t> bp_l0_decoding;
        std::vector<uint8_t> bp_l1_decoding;
        std::vector<uint8_t> final_decoding;
        std::vector<std::vector<uint8_t>> osd_candidate_strings;
        std::vector<int> column_ordering;
        ldpc_modified::gf2sparse_linalg::RowReduce<ldpc_modified::adaptive_bp::BpEntry>* LuDecomposition;
        int auxiliary_candidate_pair_count;
        int auxiliary_selected_row_count;
        int auxiliary_view_disagreement_weight;
        int merged_candidate_count;
        double best_candidate_cost;
        bool auxiliary_view_disagreement;
        bool auxiliary_bp_converge;
        
        AdaptiveOsdDecoder(
            ldpc_modified::adaptive_bp::BpSparse& parity_check_matrix,
            OsdMethod osd_method,
            int osd_order,
            std::vector<double>& channel_probs):
            pcm(parity_check_matrix),
            channel_probabilities(channel_probs),
            LuDecomposition(nullptr)
        {

            this->bit_count = this->pcm.n;
            this->check_count = this->pcm.m;

            this->osd_order = osd_order;
            this->osd_method = osd_method;
            this->auxiliary_candidate_pair_count = 0;
            this->auxiliary_selected_row_count = 0;
            this->auxiliary_view_disagreement_weight = 0;
            this->merged_candidate_count = 0;
            this->best_candidate_cost = std::numeric_limits<double>::infinity();
            this->auxiliary_view_disagreement = false;
            this->auxiliary_bp_converge = false;
            this->final_decoding.resize(this->bit_count, 0);
            this->osd0_l0_decoding.resize(this->bit_count, 0);
            this->osdw_l0_decoding.resize(this->bit_count, 0);
            this->osd0_l1_decoding.resize(this->bit_count, 0);
            this->osdw_l1_decoding.resize(this->bit_count, 0);
            this->bp_l0_decoding.resize(this->bit_count, 0);
            this->bp_l1_decoding.resize(this->bit_count, 0);

            this->osd_setup();


        }

        int osd_setup(){

            this->osd_candidate_strings.clear();

            if (this->LuDecomposition != nullptr) {
                delete this->LuDecomposition;
                this->LuDecomposition = nullptr;
            }

            if(this->osd_method == OSD_OFF) return 0;

            this->LuDecomposition = new ldpc_modified::gf2sparse_linalg::RowReduce<ldpc_modified::adaptive_bp::BpEntry>(this->pcm);
            this->column_ordering.resize(this->pcm.n);
            int osd_candidate_string_count;
            this->LuDecomposition->rref(false,true); 
            this->k = this->pcm.n - this->LuDecomposition->rank;

            if(this->osd_method == OSD_0 || this->osd_order==0){
                return 1;
            }

            if(this->osd_method == EXHAUSTIVE){
                osd_candidate_string_count = pow(2,this->osd_order);
                for(int i=1; i<osd_candidate_string_count; i++){
                    this->osd_candidate_strings.push_back(ldpc_modified::util::decimal_to_binary_reverse(i,k));
                }
            }

            if(this->osd_method == COMBINATION_SWEEP){
                for(int i=0; i<k; i++) {
                    std::vector<uint8_t> osd_candidate;
                    osd_candidate.resize(k,0);
                    osd_candidate[i]=1; 
                    this->osd_candidate_strings.push_back(osd_candidate);
                }

                for(int i = 0; i<this->osd_order;i++){
                    for(int j = 0; j<this->osd_order; j++){
                        if(j<=i) continue;
                        std::vector<uint8_t> osd_candidate;
                        osd_candidate.resize(k,0);
                        osd_candidate[i]=1;
                        osd_candidate[j]=1; 
                        this->osd_candidate_strings.push_back(osd_candidate);
                    }
                }

            }
            return 1;
        }

        ~AdaptiveOsdDecoder(){
            if (this->LuDecomposition != nullptr) {
                delete this->LuDecomposition;
            }
        };

        void reset_zero_syndrome_stats() {
            this->final_decoding.assign(this->bit_count, 0);
            this->osd0_decoding.assign(this->bit_count, 0);
            this->osdw_decoding.assign(this->bit_count, 0);
            this->osd0_l0_decoding.assign(this->bit_count, 0);
            this->osdw_l0_decoding.assign(this->bit_count, 0);
            this->osd0_l1_decoding.assign(this->bit_count, 0);
            this->osdw_l1_decoding.assign(this->bit_count, 0);
            this->bp_l0_decoding.assign(this->bit_count, 0);
            this->bp_l1_decoding.assign(this->bit_count, 0);
            this->auxiliary_candidate_pair_count = 0;
            this->auxiliary_selected_row_count = 0;
            this->auxiliary_view_disagreement_weight = 0;
            this->merged_candidate_count = 1;
            this->auxiliary_view_disagreement = false;
            this->auxiliary_bp_converge = true;

            this->best_candidate_cost = 0.0;
            for (int i = 0; i < this->bit_count; i++) {
                double p = this->channel_probabilities[i];
                assert(p > 0.0 && p < 1.0);
                this->best_candidate_cost += -std::log(1.0 - p);
            }
        }


        std::vector<uint8_t>& decode(std::vector<uint8_t>& syndrome, std::vector<double>& log_prob_ratios) {

            ldpc_modified::sort::soft_decision_col_sort(log_prob_ratios, this->column_ordering,bit_count);

            if(this->osd_order == 0){
                this->osd0_decoding = this->osdw_decoding =  this->LuDecomposition->fast_solve(syndrome,this->column_ordering);
                return this->osd0_decoding;
            }

            //row reduce the matrix according to the new column ordering
            this->LuDecomposition->rref(false,true,this->column_ordering);

            // find the OSD0 solution
            this->osd0_decoding = this->osdw_decoding = LuDecomposition->lu_solve(syndrome);

            // this->osd0_decoding = this->LuDecomposition->fast_solve(syndrome,this->column_ordering);

            // if(osd_order==0){
            //     return this->osd0_decoding;
            // }

            double candidate_weight, osd_min_weight;

            osd_min_weight=0;
            for(int i=0; i<this->pcm.n; i++){
                if(this->osd0_decoding[i]==1){
                    osd_min_weight+=log(1/this->channel_probabilities[i]);
                }
            }

            std::vector<int> non_pivot_columns;
            std::vector<ldpc_modified::gf2sparse::GF2Entry*> delete_entries;
            for(int i = this->LuDecomposition->rank; i<this->pcm.n; i++){
                int col = this->LuDecomposition->cols[i];
                non_pivot_columns.push_back(col);
                
                for(auto& e: this->LuDecomposition->U.iterate_column(col)){
                    delete_entries.push_back(&e);
                }
            }

            for (auto e: delete_entries){
                this->LuDecomposition->U.remove(*e);
            }
            
            for(auto& candidate_string: this->osd_candidate_strings){

                auto t_syndrome = syndrome;
                int col_index = 0;
                for(int col: non_pivot_columns){
                    if(candidate_string[col_index]==1){
                        for(auto& e: this->pcm.iterate_column(col)){
                            t_syndrome[e.row_index] ^= 1;
                        }
                    }
                    col_index++;
                }

                auto candidate_solution = this->LuDecomposition->lu_solve(t_syndrome);
                for(int i=0; i<k; i++){
                    candidate_solution[non_pivot_columns[i]]=candidate_string[i];
                }
                candidate_weight=0;
                for(int i=0; i<this->pcm.n; i++){
                    if(candidate_solution[i]==1){
                        candidate_weight+=log(1/this->channel_probabilities[i]);
                    }
                }
                if(candidate_weight<osd_min_weight){
                    osd_min_weight = candidate_weight;
                    
                    this->osdw_decoding = candidate_solution;

                }

            }

            return this->osdw_decoding;

        }

        std::vector<uint8_t>& adaptive_decode_from_llrs(
                std::vector<uint8_t>& syndrome,
                std::vector<double>& l0,
                std::vector<double>& l1,
                std::vector<uint8_t>& bp_l0_decoding,
                std::vector<uint8_t>& bp_l1_decoding,
                int auxiliary_candidate_pair_count,
                int auxiliary_selected_row_count,
                bool auxiliary_bp_converge) {

            if (syndrome.size() != size_t(this->check_count) ||
                    l0.size() != size_t(this->bit_count) ||
                    l1.size() != size_t(this->bit_count) ||
                    bp_l0_decoding.size() != size_t(this->bit_count) ||
                    bp_l1_decoding.size() != size_t(this->bit_count)) {
                throw std::runtime_error("Adaptive OSD received vectors with incompatible dimensions.");
            }

            this->merged_candidate_count = 0;
            this->best_candidate_cost = std::numeric_limits<double>::infinity();
            this->auxiliary_candidate_pair_count = auxiliary_candidate_pair_count;
            this->auxiliary_selected_row_count = auxiliary_selected_row_count;
            this->auxiliary_view_disagreement_weight = 0;
            this->auxiliary_view_disagreement = false;
            this->auxiliary_bp_converge = auxiliary_bp_converge;
            this->bp_l0_decoding = bp_l0_decoding;
            this->bp_l1_decoding = bp_l1_decoding;

            if (this->osd_method != OSD_OFF) {
                // OSD is still solved against the original H and original syndrome.
                // Only the reliability ordering changes between L(0) and L(1).
                this->decode(syndrome, l0);
                this->osd0_l0_decoding = this->osd0_decoding;
                this->osdw_l0_decoding = this->osdw_decoding;

                this->decode(syndrome, l1);
                this->osd0_l1_decoding = this->osd0_decoding;
                this->osdw_l1_decoding = this->osdw_decoding;
            } else {
                this->osd0_l0_decoding = this->bp_l0_decoding;
                this->osdw_l0_decoding = this->bp_l0_decoding;
                this->osd0_l1_decoding = this->bp_l1_decoding;
                this->osdw_l1_decoding = this->bp_l1_decoding;
            }

            this->auxiliary_view_disagreement_weight = 0;
            for (int i = 0; i < this->bit_count; i++) {
                if (this->osdw_l0_decoding[i] != this->osdw_l1_decoding[i]) {
                    this->auxiliary_view_disagreement_weight += 1;
                }
            }
            this->auxiliary_view_disagreement = this->auxiliary_view_disagreement_weight > 0;

            std::vector<std::vector<uint8_t>> merged_candidates;
            merged_candidates.push_back(this->bp_l0_decoding);
            merged_candidates.push_back(this->bp_l1_decoding);
            merged_candidates.push_back(this->osd0_l0_decoding);
            merged_candidates.push_back(this->osdw_l0_decoding);
            merged_candidates.push_back(this->osd0_l1_decoding);
            merged_candidates.push_back(this->osdw_l1_decoding);

            bool found_valid_candidate = false;
            for (auto &candidate: merged_candidates) {
                if (candidate.size() != this->bit_count) continue;

                auto candidate_syndrome = this->pcm.mulvec(candidate);
                if (!std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                    continue;
                }

                this->merged_candidate_count += 1;
                double candidate_cost = 0.0;
                for (int i = 0; i < this->bit_count; i++) {
                    double p = this->channel_probabilities[i];
                    assert(p > 0.0 && p < 1.0);
                    candidate_cost += candidate[i] ? -std::log(p) : -std::log(1.0 - p);
                }

                if (!found_valid_candidate || candidate_cost < this->best_candidate_cost) {
                    found_valid_candidate = true;
                    this->best_candidate_cost = candidate_cost;
                    this->final_decoding = candidate;
                }
            }

            // OSD should normally supply at least one valid candidate.  If OSD is
            // deliberately disabled, preserve the first BP view rather than hiding
            // the fact that no valid candidate was available.
            if (!found_valid_candidate) {
                this->final_decoding = this->bp_l0_decoding;
            }

            this->osd0_decoding = this->osd0_l0_decoding;
            this->osdw_decoding = this->final_decoding;
            return this->final_decoding;
        }

	};

}//end adaptive_osd namespace


#endif
