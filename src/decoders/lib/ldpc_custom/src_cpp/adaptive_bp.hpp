#ifndef ADAPTIVE_BP_H
#define ADAPTIVE_BP_H

#include <utility>
#include <vector>
#include <memory>
#include <iterator>
#include <cmath>
#include <limits>
#include <random>
#include <chrono>
#include <stdexcept> // required for std::runtime_error
#include <set>
#include <algorithm>
#include <cassert>
#include <tuple>

#include "math.h"
#include "sparse_matrix_base.hpp"
#include "gf2sparse.hpp"
#include "rng.hpp"

namespace ldpc_modified {
    namespace adaptive_bp {

        enum BpMethod {
            PRODUCT_SUM = 0,
            MINIMUM_SUM = 1
        };

        enum BpSchedule {
            SERIAL = 0,
            PARALLEL = 1,
            SERIAL_RELATIVE = 2
        };

        enum BpInputType {
            SYNDROME = 0,
            RECEIVED_VECTOR = 1,
            AUTO = 2
        };

        const std::vector<int> NULL_INT_VECTOR = {};

        class BpEntry : public ldpc_modified::sparse_matrix_base::EntryBase<BpEntry> {
        public:
            double bit_to_check_msg = 0.0;
            double check_to_bit_msg = 0.0;

            ~BpEntry() = default;
        };
        using BpSparse = ldpc_modified::gf2sparse::GF2Sparse<BpEntry>;

        class AdaptiveBpDecoder {
            // TODO properties should be private and only accessible via getters and setters
        public:
            BpSparse &pcm;
            std::vector<double> channel_probabilities;
            int check_count;
            int bit_count;
            int maximum_iterations;
            BpMethod bp_method;
            BpSchedule schedule;
            BpInputType bp_input_type;
            double ms_scaling_factor;
            std::vector<uint8_t> decoding;
            std::vector<uint8_t> candidate_syndrome;

            std::vector<double> log_prob_ratios;
            std::vector<double> initial_log_prob_ratios;
            std::vector<double> soft_syndrome;
            std::vector<int> serial_schedule_order;
            int iterations;
            int omp_thread_count;
            bool converge;
            int random_schedule_seed;
            bool random_serial_schedule;
            ldpc_modified::rng::RandomListShuffle<int> rng_list_shuffle;

            std::vector<int> sign_flip_count;
            std::vector<double> llr_min;
            std::vector<double> llr_max;
            std::vector<double> instability_scores;
            std::vector<int> auxiliary_row_a;
            std::vector<int> auxiliary_row_b;
            int auxiliary_candidate_pair_count;
            int auxiliary_selected_row_count;

            AdaptiveBpDecoder(
                    BpSparse &parity_check_matrix,
                    std::vector<double> channel_probabilities,
                    int maximum_iterations = 0,
                    BpMethod bp_method = PRODUCT_SUM,
                    BpSchedule schedule = PARALLEL,
                    double min_sum_scaling_factor = 0.625,
                    int omp_threads = 1,
                    const std::vector<int> &serial_schedule = NULL_INT_VECTOR,
                    int random_schedule_seed = 0,
                    bool random_serial_schedule = false,
                    BpInputType bp_input_type = AUTO) :
                    pcm(parity_check_matrix), channel_probabilities(std::move(channel_probabilities)),
                    check_count(pcm.m), bit_count(pcm.n), maximum_iterations(maximum_iterations), bp_method(bp_method),
                    schedule(schedule), ms_scaling_factor(min_sum_scaling_factor),
                    iterations(0) //the parity check matrix is passed in by reference
            {

                this->initial_log_prob_ratios.resize(bit_count);
                this->log_prob_ratios.resize(bit_count);
                this->candidate_syndrome.resize(check_count);
                this->decoding.resize(bit_count);
                this->converge = 0;
                this->omp_thread_count = omp_threads;
                this->random_schedule_seed = random_schedule_seed;
                this->random_serial_schedule = random_serial_schedule;
                this->bp_input_type = bp_input_type;

                this->sign_flip_count.resize(bit_count);
                this->llr_min.resize(bit_count);
                this->llr_max.resize(bit_count);
                this->instability_scores.resize(bit_count);
                this->auxiliary_candidate_pair_count = 0;
                this->auxiliary_selected_row_count = 0;


                if (this->channel_probabilities.size() != this->bit_count) {
                    throw std::runtime_error(
                            "Channel probabilities vector must have length equal to the number of bits");
                }
                this->refresh_initial_log_prob_ratios();
                if (serial_schedule != NULL_INT_VECTOR) {
                    this->serial_schedule_order = serial_schedule;
                    if (this->random_serial_schedule) {
                        throw std::runtime_error("Random schedule cannot be used with a fixed serial schedule. Set `random_serial_schedule` input parameter to false.");
                    }
                } else if (this->random_serial_schedule) {
                    this->serial_schedule_order.resize(bit_count);
                    for (int i = 0; i < bit_count; i++) {
                        this->serial_schedule_order[i] = i;
                    }
                    this->rng_list_shuffle.seed(this->random_schedule_seed);
                }
                else {
                    this->serial_schedule_order.resize(bit_count);
                    for (int i = 0; i < bit_count; i++) {
                        this->serial_schedule_order[i] = i;
                    }
                }

                //Initialise OMP thread pool
                // this->omp_thread_count = omp_threads;
                // this->set_omp_thread_count(this->omp_thread_count);
            }

            ~AdaptiveBpDecoder() = default;

            void refresh_initial_log_prob_ratios() {
                this->initial_log_prob_ratios.resize(this->bit_count);
                this->log_prob_ratios.resize(this->bit_count);
                for (int i = 0; i < this->bit_count; i++) {
                    this->initial_log_prob_ratios[i] = std::log(
                            (1 - this->channel_probabilities[i]) / this->channel_probabilities[i]);
                }
            }

            void update_decoder_parameters(
                    const std::vector<double> &channel_probabilities,
                    int maximum_iterations,
                    BpMethod bp_method,
                    BpSchedule schedule,
                    double min_sum_scaling_factor,
                    int omp_threads,
                    const std::vector<int> &serial_schedule = NULL_INT_VECTOR,
                    int random_schedule_seed = 0,
                    bool random_serial_schedule = false,
                    BpInputType bp_input_type = AUTO) {

                this->check_count = this->pcm.m;
                this->bit_count = this->pcm.n;
                this->maximum_iterations = maximum_iterations;
                this->bp_method = bp_method;
                this->schedule = schedule;
                this->ms_scaling_factor = min_sum_scaling_factor;
                this->omp_thread_count = omp_threads;
                this->random_schedule_seed = random_schedule_seed;
                this->random_serial_schedule = random_serial_schedule;
                this->bp_input_type = bp_input_type;
                this->iterations = 0;
                this->converge = false;

                if (channel_probabilities.size() != size_t(this->bit_count)) {
                    throw std::runtime_error(
                            "Channel probabilities vector must have length equal to the number of bits");
                }

                bool channel_changed =
                        this->channel_probabilities.size() != channel_probabilities.size() ||
                        !std::equal(
                                this->channel_probabilities.begin(),
                                this->channel_probabilities.end(),
                                channel_probabilities.begin());
                if (channel_changed) {
                    this->channel_probabilities = channel_probabilities;
                }

                this->initial_log_prob_ratios.resize(this->bit_count);
                this->log_prob_ratios.resize(this->bit_count);
                this->candidate_syndrome.assign(this->check_count, 0);
                this->decoding.assign(this->bit_count, 0);
                this->soft_syndrome.resize(this->check_count);
                this->sign_flip_count.assign(this->bit_count, 0);
                this->llr_min.resize(this->bit_count);
                this->llr_max.resize(this->bit_count);
                this->instability_scores.assign(this->bit_count, 0.0);
                if (channel_changed) {
                    this->refresh_initial_log_prob_ratios();
                }

                if (serial_schedule != NULL_INT_VECTOR) {
                    this->serial_schedule_order = serial_schedule;
                    if (this->random_serial_schedule) {
                        throw std::runtime_error("Random schedule cannot be used with a fixed serial schedule. Set `random_serial_schedule` input parameter to false.");
                    }
                } else {
                    this->serial_schedule_order.resize(this->bit_count);
                    for (int i = 0; i < this->bit_count; i++) {
                        this->serial_schedule_order[i] = i;
                    }
                    if (this->random_serial_schedule) {
                        this->rng_list_shuffle.seed(this->random_schedule_seed);
                    }
                }
            }

            void set_omp_thread_count(int count) {
                this->omp_thread_count = count;
                // omp_set_num_threads(this->omp_thread_count);
                // NotImplemented
            }

            void set_random_schedule_seed(int seed) {
                this->random_schedule_seed = seed;
                this->rng_list_shuffle.seed(seed);
            }

            void initialise_log_domain_bp() {
                // initialise BP
                for (int i = 0; i < this->bit_count; i++) {
                    for (auto &e: this->pcm.iterate_column(i)) {
                        e.bit_to_check_msg = this->initial_log_prob_ratios[i];
                    }
                }

                // These per-variable counters describe one BP run.  They are later
                // reused by the auxiliary-view constructor to find variables whose
                // belief is oscillatory or weak.
                this->decoding.assign(this->bit_count, 0);
                this->sign_flip_count.assign(this->bit_count, 0);
                this->llr_max.assign(this->bit_count, std::numeric_limits<double>::lowest());
                this->llr_min.assign(this->bit_count, std::numeric_limits<double>::max());
            }

            void reset_zero_syndrome_state() {
                this->converge = true;
                this->iterations = 0;
                this->candidate_syndrome.assign(this->check_count, 0);
                this->decoding.assign(this->bit_count, 0);
                this->sign_flip_count.assign(this->bit_count, 0);

                for (int i = 0; i < this->bit_count; i++) {
                    double prior_llr = this->initial_log_prob_ratios[i];
                    this->initial_log_prob_ratios[i] = prior_llr;
                    this->log_prob_ratios[i] = prior_llr;
                    this->llr_min[i] = prior_llr;
                    this->llr_max[i] = prior_llr;
                    this->instability_scores[i] = 0.0;
                    for (auto &e: this->pcm.iterate_column(i)) {
                        e.bit_to_check_msg = prior_llr;
                        e.check_to_bit_msg = 0.0;
                    }
                }
            }

            void build_adaptive_redundant_view(
                    std::vector<uint8_t> &syndrome,
                    int adaptive_k,
                    int adaptive_w_max,
                    int adaptive_B,
                    BpSparse &augmented_pcm,
                    std::vector<uint8_t> &augmented_syndrome) {

                if (adaptive_k < 0 || adaptive_w_max < 0 || adaptive_B < 0) {
                    throw std::runtime_error("Adaptive hyperparameters must be non-negative.");
                }

                this->auxiliary_row_a.clear();
                this->auxiliary_row_b.clear();
                this->auxiliary_candidate_pair_count = 0;
                this->auxiliary_selected_row_count = 0;
                augmented_syndrome.assign(this->check_count, 0);

                // score[v] combines three signs of BP instability:
                //   1) hard-decision sign flips across iterations;
                //   2) LLR swing range, which catches high-amplitude oscillation;
                //   3) weak final belief 1/(|L(0)|+eps), which catches near-ties.
                // The sign-flip term is multiplicative so repeatedly changing
                // variables dominate merely low-confidence but stable variables.
                for (int v = 0; v < this->bit_count; v++) {
                    double llr_range = this->llr_max[v] - this->llr_min[v];
                    double weak_belief = 1.0 / (std::abs(this->log_prob_ratios[v]) + 1e-12);
                    this->instability_scores[v] =
                            (1.0 + double(this->sign_flip_count[v])) * (1.0 + llr_range) + weak_belief;
                }

                int unstable_count = std::min(adaptive_k, this->bit_count);
                std::vector<int> variable_order(this->bit_count);
                for (int v = 0; v < this->bit_count; v++) variable_order[v] = v;
                auto variable_is_more_unstable = [this](int a, int b) {
                    if (this->instability_scores[a] == this->instability_scores[b]) {
                        return std::abs(this->log_prob_ratios[a]) < std::abs(this->log_prob_ratios[b]);
                    }
                    return this->instability_scores[a] > this->instability_scores[b];
                };
                if (unstable_count > 0 && unstable_count < this->bit_count) {
                    std::nth_element(
                            variable_order.begin(),
                            variable_order.begin() + unstable_count,
                            variable_order.end(),
                            variable_is_more_unstable);
                    std::sort(
                            variable_order.begin(),
                            variable_order.begin() + unstable_count,
                            variable_is_more_unstable);
                } else if (unstable_count == this->bit_count) {
                    std::sort(variable_order.begin(), variable_order.end(), variable_is_more_unstable);
                }

                std::vector<uint8_t> unstable_variable_flag(this->bit_count, 0);
                std::vector<uint8_t> neighboring_check_flag(this->check_count, 0);
                std::vector<int> neighboring_checks;
                for (int idx = 0; idx < unstable_count; idx++) {
                    int v = variable_order[idx];
                    unstable_variable_flag[v] = 1;
                    for (auto &e: this->pcm.iterate_column(v)) {
                        if (!neighboring_check_flag[e.row_index]) {
                            neighboring_check_flag[e.row_index] = 1;
                            neighboring_checks.push_back(e.row_index);
                        }
                    }
                }

                struct RedundantRowCandidate {
                    double score;
                    int a;
                    int b;
                    std::vector<int> support;
                };

                std::vector<RedundantRowCandidate> candidates;
                std::vector<uint8_t> row_marks(this->bit_count, 0);
                std::vector<int> touched_columns;

                for (int ia = 0; ia < int(neighboring_checks.size()); ia++) {
                    for (int ib = ia + 1; ib < int(neighboring_checks.size()); ib++) {
                        int a = neighboring_checks[ia];
                        int b = neighboring_checks[ib];
                        touched_columns.clear();

                        // Construct r_ab = H_a xor H_b by toggling sparse support marks.
                        for (auto &e: this->pcm.iterate_row(a)) {
                            if (!row_marks[e.col_index]) touched_columns.push_back(e.col_index);
                            row_marks[e.col_index] ^= 1;
                        }
                        for (auto &e: this->pcm.iterate_row(b)) {
                            if (!row_marks[e.col_index]) touched_columns.push_back(e.col_index);
                            row_marks[e.col_index] ^= 1;
                        }

                        std::vector<int> support;
                        double score = 0.0;
                        double support_tie_break_score = 0.0;
                        for (int col: touched_columns) {
                            if (row_marks[col]) {
                                support.push_back(col);
                                support_tie_break_score += this->instability_scores[col];
                                if (unstable_variable_flag[col]) {
                                    // The main score is the overlap with the top-k unstable
                                    // variables that survive in r_ab.  Instability outside
                                    // that overlap is only a tiny deterministic tie-break.
                                    score += this->instability_scores[col];
                                }
                            }
                            row_marks[col] = 0;
                        }
                        score += 1e-9 * support_tie_break_score;

                        int row_weight = int(support.size());
                        if (row_weight > 0 && row_weight <= adaptive_w_max) {
                            candidates.push_back({score, a, b, std::move(support)});
                        }
                    }
                }

                this->auxiliary_candidate_pair_count = int(candidates.size());
                int selected_count = std::min(adaptive_B, int(candidates.size()));
                this->auxiliary_selected_row_count = selected_count;
                auto candidate_is_better = [](const RedundantRowCandidate &x, const RedundantRowCandidate &y) {
                    if (x.score == y.score) return x.support.size() < y.support.size();
                    return x.score > y.score;
                };
                if (selected_count > 0 && selected_count < int(candidates.size())) {
                    std::nth_element(
                            candidates.begin(),
                            candidates.begin() + selected_count,
                            candidates.end(),
                            candidate_is_better);
                    std::sort(
                            candidates.begin(),
                            candidates.begin() + selected_count,
                            candidate_is_better);
                } else if (selected_count == int(candidates.size())) {
                    std::sort(candidates.begin(), candidates.end(), candidate_is_better);
                }

                int augmented_entry_count = this->pcm.entry_count();
                for (int i = 0; i < selected_count; i++) {
                    augmented_entry_count += int(candidates[i].support.size());
                    this->auxiliary_row_a.push_back(candidates[i].a);
                    this->auxiliary_row_b.push_back(candidates[i].b);
                }

                augmented_pcm.allocate(this->check_count + selected_count, this->bit_count, augmented_entry_count);
                augmented_syndrome.assign(this->check_count + selected_count, 0);

                // Copy the original checks unchanged.
                for (int row = 0; row < this->check_count; row++) {
                    augmented_syndrome[row] = syndrome[row];
                    for (auto &e: this->pcm.iterate_row(row)) {
                        augmented_pcm.insert_entry(row, e.col_index);
                    }
                }

                // Append the selected redundant checks.  The syndrome of H_a xor H_b
                // is syndrome[a] xor syndrome[b], so no extra matrix-vector product is
                // needed here.
                for (int i = 0; i < selected_count; i++) {
                    int row = this->check_count + i;
                    augmented_syndrome[row] = syndrome[candidates[i].a] ^ syndrome[candidates[i].b];
                    for (int col: candidates[i].support) {
                        augmented_pcm.insert_entry(row, col);
                    }
                }
            }

            std::vector<uint8_t> decode(std::vector<uint8_t> &input_vector) {


                if ((this->bp_input_type == AUTO && input_vector.size() == this->bit_count) ||
                    this->bp_input_type == RECEIVED_VECTOR) {
                    auto syndrome = pcm.mulvec(input_vector);
                    std::vector<uint8_t> rv_decoding;
                    if (schedule == PARALLEL) {
                        rv_decoding = bp_decode_parallel(syndrome);
                    } else if (schedule == SERIAL || schedule == SERIAL_RELATIVE) {
                        rv_decoding = bp_decode_serial(syndrome);
                    } else {
                        throw std::runtime_error("Invalid BP schedule");
                    }

                    for (int i = 0; i < this->bit_count; i++) {
                        this->decoding[i] = rv_decoding[i] ^ input_vector[i];
                    }

                    return this->decoding;

                }


                if (schedule == PARALLEL) {
                    return bp_decode_parallel(input_vector);
                }
                if (schedule == SERIAL || schedule == SERIAL_RELATIVE) {
                    return bp_decode_serial(input_vector);
                } else { throw std::runtime_error("Invalid BP schedule"); }

            }

            std::vector<uint8_t> &bp_decode_parallel(std::vector<uint8_t> &syndrome) {

                this->converge = 0;

                this->initialise_log_domain_bp();

                //main interation loop
                for (int it = 1; it <= this->maximum_iterations; it++) {

                    if (this->bp_method == PRODUCT_SUM) {
                        for (int i = 0; i < this->check_count; i++) {
                            this->candidate_syndrome[i] = 0;

                            double temp = 1.0;
                            for (auto &e: this->pcm.iterate_row(i)) {
                                e.check_to_bit_msg = temp;
                                temp *= std::tanh(e.bit_to_check_msg / 2);
                            }

                            temp = 1;
                            for (auto &e: this->pcm.reverse_iterate_row(i)) {
                                e.check_to_bit_msg *= temp;
                                int message_sign = syndrome[i] != 0u ? -1.0 : 1.0;
                                e.check_to_bit_msg =
                                        message_sign * std::log((1 + e.check_to_bit_msg) / (1 - e.check_to_bit_msg));
                                temp *= std::tanh(e.bit_to_check_msg / 2);
                            }
                        }
                    } else if (this->bp_method == MINIMUM_SUM) {

                        double alpha;
                        if(this->ms_scaling_factor == 0.0) {
                            alpha = 1.0 - std::pow(2.0, -1.0*it);
                        }
                        else {
                            alpha = this->ms_scaling_factor;
                        }

                        //check to bit updates
                        for (int i = 0; i < check_count; i++) {

                            this->candidate_syndrome[i] = 0;
                            int total_sgn = 0;
                            int sgn = 0;
                            total_sgn = syndrome[i];
                            double temp = std::numeric_limits<double>::max();

                            for (auto &e: this->pcm.iterate_row(i)) {
                                if (e.bit_to_check_msg <= 0) {
                                    total_sgn += 1;
                                }
                                e.check_to_bit_msg = temp;
                                double abs_bit_to_check_msg = std::abs(e.bit_to_check_msg);
                                if (abs_bit_to_check_msg < temp) {
                                    temp = abs_bit_to_check_msg;
                                }
                            }

                            temp = std::numeric_limits<double>::max();
                            for (auto &e: this->pcm.reverse_iterate_row(i)) {
                                sgn = total_sgn;
                                if (e.bit_to_check_msg <= 0) {
                                    sgn += 1;
                                }
                                if (temp < e.check_to_bit_msg) {
                                    e.check_to_bit_msg = temp;
                                }

                                int message_sign = (sgn % 2 == 0) ? 1.0 : -1.0;
                                
                                e.check_to_bit_msg *= message_sign * alpha;

                                
                                double abs_bit_to_check_msg = std::abs(e.bit_to_check_msg);
                                if (abs_bit_to_check_msg < temp) {
                                    temp = abs_bit_to_check_msg;
                                }

                            }

                        }
                    }


                    //compute log probability ratios
                    for (int i = 0; i < this->bit_count; i++) {
                        double temp = initial_log_prob_ratios[i];
                        for (auto &e: this->pcm.iterate_column(i)) {
                            e.bit_to_check_msg = temp;
                            temp += e.check_to_bit_msg;
                            // if(isnan(temp)) temp = e.bit_to_check_msg;


                        }

                        //make hard decision on basis of log probability ratio for bit i
                        this->log_prob_ratios[i] = temp;
                        this->llr_max[i] = std::max(this->llr_max[i], temp);
                        this->llr_min[i] = std::min(this->llr_min[i], temp);
                        // if(isnan(log_prob_ratios[i])) log_prob_ratios[i] = initial_log_prob_ratios[i];
                        if (temp <= 0) {
                            if (this->decoding[i] == 0) this->sign_flip_count[i] += 1;
                            this->decoding[i] = 1;
                            for (auto &e: this->pcm.iterate_column(i)) {
                                this->candidate_syndrome[e.row_index] ^= 1;
                            }
                        } else {
                            if (this->decoding[i] == 1) this->sign_flip_count[i] += 1;
                            this->decoding[i] = 0;
                        }
                    }

                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        this->converge = true;
                    }

                    this->iterations = it;

                    if (this->converge) {
                        return this->decoding;
                    }


                    //compute bit to check update
                    for (int i = 0; i < bit_count; i++) {
                        double temp = 0;
                        for (auto &e: this->pcm.reverse_iterate_column(i)) {
                            e.bit_to_check_msg += temp;
                            temp += e.check_to_bit_msg;
                        }
                    }

                }


                return this->decoding;

            }

            std::vector<uint8_t> &bp_decode_single_scan(std::vector<uint8_t> &syndrome) {

                converge = 0;
                int CONVERGED = 0;

                std::vector<double> log_prob_ratios_old;
                log_prob_ratios_old.resize(bit_count);

                for (int i = 0; i < bit_count; i++) {
                    this->initial_log_prob_ratios[i] = std::log(
                            (1 - this->channel_probabilities[i]) / this->channel_probabilities[i]);
                    this->log_prob_ratios[i] = this->initial_log_prob_ratios[i];

                }

                // initialise_log_domain_bp();

                //main interation loop
                for (int it = 1; it <= maximum_iterations; it++) {

                    if (CONVERGED != 0) {
                        continue;
                    }

                    // std::fill(candidate_syndrome.begin(), candidate_syndrome.end(), 0);

                    log_prob_ratios_old = this->log_prob_ratios;

                    if (it != 1) {
                        this->log_prob_ratios = this->initial_log_prob_ratios;
                    }

                    //check to bit updates
                    for (int i = 0; i < check_count; i++) {

                        this->candidate_syndrome[i] = 0;

                        int total_sgn = 0;
                        int sgn = 0;
                        total_sgn = syndrome[i];
                        double temp = std::numeric_limits<double>::max();

                        double bit_to_check_msg = NAN;

                        for (auto &e: pcm.iterate_row(i)) {
                            if (it == 1) {
                                e.check_to_bit_msg = 0;
                            }
                            bit_to_check_msg = log_prob_ratios_old[e.col_index] - e.check_to_bit_msg;
                            if (bit_to_check_msg <= 0) {
                                total_sgn += 1;
                            }
                            e.bit_to_check_msg = temp;
                            double abs_bit_to_check_msg = std::abs(bit_to_check_msg);
                            if (abs_bit_to_check_msg < temp) {
                                temp = abs_bit_to_check_msg;
                            }
                        }

                        temp = std::numeric_limits<double>::max();
                        for (auto &e: pcm.reverse_iterate_row(i)) {
                            sgn = total_sgn;
                            if (it == 1) {
                                e.check_to_bit_msg = 0;
                            }
                            bit_to_check_msg = log_prob_ratios_old[e.col_index] - e.check_to_bit_msg;
                            if (bit_to_check_msg <= 0) {
                                sgn += 1;
                            }
                            if (temp < e.bit_to_check_msg) {
                                e.bit_to_check_msg = temp;
                            }

                            int message_sign = (sgn % 2 == 0) ? 1.0 : -1.0;
                            e.check_to_bit_msg = message_sign * ms_scaling_factor * e.bit_to_check_msg;
                            this->log_prob_ratios[e.col_index] += e.check_to_bit_msg;


                            double abs_bit_to_check_msg = std::abs(bit_to_check_msg);
                            if (abs_bit_to_check_msg < temp) {
                                temp = abs_bit_to_check_msg;
                            }

                        }


                    }



                    //compute hard decisions and calculate syndrome
                    for (int i = 0; i < bit_count; i++) {
                        if (this->log_prob_ratios[i] <= 0) {
                            this->decoding[i] = 1;
                            for (auto &e: pcm.iterate_column(i)) {
                                this->candidate_syndrome[e.row_index] ^= 1;
                            }
                        } else {
                            this->decoding[i] = 0;
                        }
                    }

                    int loop_break = 0;
                    CONVERGED = 0;

                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        CONVERGED = 1;
                    }

                    iterations = it;

                    if (CONVERGED != 0) {
                        converge = (CONVERGED != 0);
                        return decoding;
                    }

                }


                converge = (CONVERGED != 0);
                return decoding;

            }

            std::vector<uint8_t> &bp_decode_serial(std::vector<uint8_t> &syndrome) {
                int check_index = 0;
                this->converge = false;
                // initialise BP
                this->initialise_log_domain_bp();

                for (int it = 1; it <= maximum_iterations; it++) {

                    double alpha;
                    if(this->ms_scaling_factor == 0.0) {
                        alpha = 1.0 - std::pow(2.0, -1.0*it);
                    }
                    else {
                        alpha = this->ms_scaling_factor;
                    }

                    if (this->random_serial_schedule) {
                        this->rng_list_shuffle.shuffle(this->serial_schedule_order);
                    } else if (this->schedule == BpSchedule::SERIAL_RELATIVE) {
                        // resort by LLRs in each iteration to ensure that the most reliable bits are considered first
                        std::sort(this->serial_schedule_order.begin(), this->serial_schedule_order.end(),
                                  [this, it](int bit1, int bit2) {
                                      if (it != 1) {
                                          return this->log_prob_ratios[bit1] > this->log_prob_ratios[bit2];
                                      } else {
                                          return std::log(
                                                  (1 - channel_probabilities[bit1]) / channel_probabilities[bit1]) >
                                                 std::log((1 - channel_probabilities[bit2]) /
                                                          channel_probabilities[bit2]);
                                      }
                                  });
                    }

                    for (int bit_index: this->serial_schedule_order) {
                        double temp = NAN;
                        this->log_prob_ratios[bit_index] = std::log(
                                (1 - channel_probabilities[bit_index]) / channel_probabilities[bit_index]);
                        if (this->bp_method == 0) {
                            for (auto &e: this->pcm.iterate_column(bit_index)) {
                                check_index = e.row_index;
                                e.check_to_bit_msg = 1.0;
                                for (auto &g: this->pcm.iterate_row(check_index)) {
                                    if (&g != &e) {
                                        e.check_to_bit_msg *= tanh(g.bit_to_check_msg / 2);
                                    }
                                }
                                e.check_to_bit_msg = pow(-1, syndrome[check_index]) *
                                                     std::log((1 + e.check_to_bit_msg) / (1 - e.check_to_bit_msg));
                                e.bit_to_check_msg = log_prob_ratios[bit_index];
                                this->log_prob_ratios[bit_index] += e.check_to_bit_msg;
                            }
                        } else if (this->bp_method == 1) {
                            for (auto &e: pcm.iterate_column(bit_index)) {
                                check_index = e.row_index;
                                int sgn = syndrome[check_index];
                                temp = std::numeric_limits<double>::max();
                                for (auto &g: this->pcm.iterate_row(check_index)) {
                                    if (&g != &e) {
                                        double abs_bit_to_check_msg = std::abs(g.bit_to_check_msg);
                                        if (abs_bit_to_check_msg < temp) {
                                            temp = abs_bit_to_check_msg;
                                        }
                                        if (g.bit_to_check_msg <= 0) {
                                            sgn += 1;
                                        }
                                    }
                                }
                                double message_sign = (sgn % 2 == 0) ? 1.0 : -1.0;
                                e.check_to_bit_msg = alpha * message_sign * temp;
                                e.bit_to_check_msg = log_prob_ratios[bit_index];
                                this->log_prob_ratios[bit_index] += e.check_to_bit_msg;
                            }
                        }
                        this->llr_max[bit_index] = std::max(this->llr_max[bit_index], this->log_prob_ratios[bit_index]);
                        this->llr_min[bit_index] = std::min(this->llr_min[bit_index], this->log_prob_ratios[bit_index]);

                        if (this->log_prob_ratios[bit_index] <= 0) {
                            if (this->decoding[bit_index] == 0) this->sign_flip_count[bit_index] += 1;
                            this->decoding[bit_index] = 1;
                        } else {
                            if (this->decoding[bit_index] == 1) this->sign_flip_count[bit_index] += 1;
                            this->decoding[bit_index] = 0;
                        }
                        temp = 0;
                        for (auto &e: this->pcm.reverse_iterate_column(bit_index)) {
                            e.bit_to_check_msg += temp;
                            temp += e.check_to_bit_msg;
                        }
                    }

                    // compute the syndrome for the current candidate decoding solution
                    this->candidate_syndrome = pcm.mulvec(decoding, candidate_syndrome);
                    this->iterations = it;
                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        this->converge = true;
                        return this->decoding;
                    }
                }
                return this->decoding;
            }

            std::vector<uint8_t> &
            soft_info_decode_serial(std::vector<double> &soft_info_syndrome, double cutoff, double sigma) {
                // compute the syndrome log-likelihoods and initialize hard syndrome
                std::vector<uint8_t> syndrome;
                this->soft_syndrome = soft_info_syndrome;
                for (int i = 0; i < this->check_count; i++) {
                    this->soft_syndrome[i] = 2 * this->soft_syndrome[i] / (sigma * sigma);
                    if (this->soft_syndrome[i] <= 0) {
                        syndrome.push_back(1);
                    } else {
                        syndrome.push_back(0);
                    }
                }

                int check_index = 0;
                this->converge = false;
                bool CONVERGED = false;
                bool loop_break = false;
                // initialise BP
                this->initialise_log_domain_bp();
                std::set<int> check_indices_updated;

                for (int it = 1; it <= maximum_iterations; it++) {
                    if (CONVERGED) {
                        continue;
                    }
                    if (this->random_serial_schedule && omp_thread_count == 1) {
                        // reorder schedule elements randomly
                        shuffle(serial_schedule_order.begin(), serial_schedule_order.end(),
                                std::default_random_engine(random_schedule_seed));
                    }

                    check_indices_updated.clear();
                    for (auto bit_index: serial_schedule_order) {
                        double temp = NAN;
                        log_prob_ratios[bit_index] = std::log(
                                (1 - channel_probabilities[bit_index]) / channel_probabilities[bit_index]);
                        for (auto &check_nbr: pcm.iterate_column(bit_index)) {
                            // first, we compute the min absolute value of neighbours excluding the current recipient
                            check_index = check_nbr.row_index;
                            int sgn = 0;
                            temp = std::numeric_limits<double>::max();
                            for (auto &g: pcm.iterate_row(check_index)) {
                                if (&g != &check_nbr) {
                                    if (std::abs(g.bit_to_check_msg) < temp) {
                                        temp = std::abs(g.bit_to_check_msg);
                                    }
                                    if (g.bit_to_check_msg <= 0) {
                                        sgn ^= 1;
                                    }
                                }
                            }
                            double min_bit_to_check_msg = temp;
                            double propagated_msg = min_bit_to_check_msg;
                            double soft_syndrome_magnitude = std::abs(this->soft_syndrome[check_index]);

                            // if the soft syndrome magnitude is below cutoff, we apply the virtual update rules
                            if (soft_syndrome_magnitude < cutoff) {
                                if (soft_syndrome_magnitude < std::abs(min_bit_to_check_msg)) {
                                    propagated_msg = soft_syndrome_magnitude;
                                    int check_node_sgn = sgn;
                                    if (check_nbr.bit_to_check_msg <= 0) {
                                        check_node_sgn ^= 1;
                                    }
                                    // now we check whether we have to update the soft syndrome magnitude and sign
                                    if (check_node_sgn == syndrome[check_index]) {
                                        if (std::abs(check_nbr.bit_to_check_msg) < min_bit_to_check_msg) {
                                            this->soft_syndrome[check_index] =
                                                    pow(-1, syndrome[check_index]) *
                                                    std::abs(check_nbr.bit_to_check_msg);
                                        } else {
                                            this->soft_syndrome[check_index] =
                                                    pow(-1, syndrome[check_index]) * min_bit_to_check_msg;
                                        }
                                    } else {
                                        syndrome[check_index] ^= 1;
                                        this->soft_syndrome[check_index] *= -1;
                                    }
                                }
                            }
                            sgn ^= syndrome[check_index];
                            check_nbr.check_to_bit_msg = ms_scaling_factor * pow(-1, sgn) * propagated_msg;
                            check_nbr.bit_to_check_msg = log_prob_ratios[bit_index];
                            log_prob_ratios[bit_index] += check_nbr.check_to_bit_msg;
                        }
                        // hard decision on bit
                        if (log_prob_ratios[bit_index] <= 0) {
                            decoding[bit_index] = 1;
                        } else {
                            decoding[bit_index] = 0;
                        }
                        temp = 0;
                        for (auto &e: pcm.reverse_iterate_column(bit_index)) {
                            e.bit_to_check_msg += temp;
                            temp += e.check_to_bit_msg;
                        }
                    }
                    // compute the syndrome for the current candidate decoding solution
                    loop_break = false;
                    CONVERGED = true;
                    for (auto i = 0; i < soft_info_syndrome.size(); i++) {
                        if (soft_info_syndrome[i] <= 0) {
                            candidate_syndrome[i] = 1;
                        } else {
                            candidate_syndrome[i] = 0;
                        }
                    }
                    candidate_syndrome = pcm.mulvec(decoding, candidate_syndrome);
                    for (auto i = 0; i < check_count && !loop_break; i++) {
                        if (candidate_syndrome[i] != syndrome[i]) {
                            CONVERGED = false;
                            loop_break = true;
                        }
                    }
                    iterations = it;
                }
                converge = CONVERGED;
                return decoding;
            }
        };
    }
}  // namespace ldpc_modified::adaptive_bp

#endif
