#ifndef EXPANDED_ADAPTIVE_BP_H
#define EXPANDED_ADAPTIVE_BP_H

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
#include <atomic>
#include <tuple>
#include <thread>
#include <mutex>
#include <exception>

#include "math.h"
#include "sparse_matrix_base.hpp"
#include "gf2sparse.hpp"
#include "rng.hpp"

#ifndef LDPC_MODIFIED_ADAPTIVE_BP_USE_MAX_OVERLAP
#define LDPC_MODIFIED_ADAPTIVE_BP_USE_MAX_OVERLAP 1
#endif

#ifndef LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE
#define LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE 1
#endif

#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE < 0 || LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE > 1
#error "LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE must be 0 (legacy) or 1 (dynamic-spatial)."
#endif

namespace ldpc_modified {
    namespace expanded_adaptive_bp {

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

        inline std::atomic<int> &expanded_adaptive_bp_live_instance_count_ref() {
            static std::atomic<int> count{0};
            return count;
        }

        inline int expanded_adaptive_bp_live_instance_count() {
            return expanded_adaptive_bp_live_instance_count_ref().load(std::memory_order_relaxed);
        }

        class BpEntry : public ldpc_modified::sparse_matrix_base::EntryBase<BpEntry> {
        public:
            double bit_to_check_msg = 0.0;
            double check_to_bit_msg = 0.0;

            ~BpEntry() = default;
        };
        using BpSparse = ldpc_modified::gf2sparse::GF2Sparse<BpEntry>;

        class ExpandedAdaptiveBpDecoder {
            // TODO properties should be private and only accessible via getters and setters
        public:
            struct RedundantRowCandidate {
                int row_a;
                int row_b;
                bool residual_relevant;
                std::vector<int> support;
                std::vector<int> unstable_support;
            };

            struct PreparedAdaptiveView {
                std::vector<int> candidate_indices;
                int automatic_weight_limit;
            };

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
            std::vector<double> instability_flip_component;
            std::vector<double> instability_movement_component;
            std::vector<double> instability_message_conflict_component;
            std::vector<double> instability_residual_component;
            std::vector<double> instability_weak_belief_component;
            std::vector<double> instability_spatial_component;
            std::vector<double> instability_seed_scores;
            std::vector<uint8_t> previous_hard_decision;
            std::vector<double> previous_iteration_llr;
            std::vector<double> weighted_flip_sum;
            std::vector<double> weighted_llr_movement_sum;
            std::vector<double> weighted_message_conflict_sum;
            std::vector<double> weighted_residual_pressure_sum;
            std::vector<uint8_t> latest_check_residual;
            double dynamics_transition_weight_sum;
            double dynamics_sample_weight_sum;
            bool has_previous_dynamics_sample;
            bool instability_scores_pending;

            std::vector<RedundantRowCandidate> adaptive_candidates;
            std::vector<PreparedAdaptiveView> adaptive_views;
            std::vector<int> adaptive_region_variables;
            std::vector<uint8_t> adaptive_final_decoding;
            int adaptive_executed_view_count;
            int adaptive_main_bp_run_count;
            int adaptive_core_variable_count;
            int adaptive_unstable_variable_count;
            int adaptive_feasible_unstable_variable_count;
            int adaptive_covered_unstable_variable_count;
            int adaptive_candidate_count;
            int adaptive_valid_candidate_count;
            int adaptive_selected_view;
            int adaptive_fallback_residual_weight;
            double adaptive_best_candidate_cost;
            bool adaptive_selected_candidate_valid;

            ExpandedAdaptiveBpDecoder(
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
                this->resize_instability_state();
                this->instability_scores_pending = false;
                this->adaptive_final_decoding.assign(bit_count, 0);
                this->adaptive_executed_view_count = 0;
                this->adaptive_main_bp_run_count = 0;
                this->adaptive_core_variable_count = 0;
                this->adaptive_unstable_variable_count = 0;
                this->adaptive_feasible_unstable_variable_count = 0;
                this->adaptive_covered_unstable_variable_count = 0;
                this->adaptive_candidate_count = 0;
                this->adaptive_valid_candidate_count = 0;
                this->adaptive_selected_view = -1;
                this->adaptive_fallback_residual_weight = 0;
                this->adaptive_best_candidate_cost = std::numeric_limits<double>::infinity();
                this->adaptive_selected_candidate_valid = false;


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
                expanded_adaptive_bp_live_instance_count_ref().fetch_add(1, std::memory_order_relaxed);
            }

            ~ExpandedAdaptiveBpDecoder() {
                expanded_adaptive_bp_live_instance_count_ref().fetch_sub(1, std::memory_order_relaxed);
            }

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
                this->sign_flip_count.assign(this->bit_count, 0);
                this->llr_min.resize(this->bit_count);
                this->llr_max.resize(this->bit_count);
                this->instability_scores.assign(this->bit_count, 0.0);
                this->resize_instability_state();
                this->instability_scores_pending = false;
                this->adaptive_candidates.clear();
                this->adaptive_views.clear();
                this->adaptive_region_variables.clear();
                this->adaptive_final_decoding.assign(this->bit_count, 0);
                this->adaptive_executed_view_count = 0;
                this->adaptive_main_bp_run_count = 0;
                this->adaptive_core_variable_count = 0;
                this->adaptive_unstable_variable_count = 0;
                this->adaptive_feasible_unstable_variable_count = 0;
                this->adaptive_covered_unstable_variable_count = 0;
                this->adaptive_candidate_count = 0;
                this->adaptive_valid_candidate_count = 0;
                this->adaptive_selected_view = -1;
                this->adaptive_fallback_residual_weight = 0;
                this->adaptive_best_candidate_cost = std::numeric_limits<double>::infinity();
                this->adaptive_selected_candidate_valid = false;
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
                this->reset_instability_dynamics();
                this->instability_scores_pending = true;
            }

            void reset_zero_syndrome_state() {
                this->converge = true;
                this->iterations = 0;
                this->candidate_syndrome.assign(this->check_count, 0);
                this->decoding.assign(this->bit_count, 0);
                this->reset_instability_dynamics();

                for (int i = 0; i < this->bit_count; i++) {
                    double prior_llr = this->initial_log_prob_ratios[i];
                    this->log_prob_ratios[i] = prior_llr;
                    this->llr_min[i] = prior_llr;
                    this->llr_max[i] = prior_llr;
                    this->instability_scores[i] = 0.0;
                    for (auto &e: this->pcm.iterate_column(i)) {
                        e.bit_to_check_msg = prior_llr;
                        e.check_to_bit_msg = 0.0;
                    }
                }
                this->instability_scores_pending = false;
            }

            void reset_adaptive_multi_view_state() {
                this->adaptive_candidates.clear();
                this->adaptive_views.clear();
                this->adaptive_region_variables.clear();
                this->adaptive_executed_view_count = 0;
                this->adaptive_main_bp_run_count = 0;
                this->adaptive_core_variable_count = 0;
                this->adaptive_unstable_variable_count = 0;
                this->adaptive_feasible_unstable_variable_count = 0;
                this->adaptive_covered_unstable_variable_count = 0;
                this->adaptive_candidate_count = 0;
                this->adaptive_valid_candidate_count = 0;
                this->adaptive_selected_view = -1;
                this->adaptive_fallback_residual_weight = 0;
                this->adaptive_best_candidate_cost = std::numeric_limits<double>::infinity();
                this->adaptive_selected_candidate_valid = false;
                this->adaptive_final_decoding.assign(this->bit_count, 0);
            }

            void resize_instability_state() {
                this->instability_flip_component.assign(this->bit_count, 0.0);
                this->instability_movement_component.assign(this->bit_count, 0.0);
                this->instability_message_conflict_component.assign(this->bit_count, 0.0);
                this->instability_residual_component.assign(this->bit_count, 0.0);
                this->instability_weak_belief_component.assign(this->bit_count, 0.0);
                this->instability_spatial_component.assign(this->bit_count, 0.0);
                this->instability_seed_scores.assign(this->bit_count, 0.0);
                this->previous_hard_decision.assign(this->bit_count, 0);
                this->previous_iteration_llr.assign(this->bit_count, 0.0);
                this->weighted_flip_sum.assign(this->bit_count, 0.0);
                this->weighted_llr_movement_sum.assign(this->bit_count, 0.0);
                this->weighted_message_conflict_sum.assign(this->bit_count, 0.0);
                this->weighted_residual_pressure_sum.assign(this->bit_count, 0.0);
                this->latest_check_residual.assign(this->check_count, 0);
                this->dynamics_transition_weight_sum = 0.0;
                this->dynamics_sample_weight_sum = 0.0;
                this->has_previous_dynamics_sample = false;
            }

            void reset_instability_dynamics() {
                this->resize_instability_state();
                this->sign_flip_count.assign(this->bit_count, 0);
            }

            static std::vector<double> percentile_ranks(const std::vector<double> &values) {
                std::vector<double> ranks(values.size(), 0.0);
                if (values.size() <= 1) return ranks;
                std::vector<double> sorted = values;
                std::sort(sorted.begin(), sorted.end());
                double denominator = double(values.size() - 1);
                for (size_t i = 0; i < values.size(); i++) {
                    auto first = std::lower_bound(sorted.begin(), sorted.end(), values[i]);
                    ranks[i] = double(std::distance(sorted.begin(), first)) / denominator;
                }
                return ranks;
            }

            static double bounded_signed_value(double value) {
                if (std::isnan(value)) return 0.0;
                if (std::isinf(value)) return std::signbit(value) ? -1.0 : 1.0;
                return std::tanh(0.5 * value);
            }

            void record_iteration_dynamics(
                    const std::vector<uint8_t> &syndrome,
                    int iteration) {
#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE == 1
                assert(syndrome.size() == size_t(this->check_count));
                assert(iteration > 0);
                double sample_weight = double(iteration);
                if (this->has_previous_dynamics_sample) {
                    this->dynamics_transition_weight_sum += sample_weight;
                }
                this->dynamics_sample_weight_sum += sample_weight;

                for (int row = 0; row < this->check_count; row++) {
                    this->latest_check_residual[row] =
                            this->candidate_syndrome[row] ^ syndrome[row];
                }

                for (int v = 0; v < this->bit_count; v++) {
                    if (this->has_previous_dynamics_sample) {
                        bool flipped = this->decoding[v] != this->previous_hard_decision[v];
                        this->weighted_flip_sum[v] += sample_weight * double(flipped);
                        double current_llr = this->log_prob_ratios[v];
                        double previous_llr = this->previous_iteration_llr[v];
                        double movement = std::isfinite(current_llr) && std::isfinite(previous_llr)
                                ? std::abs(current_llr - previous_llr)
                                : std::abs(
                                        bounded_signed_value(current_llr) -
                                        bounded_signed_value(previous_llr));
                        this->weighted_llr_movement_sum[v] += sample_weight * movement;
                        this->sign_flip_count[v] += int(flipped);
                    }

                    double message_sum = 0.0;
                    double message_abs_sum = 0.0;
                    bool messages_are_finite = true;
                    int residual_count = 0;
                    int degree = 0;
                    for (auto &entry: this->pcm.iterate_column(v)) {
                        message_sum += entry.check_to_bit_msg;
                        message_abs_sum += std::abs(entry.check_to_bit_msg);
                        messages_are_finite &= std::isfinite(entry.check_to_bit_msg);
                        residual_count += this->latest_check_residual[entry.row_index];
                        degree += 1;
                    }
                    double conflict_denominator = message_abs_sum;
                    double conflict_numerator = std::abs(message_sum);
                    if (!messages_are_finite) {
                        double bounded_message_sum = 0.0;
                        double bounded_message_abs_sum = 0.0;
                        for (auto &entry: this->pcm.iterate_column(v)) {
                            double bounded_message = bounded_signed_value(entry.check_to_bit_msg);
                            bounded_message_sum += bounded_message;
                            bounded_message_abs_sum += std::abs(bounded_message);
                        }
                        conflict_denominator = bounded_message_abs_sum;
                        conflict_numerator = std::abs(bounded_message_sum);
                    }
                    double conflict = conflict_denominator > 1e-12
                            ? 1.0 - conflict_numerator / conflict_denominator
                            : 0.0;
                    double residual_pressure = degree > 0
                            ? double(residual_count) / double(degree)
                            : 0.0;
                    this->weighted_message_conflict_sum[v] += sample_weight * conflict;
                    this->weighted_residual_pressure_sum[v] += sample_weight * residual_pressure;
                    this->previous_hard_decision[v] = this->decoding[v];
                    this->previous_iteration_llr[v] = this->log_prob_ratios[v];
                }
                this->has_previous_dynamics_sample = true;
#else
                (void) syndrome;
                (void) iteration;
#endif
            }

            void ensure_instability_scores() {
                if (!this->instability_scores_pending) return;
#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE == 0
                for (int v = 0; v < this->bit_count; v++) {
                    double llr_range = this->llr_max[v] - this->llr_min[v];
                    double weak_belief = 1.0 / (std::abs(this->log_prob_ratios[v]) + 1e-12);
                    this->instability_scores[v] =
                            (1.0 + double(this->sign_flip_count[v])) * (1.0 + llr_range) + weak_belief;
                    this->instability_seed_scores[v] = this->instability_scores[v];
                }
#else
                std::vector<double> raw_flip(this->bit_count, 0.0);
                std::vector<double> raw_movement(this->bit_count, 0.0);
                std::vector<double> raw_conflict(this->bit_count, 0.0);
                std::vector<double> raw_residual(this->bit_count, 0.0);
                std::vector<double> raw_weak(this->bit_count, 0.0);
                for (int v = 0; v < this->bit_count; v++) {
                    if (this->dynamics_transition_weight_sum > 0.0) {
                        raw_flip[v] = this->weighted_flip_sum[v] /
                                this->dynamics_transition_weight_sum;
                        raw_movement[v] = this->weighted_llr_movement_sum[v] /
                                this->dynamics_transition_weight_sum;
                    }
                    if (this->dynamics_sample_weight_sum > 0.0) {
                        raw_conflict[v] = this->weighted_message_conflict_sum[v] /
                                this->dynamics_sample_weight_sum;
                        raw_residual[v] = this->weighted_residual_pressure_sum[v] /
                                this->dynamics_sample_weight_sum;
                    }
                    double final_llr = this->log_prob_ratios[v];
                    raw_weak[v] = std::isnan(final_llr)
                            ? 1.0 / 1e-12
                            : 1.0 / (std::abs(final_llr) + 1e-12);
                }

                this->instability_flip_component = percentile_ranks(raw_flip);
                this->instability_movement_component = percentile_ranks(raw_movement);
                this->instability_message_conflict_component = percentile_ranks(raw_conflict);
                this->instability_residual_component = percentile_ranks(raw_residual);
                this->instability_weak_belief_component = percentile_ranks(raw_weak);
                for (int v = 0; v < this->bit_count; v++) {
                    double strongest_extrinsic_anomaly = std::max({
                            this->instability_message_conflict_component[v],
                            this->instability_residual_component[v]});
                    this->instability_seed_scores[v] =
                            (1.0 + raw_flip[v]) *
                            (1.0 + raw_movement[v]) *
                            (1.0 + strongest_extrinsic_anomaly) +
                            raw_weak[v];
                    assert(std::isfinite(this->instability_seed_scores[v]));
                }

                std::vector<double> raw_spatial(this->bit_count, 0.0);
                for (int v = 0; v < this->bit_count; v++) {
                    double corroboration_sum = 0.0;
                    int residual_neighbor_count = 0;
                    for (auto &column_entry: this->pcm.iterate_column(v)) {
                        int row = column_entry.row_index;
                        if (!this->latest_check_residual[row]) continue;
                        double neighbor_seed_sum = 0.0;
                        int neighbor_count = 0;
                        for (auto &row_entry: this->pcm.iterate_row(row)) {
                            if (row_entry.col_index == v) continue;
                            neighbor_seed_sum += this->instability_seed_scores[row_entry.col_index];
                            neighbor_count += 1;
                        }
                        if (neighbor_count > 0) {
                            corroboration_sum += neighbor_seed_sum / double(neighbor_count);
                            residual_neighbor_count += 1;
                        }
                    }
                    if (residual_neighbor_count > 0) {
                        raw_spatial[v] = corroboration_sum / double(residual_neighbor_count);
                    }
                }
                this->instability_spatial_component = percentile_ranks(raw_spatial);
                for (int v = 0; v < this->bit_count; v++) {
                    this->instability_scores[v] =
                            this->instability_seed_scores[v] *
                            (1.0 + this->instability_spatial_component[v]);
                    assert(std::isfinite(this->instability_scores[v]));
                }
#endif
                this->instability_scores_pending = false;
            }

            int adaptive_prepared_view_count() const {
                return int(this->adaptive_views.size());
            }

            int adaptive_candidate_pair_count() const {
                return int(this->adaptive_candidates.size());
            }

            int adaptive_selected_row_count(int view_index) const {
                return int(this->adaptive_views.at(view_index).candidate_indices.size());
            }

            int adaptive_weight_limit(int view_index) const {
                return this->adaptive_views.at(view_index).automatic_weight_limit;
            }

            int adaptive_selected_row_weight(int view_index, int row_index) const {
                int candidate_index = this->adaptive_views.at(view_index).candidate_indices.at(row_index);
                return int(this->adaptive_candidates[candidate_index].support.size());
            }

            void prepare_adaptive_redundant_views(
                    std::vector<uint8_t> &syndrome,
                    int adaptive_max_views) {
                if (adaptive_max_views < 0) {
                    throw std::runtime_error("Adaptive hyperparameters must be non-negative.");
                }
                if (syndrome.size() != size_t(this->check_count)) {
                    throw std::runtime_error("Syndrome length must equal the number of checks.");
                }

                this->reset_adaptive_multi_view_state();
                if (this->bit_count == 0 || adaptive_max_views == 0) {
                    this->instability_scores_pending = true;
                    return;
                }
                this->instability_scores_pending = true;
                this->ensure_instability_scores();

                std::vector<int> ranked_variables(this->bit_count);
                for (int v = 0; v < this->bit_count; v++) ranked_variables[v] = v;
                auto variable_is_more_unstable = [this](int a, int b) {
                    if (this->instability_scores[a] == this->instability_scores[b]) {
                        double abs_a = std::abs(this->log_prob_ratios[a]);
                        double abs_b = std::abs(this->log_prob_ratios[b]);
                        if (abs_a == abs_b) return a < b;
                        return abs_a < abs_b;
                    }
                    return this->instability_scores[a] > this->instability_scores[b];
                };
                std::sort(ranked_variables.begin(), ranked_variables.end(), variable_is_more_unstable);

                std::vector<double> sorted_scores = this->instability_scores;
                std::sort(sorted_scores.begin(), sorted_scores.end());
                double baseline = sorted_scores[this->bit_count / 2];
                if (this->bit_count % 2 == 0) {
                    baseline = 0.5 * (baseline + sorted_scores[this->bit_count / 2 - 1]);
                }

                double excess_mass = 0.0;
                for (double score: this->instability_scores) {
                    excess_mass += std::max(0.0, score - baseline);
                }
                int minimum_core_count = std::min(
                        this->bit_count,
                        std::max(4, int(std::ceil(std::log2(double(this->bit_count) + 1.0)))));
                int maximum_region_count = std::min(
                        this->bit_count,
                        std::max(minimum_core_count,
                                 int(std::ceil(2.0 * std::sqrt(double(this->bit_count))))));
                int core_count = minimum_core_count;
                if (excess_mass > 0.0) {
                    double cumulative_mass = 0.0;
                    for (int rank = 0; rank < maximum_region_count; rank++) {
                        cumulative_mass += std::max(
                                0.0, this->instability_scores[ranked_variables[rank]] - baseline);
                        if (rank + 1 >= minimum_core_count &&
                            cumulative_mass >= 0.8 * excess_mass) {
                            core_count = rank + 1;
                            break;
                        }
                        core_count = rank + 1;
                    }
                }
                core_count = std::min(core_count, maximum_region_count);
                this->adaptive_core_variable_count = core_count;

                std::vector<int> variable_order(
                        ranked_variables.begin(), ranked_variables.begin() + core_count);
                std::vector<uint8_t> selected_variable_flag(this->bit_count, 0);
                std::vector<uint8_t> core_neighbor_check_flag(this->check_count, 0);
                for (int variable: variable_order) {
                    selected_variable_flag[variable] = 1;
                    for (auto &entry: this->pcm.iterate_column(variable)) {
                        core_neighbor_check_flag[entry.row_index] = 1;
                    }
                }

                double boundary_score = this->instability_scores[ranked_variables[core_count - 1]];
                double expansion_threshold = baseline + 0.5 * std::max(0.0, boundary_score - baseline);
                std::vector<int> frontier;
                for (int variable: ranked_variables) {
                    if (selected_variable_flag[variable] ||
                        this->instability_scores[variable] <= expansion_threshold) {
                        continue;
                    }
                    for (auto &entry: this->pcm.iterate_column(variable)) {
                        if (core_neighbor_check_flag[entry.row_index]) {
                            frontier.push_back(variable);
                            break;
                        }
                    }
                }
                for (int variable: frontier) {
                    if (int(variable_order.size()) >= maximum_region_count) break;
                    selected_variable_flag[variable] = 1;
                    variable_order.push_back(variable);
                }

                int unstable_count = int(variable_order.size());
                this->adaptive_unstable_variable_count = unstable_count;
                this->adaptive_region_variables = variable_order;

                std::vector<int> unstable_rank(this->bit_count, -1);
                std::vector<uint8_t> neighboring_check_flag(this->check_count, 0);
                std::vector<int> neighboring_checks;
                for (int idx = 0; idx < unstable_count; idx++) {
                    int v = variable_order[idx];
                    unstable_rank[v] = idx;
                    for (auto &e: this->pcm.iterate_column(v)) {
                        if (!neighboring_check_flag[e.row_index]) {
                            neighboring_check_flag[e.row_index] = 1;
                            neighboring_checks.push_back(e.row_index);
                        }
                    }
                }
                std::sort(neighboring_checks.begin(), neighboring_checks.end());

                std::vector<uint8_t> residual(this->check_count, 0);
                for (int row = 0; row < this->check_count; row++) {
                    residual[row] = this->candidate_syndrome[row] ^ syndrome[row];
                }

                {
                    std::vector<uint8_t> row_marks(this->bit_count, 0);
                    std::vector<int> touched_columns;
                    std::set<std::vector<int>> seen_supports;
                    for (int ia = 0; ia < int(neighboring_checks.size()); ia++) {
                        for (int ib = ia + 1; ib < int(neighboring_checks.size()); ib++) {
                            int row_a = neighboring_checks[ia];
                            int row_b = neighboring_checks[ib];
                            touched_columns.clear();

                            for (auto &e: this->pcm.iterate_row(row_a)) {
                                if (!row_marks[e.col_index]) touched_columns.push_back(e.col_index);
                                row_marks[e.col_index] ^= 1;
                            }
                            for (auto &e: this->pcm.iterate_row(row_b)) {
                                if (!row_marks[e.col_index]) touched_columns.push_back(e.col_index);
                                row_marks[e.col_index] ^= 1;
                            }

                            std::vector<int> support;
                            for (int col: touched_columns) {
                                if (row_marks[col]) support.push_back(col);
                                row_marks[col] = 0;
                            }
                            std::sort(support.begin(), support.end());
                            if (support.empty()) continue;
                            if (!seen_supports.insert(support).second) continue;

                            std::vector<int> unstable_support;
                            for (int col: support) {
                                if (unstable_rank[col] != -1) {
                                    unstable_support.push_back(unstable_rank[col]);
                                }
                            }
                            if (unstable_support.empty()) continue;

                            this->adaptive_candidates.push_back({
                                    row_a,
                                    row_b,
                                    bool(residual[row_a] ^ residual[row_b]),
                                    std::move(support),
                                    std::move(unstable_support)});
                        }
                    }
                }

                std::vector<uint8_t> feasible_unstable_flag(unstable_count, 0);
                for (const auto &candidate: this->adaptive_candidates) {
                    for (int rank: candidate.unstable_support) feasible_unstable_flag[rank] = 1;
                }
                for (int rank = 0; rank < unstable_count; rank++) {
                    this->adaptive_feasible_unstable_variable_count += feasible_unstable_flag[rank];
                }

                std::vector<uint8_t> candidate_used(this->adaptive_candidates.size(), 0);
                std::vector<int> global_coverage(unstable_count, 0);
                for (int view = 0; view < adaptive_max_views; view++) {
                    std::vector<int> minimum_unused_weight(unstable_count, std::numeric_limits<int>::max());
                    for (int candidate_index = 0;
                         candidate_index < int(this->adaptive_candidates.size());
                         candidate_index++) {
                        if (candidate_used[candidate_index]) continue;
                        const auto &candidate = this->adaptive_candidates[candidate_index];
                        int weight = int(candidate.support.size());
                        for (int rank: candidate.unstable_support) {
                            minimum_unused_weight[rank] = std::min(minimum_unused_weight[rank], weight);
                        }
                    }

                    int automatic_weight_limit = 0;
                    int uncovered_variable_count = 0;
                    for (int rank = 0; rank < unstable_count; rank++) {
                        if (minimum_unused_weight[rank] != std::numeric_limits<int>::max()) {
                            automatic_weight_limit = std::max(
                                    automatic_weight_limit, minimum_unused_weight[rank]);
                            uncovered_variable_count += 1;
                        }
                    }
                    if (automatic_weight_limit == 0) break;

                    std::vector<uint8_t> local_coverage(unstable_count, 0);
                    std::vector<int> selected_indices;

                    while (uncovered_variable_count > 0) {
                        int best_index = -1;
                        int best_new_count = -1;
                        double best_new_instability = -1.0;
                        double best_balanced_instability = -1.0;
                        int best_residual_relevance = -1;
#if LDPC_MODIFIED_ADAPTIVE_BP_USE_MAX_OVERLAP
                        int best_max_overlap = std::numeric_limits<int>::max();
#endif
                        int best_weight = std::numeric_limits<int>::max();

                        for (int candidate_index = 0;
                             candidate_index < int(this->adaptive_candidates.size());
                             candidate_index++) {
                            if (candidate_used[candidate_index]) continue;
                            const auto &candidate = this->adaptive_candidates[candidate_index];
                            if (int(candidate.support.size()) > automatic_weight_limit) continue;

                            int new_count = 0;
                            double new_instability = 0.0;
                            double balanced_instability = 0.0;
                            for (int rank: candidate.unstable_support) {
                                int v = variable_order[rank];
                                if (!local_coverage[rank]) {
                                    new_count += 1;
                                    new_instability += this->instability_scores[v];
                                }
                                balanced_instability +=
                                        this->instability_scores[v] / (1.0 + double(global_coverage[rank]));
                            }
                            if (new_count == 0) continue;

                            int residual_relevance = candidate.residual_relevant ? 1 : 0;
                            int weight = int(candidate.support.size());
                            int prefix_order = 0;
                            if (best_index == -1) {
                                prefix_order = 1;
                            } else if (new_count != best_new_count) {
                                prefix_order = new_count > best_new_count ? 1 : -1;
                            } else if (new_instability != best_new_instability) {
                                prefix_order = new_instability > best_new_instability ? 1 : -1;
                            } else if (balanced_instability != best_balanced_instability) {
                                prefix_order = balanced_instability > best_balanced_instability ? 1 : -1;
                            } else if (residual_relevance != best_residual_relevance) {
                                prefix_order = residual_relevance > best_residual_relevance ? 1 : -1;
                            }
                            if (prefix_order < 0) continue;

#if LDPC_MODIFIED_ADAPTIVE_BP_USE_MAX_OVERLAP
                            int max_overlap = 0;
                            for (int selected_index: selected_indices) {
                                const auto &selected_support =
                                        this->adaptive_candidates[selected_index].support;
                                int overlap = 0;
                                size_t left = 0;
                                size_t right = 0;
                                while (left < candidate.support.size() && right < selected_support.size()) {
                                    if (candidate.support[left] == selected_support[right]) {
                                        overlap += 1;
                                        left += 1;
                                        right += 1;
                                    } else if (candidate.support[left] < selected_support[right]) {
                                        left += 1;
                                    } else {
                                        right += 1;
                                    }
                                }
                                max_overlap = std::max(max_overlap, overlap);
                            }
#endif

                            bool better = prefix_order > 0;
#if LDPC_MODIFIED_ADAPTIVE_BP_USE_MAX_OVERLAP
                            if (prefix_order == 0 && max_overlap != best_max_overlap) {
                                better = max_overlap < best_max_overlap;
                            } else if (prefix_order == 0 && weight != best_weight) {
#else
                            if (prefix_order == 0 && weight != best_weight) {
#endif
                                better = weight < best_weight;
                            }

                            if (better) {
                                best_index = candidate_index;
                                best_new_count = new_count;
                                best_new_instability = new_instability;
                                best_balanced_instability = balanced_instability;
                                best_residual_relevance = residual_relevance;
#if LDPC_MODIFIED_ADAPTIVE_BP_USE_MAX_OVERLAP
                                best_max_overlap = max_overlap;
#endif
                                best_weight = weight;
                            }
                        }

                        if (best_index == -1) break;
                        candidate_used[best_index] = 1;
                        selected_indices.push_back(best_index);
                        for (int rank: this->adaptive_candidates[best_index].unstable_support) {
                            if (!local_coverage[rank]) {
                                local_coverage[rank] = 1;
                                uncovered_variable_count -= 1;
                            }
                            global_coverage[rank] += 1;
                        }
                    }

                    if (selected_indices.empty()) break;
                    assert(int(selected_indices.size()) <= unstable_count);
                    this->adaptive_views.push_back({std::move(selected_indices), automatic_weight_limit});
                }

                for (int rank = 0; rank < unstable_count; rank++) {
                    this->adaptive_covered_unstable_variable_count += global_coverage[rank] > 0;
                }
            }

            void build_prepared_adaptive_redundant_view(
                    int view_index,
                    std::vector<uint8_t> &syndrome,
                    BpSparse &augmented_pcm,
                    std::vector<uint8_t> &augmented_syndrome) {
                if (view_index < 0 || view_index >= int(this->adaptive_views.size())) {
                    throw std::runtime_error("Adaptive view index is out of range.");
                }

                const auto &selected_indices = this->adaptive_views[view_index].candidate_indices;
                this->adaptive_executed_view_count = std::max(
                        this->adaptive_executed_view_count, view_index + 1);

                int augmented_entry_count = this->pcm.entry_count();
                for (int candidate_index: selected_indices) {
                    const auto &candidate = this->adaptive_candidates[candidate_index];
                    augmented_entry_count += int(candidate.support.size());
                }

                augmented_pcm.allocate(
                        this->check_count + int(selected_indices.size()),
                        this->bit_count,
                        augmented_entry_count);
                augmented_syndrome.assign(
                        this->check_count + int(selected_indices.size()), 0);

                for (int row = 0; row < this->check_count; row++) {
                    augmented_syndrome[row] = syndrome[row];
                    for (auto &e: this->pcm.iterate_row(row)) {
                        augmented_pcm.insert_entry(row, e.col_index);
                    }
                }

                for (int i = 0; i < int(selected_indices.size()); i++) {
                    const auto &candidate = this->adaptive_candidates[selected_indices[i]];
                    int row = this->check_count + i;
                    augmented_syndrome[row] = syndrome[candidate.row_a] ^ syndrome[candidate.row_b];
                    for (int col: candidate.support) augmented_pcm.insert_entry(row, col);
                }
            }

            void consider_adaptive_candidate(
                    std::vector<uint8_t> &syndrome,
                    const std::vector<uint8_t> &candidate,
                    const std::vector<uint8_t> &candidate_check_values,
                    int source_view) {
                if (syndrome.size() != size_t(this->check_count) ||
                    candidate.size() != size_t(this->bit_count) ||
                    candidate_check_values.size() < size_t(this->check_count)) {
                    throw std::runtime_error("Adaptive candidate dimensions do not match the original PCM.");
                }

                int residual_weight = 0;
                for (int row = 0; row < this->check_count; row++) {
                    residual_weight += candidate_check_values[row] != syndrome[row];
                }

                this->adaptive_candidate_count += 1;
                if (source_view == 0) {
                    this->adaptive_final_decoding = candidate;
                    this->adaptive_selected_view = 0;
                    this->adaptive_fallback_residual_weight = residual_weight;
                }

                if (residual_weight != 0) return;
                this->adaptive_valid_candidate_count += 1;

                double candidate_cost = 0.0;
                for (int v = 0; v < this->bit_count; v++) {
                    double p = this->channel_probabilities[v];
                    assert(p > 0.0 && p < 1.0);
                    candidate_cost += candidate[v] ? -std::log(p) : -std::log(1.0 - p);
                }
                if (!this->adaptive_selected_candidate_valid ||
                    candidate_cost < this->adaptive_best_candidate_cost) {
                    this->adaptive_selected_candidate_valid = true;
                    this->adaptive_best_candidate_cost = candidate_cost;
                    this->adaptive_selected_view = source_view;
                    this->adaptive_final_decoding = candidate;
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
                } else if (schedule == SERIAL || schedule == SERIAL_RELATIVE) {
                    return bp_decode_serial(input_vector);
                } else {
                    throw std::runtime_error("Invalid BP schedule");
                }
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
#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE == 0
                            if (this->decoding[i] == 0) this->sign_flip_count[i] += 1;
#endif
                            this->decoding[i] = 1;
                            for (auto &e: this->pcm.iterate_column(i)) {
                                this->candidate_syndrome[e.row_index] ^= 1;
                            }
                        } else {
#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE == 0
                            if (this->decoding[i] == 1) this->sign_flip_count[i] += 1;
#endif
                            this->decoding[i] = 0;
                        }
                    }

                    this->record_iteration_dynamics(syndrome, it);

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
#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE == 0
                            if (this->decoding[bit_index] == 0) this->sign_flip_count[bit_index] += 1;
#endif
                            this->decoding[bit_index] = 1;
                        } else {
#if LDPC_MODIFIED_EXPANDED_BP_INSTABILITY_MODE == 0
                            if (this->decoding[bit_index] == 1) this->sign_flip_count[bit_index] += 1;
#endif
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
                    this->record_iteration_dynamics(syndrome, it);
                    this->iterations = it;
                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        this->converge = true;
                        return this->decoding;
                    }
                }
                return this->decoding;
            }

        };

        class AdaptiveViewDecoderPool {
        public:
            AdaptiveViewDecoderPool() = default;
            AdaptiveViewDecoderPool(const AdaptiveViewDecoderPool &) = delete;
            AdaptiveViewDecoderPool &operator=(const AdaptiveViewDecoderPool &) = delete;

            void run(
                    ExpandedAdaptiveBpDecoder &main_decoder,
                    std::vector<uint8_t> &syndrome,
                    int view_thread_count) {
                if (view_thread_count < 0) {
                    throw std::runtime_error("view_thread_count must be non-negative.");
                }
                active_view_count_ = main_decoder.adaptive_prepared_view_count();
                ensure_slots(main_decoder, active_view_count_);
                const std::vector<int> &serial_schedule = main_decoder.random_serial_schedule
                        ? NULL_INT_VECTOR
                        : main_decoder.serial_schedule_order;
                for (int view_index = 0; view_index < active_view_count_; view_index++) {
                    main_decoder.build_prepared_adaptive_redundant_view(
                            view_index,
                            syndrome,
                            *view_pcms_[view_index],
                            view_syndromes_[view_index]);
                    view_decoders_[view_index]->update_decoder_parameters(
                            main_decoder.channel_probabilities,
                            main_decoder.maximum_iterations,
                            main_decoder.bp_method,
                            main_decoder.schedule,
                            main_decoder.ms_scaling_factor,
                            1,
                            serial_schedule,
                            main_decoder.random_schedule_seed,
                            main_decoder.random_serial_schedule,
                            SYNDROME);
                }

                last_worker_count_ = resolve_worker_count(
                        active_view_count_, view_thread_count);
                if (last_worker_count_ <= 1) {
                    for (int view_index = 0; view_index < active_view_count_; view_index++) {
                        view_decoders_[view_index]->decode(view_syndromes_[view_index]);
                    }
                    return;
                }

                std::atomic<int> next_view{0};
                std::mutex error_mutex;
                std::exception_ptr worker_error;
                auto worker = [&]() {
                    while (true) {
                        int view_index = next_view.fetch_add(1, std::memory_order_relaxed);
                        if (view_index >= active_view_count_) return;
                        try {
                            view_decoders_[view_index]->decode(view_syndromes_[view_index]);
                        } catch (...) {
                            std::lock_guard<std::mutex> lock(error_mutex);
                            if (!worker_error) worker_error = std::current_exception();
                        }
                    }
                };
                {
                    std::vector<std::jthread> workers;
                    workers.reserve(last_worker_count_);
                    for (int i = 0; i < last_worker_count_; i++) {
                        workers.emplace_back(worker);
                    }
                }
                if (worker_error) std::rethrow_exception(worker_error);
            }

            void reset() {
                active_view_count_ = 0;
                last_worker_count_ = 0;
            }

            int active_view_count() const {
                return active_view_count_;
            }

            int last_worker_count() const {
                return last_worker_count_;
            }

            ExpandedAdaptiveBpDecoder &decoder(int view_index) {
                return *view_decoders_.at(view_index);
            }

            BpSparse &pcm(int view_index) {
                return *view_pcms_.at(view_index);
            }

        private:
            std::vector<std::unique_ptr<BpSparse>> view_pcms_;
            std::vector<std::unique_ptr<ExpandedAdaptiveBpDecoder>> view_decoders_;
            std::vector<std::vector<uint8_t>> view_syndromes_;
            int active_view_count_ = 0;
            int last_worker_count_ = 0;

            static int resolve_worker_count(int view_count, int configured_count) {
                if (view_count == 0) return 0;
                unsigned int hardware_threads = std::thread::hardware_concurrency();
                int available = hardware_threads == 0 ? 1 : int(hardware_threads);
                int requested = configured_count == 0 ? available : configured_count;
                return std::max(1, std::min(view_count, requested));
            }

            void ensure_slots(ExpandedAdaptiveBpDecoder &main_decoder, int count) {
                while (int(view_decoders_.size()) < count) {
                    auto pcm = std::make_unique<BpSparse>(
                            main_decoder.check_count, main_decoder.bit_count, 0);
                    BpSparse *pcm_pointer = pcm.get();
                    view_pcms_.push_back(std::move(pcm));
                    view_decoders_.push_back(std::make_unique<ExpandedAdaptiveBpDecoder>(
                            *pcm_pointer,
                            main_decoder.channel_probabilities,
                            main_decoder.maximum_iterations,
                            main_decoder.bp_method,
                            main_decoder.schedule,
                            main_decoder.ms_scaling_factor,
                            1,
                            NULL_INT_VECTOR,
                            main_decoder.random_schedule_seed,
                            main_decoder.random_serial_schedule,
                            SYNDROME));
                    view_syndromes_.emplace_back();
                }
            }
        };
    }
}  // namespace ldpc_modified::expanded_adaptive_bp

#endif
