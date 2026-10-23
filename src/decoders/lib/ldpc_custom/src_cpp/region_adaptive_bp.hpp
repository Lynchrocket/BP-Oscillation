#ifndef REGION_ADAPTIVE_BP_H
#define REGION_ADAPTIVE_BP_H

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
#include <cstdint>
#include <thread>
#include <mutex>
#include <exception>

#include "math.h"
#include "sparse_matrix_base.hpp"
#include "gf2sparse.hpp"
#include "rng.hpp"

#ifndef LDPC_MODIFIED_REGION_BP_SELECTION_MODE
#define LDPC_MODIFIED_REGION_BP_SELECTION_MODE 0
#endif

#ifndef LDPC_MODIFIED_REGION_BP_RANDOM_SEED
#define LDPC_MODIFIED_REGION_BP_RANDOM_SEED 0
#endif

namespace ldpc_modified {
    namespace region_adaptive_bp {

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

        inline std::atomic<int> &region_adaptive_bp_live_instance_count_ref() {
            static std::atomic<int> count{0};
            return count;
        }

        inline int region_adaptive_bp_live_instance_count() {
            return region_adaptive_bp_live_instance_count_ref().load(std::memory_order_relaxed);
        }

        class BpEntry : public ldpc_modified::sparse_matrix_base::EntryBase<BpEntry> {
        public:
            double bit_to_check_msg = 0.0;
            double check_to_bit_msg = 0.0;

            ~BpEntry() = default;
        };
        using BpSparse = ldpc_modified::gf2sparse::GF2Sparse<BpEntry>;

        class RegionAdaptiveBpDecoder {
            // TODO properties should be private and only accessible via getters and setters
        public:
            struct RedundantRowCandidate {
                int row_a;
                int row_b;
                std::vector<int> support;
                int primary_target;
                double actionability;
                double utility;
                double cycle_density;
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
            bool instability_scores_pending;
            std::vector<int> check_residual_count;
            std::vector<int> check_residual_flip_count;
            std::vector<uint8_t> previous_check_residual;
            int trajectory_sample_count;

            std::vector<int> region_core_variables;
            std::vector<int> region_conflict_checks;
            std::vector<int> region_anchor_variables;

            std::vector<RedundantRowCandidate> adaptive_candidates;
            std::vector<PreparedAdaptiveView> adaptive_views;
            std::vector<uint8_t> adaptive_final_decoding;
            int adaptive_executed_view_count;
            int adaptive_main_bp_run_count;
            int adaptive_unstable_variable_count;
            int adaptive_feasible_unstable_variable_count;
            int adaptive_covered_unstable_variable_count;
            int adaptive_candidate_count;
            int adaptive_valid_candidate_count;
            int adaptive_selected_view;
            int adaptive_fallback_residual_weight;
            double adaptive_best_candidate_cost;
            bool adaptive_selected_candidate_valid;

            RegionAdaptiveBpDecoder(
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
                this->instability_scores_pending = false;
                this->check_residual_count.resize(check_count);
                this->check_residual_flip_count.resize(check_count);
                this->previous_check_residual.resize(check_count);
                this->trajectory_sample_count = 0;
                this->adaptive_final_decoding.assign(bit_count, 0);
                this->adaptive_executed_view_count = 0;
                this->adaptive_main_bp_run_count = 0;
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
                region_adaptive_bp_live_instance_count_ref().fetch_add(1, std::memory_order_relaxed);
            }

            ~RegionAdaptiveBpDecoder() {
                region_adaptive_bp_live_instance_count_ref().fetch_sub(1, std::memory_order_relaxed);
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
                this->instability_scores_pending = false;
                this->adaptive_candidates.clear();
                this->adaptive_views.clear();
                this->adaptive_final_decoding.assign(this->bit_count, 0);
                this->adaptive_executed_view_count = 0;
                this->adaptive_main_bp_run_count = 0;
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
                this->instability_scores_pending = true;
                this->check_residual_count.assign(this->check_count, 0);
                this->check_residual_flip_count.assign(this->check_count, 0);
                this->previous_check_residual.assign(this->check_count, 0);
                this->trajectory_sample_count = 0;
            }

            void reset_zero_syndrome_state() {
                this->converge = true;
                this->iterations = 0;
                this->candidate_syndrome.assign(this->check_count, 0);
                this->decoding.assign(this->bit_count, 0);
                this->sign_flip_count.assign(this->bit_count, 0);
                this->check_residual_count.assign(this->check_count, 0);
                this->check_residual_flip_count.assign(this->check_count, 0);
                this->previous_check_residual.assign(this->check_count, 0);
                this->trajectory_sample_count = 0;

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
                this->adaptive_executed_view_count = 0;
                this->adaptive_main_bp_run_count = 0;
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
                this->region_core_variables.clear();
                this->region_conflict_checks.clear();
                this->region_anchor_variables.clear();
            }

            void ensure_instability_scores() {
                if (!this->instability_scores_pending) return;
                std::vector<double> oscillation_scores(this->bit_count, 0.0);
                std::vector<double> excursion_scores(this->bit_count, 0.0);
                std::vector<double> weak_belief_scores(this->bit_count, 0.0);
                for (int v = 0; v < this->bit_count; v++) {
                    double llr_range = this->llr_max[v] - this->llr_min[v];
                    double prior_scale = std::abs(this->initial_log_prob_ratios[v]);
                    oscillation_scores[v] = double(this->sign_flip_count[v]) /
                            (1.0 + double(this->sign_flip_count[v]));
                    excursion_scores[v] = llr_range /
                            (llr_range + prior_scale + 1e-12);
                    weak_belief_scores[v] = prior_scale /
                            (prior_scale + std::abs(this->log_prob_ratios[v]) + 1e-12);
                }
#if LDPC_MODIFIED_REGION_BP_SELECTION_MODE == 1
                this->instability_scores = weak_belief_scores;
#else
                double oscillation_median = robust_median(oscillation_scores);
                double excursion_median = robust_median(excursion_scores);
                std::vector<double> oscillation_deviations(this->bit_count, 0.0);
                std::vector<double> excursion_deviations(this->bit_count, 0.0);
                for (int v = 0; v < this->bit_count; v++) {
                    oscillation_deviations[v] =
                            std::abs(oscillation_scores[v] - oscillation_median);
                    excursion_deviations[v] =
                            std::abs(excursion_scores[v] - excursion_median);
                }
                double oscillation_scale = robust_median(oscillation_deviations) + 1e-12;
                double excursion_scale = robust_median(excursion_deviations) + 1e-12;
                for (int v = 0; v < this->bit_count; v++) {
                    double oscillation_excess = std::max(
                            0.0, oscillation_scores[v] - oscillation_median);
                    double excursion_excess = std::max(
                            0.0, excursion_scores[v] - excursion_median);
                    double normalized_oscillation = oscillation_excess /
                            (oscillation_excess + oscillation_scale);
                    double normalized_excursion = excursion_excess /
                            (excursion_excess + excursion_scale);
                    double anomalous_dynamics = 1.0 -
                            (1.0 - normalized_oscillation) *
                            (1.0 - normalized_excursion);
                    this->instability_scores[v] = 1.0 -
                            (1.0 - weak_belief_scores[v]) *
                            (1.0 - anomalous_dynamics);
                }
#endif
                this->instability_scores_pending = false;
            }

            void observe_check_residuals(const std::vector<uint8_t> &syndrome) {
                for (int row = 0; row < this->check_count; row++) {
                    uint8_t residual = this->candidate_syndrome[row] ^ syndrome[row];
                    this->check_residual_count[row] += residual;
                    if (this->trajectory_sample_count > 0 &&
                        residual != this->previous_check_residual[row]) {
                        this->check_residual_flip_count[row] += 1;
                    }
                    this->previous_check_residual[row] = residual;
                }
                this->trajectory_sample_count += 1;
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

            static double robust_median(std::vector<double> values) {
                assert(!values.empty());
                size_t middle = values.size() / 2;
                std::nth_element(values.begin(), values.begin() + middle, values.end());
                double upper = values[middle];
                if (values.size() % 2 != 0) return upper;
                std::nth_element(values.begin(), values.begin() + middle - 1, values.end());
                return 0.5 * (upper + values[middle - 1]);
            }

            static int sorted_support_overlap(
                    const std::vector<int> &left_support,
                    const std::vector<int> &right_support) {
                int overlap = 0;
                size_t left = 0;
                size_t right = 0;
                while (left < left_support.size() && right < right_support.size()) {
                    if (left_support[left] == right_support[right]) {
                        overlap += 1;
                        left += 1;
                        right += 1;
                    } else if (left_support[left] < right_support[right]) {
                        left += 1;
                    } else {
                        right += 1;
                    }
                }
                return overlap;
            }

            double predicted_redundant_message(
                    const std::vector<int> &support,
                    int target,
                    uint8_t redundant_syndrome) const {
                int negative_count = redundant_syndrome;
                double minimum_magnitude = std::numeric_limits<double>::max();
                double product = 1.0;
                for (int variable: support) {
                    if (variable == target) continue;
                    double llr = this->log_prob_ratios[variable];
                    negative_count += llr <= 0.0;
                    minimum_magnitude = std::min(minimum_magnitude, std::abs(llr));
                    product *= std::tanh(llr / 2.0);
                }
                double sign = negative_count % 2 == 0 ? 1.0 : -1.0;
                if (this->bp_method == MINIMUM_SUM) {
                    double alpha = this->ms_scaling_factor == 0.0
                            ? 1.0 - std::pow(2.0, -double(std::max(1, this->iterations)))
                            : this->ms_scaling_factor;
                    return sign * alpha * minimum_magnitude;
                }
                product = std::max(-1.0 + 1e-15, std::min(1.0 - 1e-15, product));
                return 2.0 * std::atanh(product) * (redundant_syndrome ? -1.0 : 1.0);
            }

            void prepare_region_redundant_views(
                    std::vector<uint8_t> &syndrome,
                    int adaptive_max_views) {
                assert(adaptive_max_views >= 0);
                assert(syndrome.size() == size_t(this->check_count));
                this->reset_adaptive_multi_view_state();
                if (adaptive_max_views == 0 || this->bit_count == 0) return;

                this->instability_scores_pending = true;
                this->ensure_instability_scores();

                std::vector<double> conflict_scores(this->check_count, 0.0);
                for (int row = 0; row < this->check_count; row++) {
                    double activity = this->trajectory_sample_count == 0 ? 0.0 :
                            double(this->check_residual_count[row]) /
                            double(this->trajectory_sample_count);
                    double flip_rate = this->trajectory_sample_count <= 1 ? 0.0 :
                            double(this->check_residual_flip_count[row]) /
                            double(this->trajectory_sample_count - 1);
                    double final_residual = this->candidate_syndrome[row] ^ syndrome[row];
                    conflict_scores[row] = 1.0 -
                            (1.0 - final_residual) * (1.0 - activity) * (1.0 - flip_rate);
                }
                double conflict_median = robust_median(conflict_scores);
                std::vector<double> conflict_deviations(this->check_count, 0.0);
                for (int row = 0; row < this->check_count; row++) {
                    conflict_deviations[row] = std::abs(conflict_scores[row] - conflict_median);
                }
                double conflict_threshold = conflict_median + robust_median(conflict_deviations);
                std::vector<uint8_t> conflict_flag(this->check_count, 0);
                for (int row = 0; row < this->check_count; row++) {
                    bool final_residual = this->candidate_syndrome[row] ^ syndrome[row];
                    if (final_residual || conflict_scores[row] > conflict_threshold) {
                        conflict_flag[row] = 1;
                        this->region_conflict_checks.push_back(row);
                    }
                }

                double instability_median = robust_median(this->instability_scores);
                std::vector<double> instability_deviations(this->bit_count, 0.0);
                for (int variable = 0; variable < this->bit_count; variable++) {
                    instability_deviations[variable] =
                            std::abs(this->instability_scores[variable] - instability_median);
                }
                double instability_threshold =
                        instability_median + robust_median(instability_deviations);

                std::vector<double> variable_relevance(this->bit_count, 0.0);
                std::vector<uint8_t> core_candidate_flag(this->bit_count, 0);
                for (int variable = 0; variable < this->bit_count; variable++) {
                    double neighboring_conflict = 0.0;
                    for (auto &entry: this->pcm.iterate_column(variable)) {
                        neighboring_conflict = std::max(
                                neighboring_conflict, conflict_scores[entry.row_index]);
                        core_candidate_flag[variable] |= conflict_flag[entry.row_index];
                    }
                    if (this->instability_scores[variable] > instability_threshold) {
                        core_candidate_flag[variable] = 1;
                    }
                    variable_relevance[variable] =
                            this->instability_scores[variable] + neighboring_conflict;
                }

                std::vector<int> variable_order;
                for (int variable = 0; variable < this->bit_count; variable++) {
                    if (core_candidate_flag[variable]) variable_order.push_back(variable);
                }
                if (variable_order.empty()) {
                    variable_order.resize(this->bit_count);
                    for (int variable = 0; variable < this->bit_count; variable++) {
                        variable_order[variable] = variable;
                    }
                }
                std::sort(variable_order.begin(), variable_order.end(),
                          [&variable_relevance](int left, int right) {
                    if (variable_relevance[left] != variable_relevance[right]) {
                        return variable_relevance[left] > variable_relevance[right];
                    }
                    return left < right;
                });
                int automatic_core_limit = std::min(
                        this->bit_count,
                        std::max(4, int(std::ceil(std::sqrt(double(this->bit_count))))));
                if (int(variable_order.size()) > automatic_core_limit) {
                    variable_order.resize(automatic_core_limit);
                }
                this->region_core_variables = variable_order;
                this->adaptive_unstable_variable_count = int(variable_order.size());

                std::vector<uint8_t> core_flag(this->bit_count, 0);
                std::vector<uint8_t> region_check_flag(this->check_count, 0);
                std::vector<int> region_checks;
                for (int variable: this->region_core_variables) {
                    core_flag[variable] = 1;
                    for (auto &entry: this->pcm.iterate_column(variable)) {
                        if (!region_check_flag[entry.row_index]) {
                            region_check_flag[entry.row_index] = 1;
                            region_checks.push_back(entry.row_index);
                        }
                    }
                }
                for (int row: this->region_conflict_checks) {
                    if (!region_check_flag[row]) {
                        region_check_flag[row] = 1;
                        region_checks.push_back(row);
                    }
                }
                std::sort(region_checks.begin(), region_checks.end());

                std::vector<uint8_t> anchor_flag(this->bit_count, 0);
                for (int row: region_checks) {
                    for (auto &entry: this->pcm.iterate_row(row)) {
                        if (!core_flag[entry.col_index] && !anchor_flag[entry.col_index]) {
                            anchor_flag[entry.col_index] = 1;
                            this->region_anchor_variables.push_back(entry.col_index);
                        }
                    }
                }
                std::sort(this->region_anchor_variables.begin(), this->region_anchor_variables.end());

                std::vector<uint8_t> row_marks(this->bit_count, 0);
                std::vector<int> touched_columns;
                std::set<std::vector<int>> seen_supports;
                for (int left_index = 0; left_index < int(region_checks.size()); left_index++) {
                    for (int right_index = left_index + 1;
                         right_index < int(region_checks.size()); right_index++) {
                        int row_a = region_checks[left_index];
                        int row_b = region_checks[right_index];
                        touched_columns.clear();
                        for (auto &entry: this->pcm.iterate_row(row_a)) {
                            if (!row_marks[entry.col_index]) touched_columns.push_back(entry.col_index);
                            row_marks[entry.col_index] ^= 1;
                        }
                        for (auto &entry: this->pcm.iterate_row(row_b)) {
                            if (!row_marks[entry.col_index]) touched_columns.push_back(entry.col_index);
                            row_marks[entry.col_index] ^= 1;
                        }
                        std::vector<int> support;
                        for (int variable: touched_columns) {
                            if (row_marks[variable]) support.push_back(variable);
                            row_marks[variable] = 0;
                        }
                        std::sort(support.begin(), support.end());
                        if (support.size() < 2 || !seen_supports.insert(support).second) continue;

                        std::vector<int> targets;
                        for (int variable: support) {
                            if (core_flag[variable]) targets.push_back(variable);
                        }
                        if (targets.empty()) continue;

                        uint8_t redundant_syndrome = syndrome[row_a] ^ syndrome[row_b];
                        int primary_target = -1;
                        double best_actionability = 0.0;
                        for (int target: targets) {
                            double anchor_quality = 1.0;
                            int ambiguous_core_count = 0;
                            for (int variable: support) {
                                if (variable == target) continue;
                                ambiguous_core_count += core_flag[variable];
                                double magnitude = std::abs(this->log_prob_ratios[variable]);
                                double reliability = (1.0 - this->instability_scores[variable]) *
                                        magnitude /
                                        (magnitude + std::abs(this->initial_log_prob_ratios[variable]) + 1e-12);
                                anchor_quality = std::min(anchor_quality, reliability);
                            }
                            double message = this->predicted_redundant_message(
                                    support, target, redundant_syndrome);
                            double target_llr = this->log_prob_ratios[target];
                            bool corrective_direction = target_llr == 0.0 ||
                                    std::copysign(1.0, target_llr) * message < 0.0;
                            double intervention = std::abs(message) *
                                    (corrective_direction ? 1.0 : 0.25);
                            double actionability = this->instability_scores[target] *
                                    intervention * (0.5 + 0.5 * anchor_quality) /
                                    (1.0 + double(ambiguous_core_count));
                            if (actionability > best_actionability ||
                                (actionability == best_actionability &&
                                 (primary_target == -1 || target < primary_target))) {
                                primary_target = target;
                                best_actionability = actionability;
                            }
                        }
                        if (primary_target == -1) primary_target = targets.front();

                        int cycle_cost = 0;
                        for (int row = 0; row < this->check_count; row++) {
                            int overlap = 0;
                            for (auto &entry: this->pcm.iterate_row(row)) {
                                overlap += std::binary_search(
                                        support.begin(), support.end(), entry.col_index);
                            }
                            cycle_cost += overlap * (overlap - 1) / 2;
                        }
                        double support_pairs = std::max(
                                1.0, double(support.size() * (support.size() - 1)) / 2.0);
                        double cycle_density = double(cycle_cost) / support_pairs;
                        double coverage_signal = 0.0;
                        for (int target: targets) {
                            coverage_signal += this->instability_scores[target];
                        }
                        coverage_signal /= double(support.size());
                        double bounded_intervention =
                                best_actionability / (1.0 + best_actionability);
                        double utility = coverage_signal * (1.0 + bounded_intervention) /
                                (1.0 + cycle_density);
#if LDPC_MODIFIED_REGION_BP_SELECTION_MODE == 2
                        uint64_t hash = uint64_t(row_a + 1) * 11400714819323198485ull ^
                                uint64_t(row_b + 1) * 14029467366897019727ull ^
                                uint64_t(LDPC_MODIFIED_REGION_BP_RANDOM_SEED + 1) *
                                1609587929392839161ull;
                        utility = double(hash % 1000003ull + 1ull) / 1000004.0;
                        primary_target = targets[
                                (hash ^ 9650029242287828579ull) % targets.size()];
#elif LDPC_MODIFIED_REGION_BP_SELECTION_MODE == 3
                        utility = coverage_signal / (1.0 + cycle_density);
                        primary_target = *std::max_element(
                                targets.begin(), targets.end(),
                                [this](int left, int right) {
                                    if (this->instability_scores[left] !=
                                        this->instability_scores[right]) {
                                        return this->instability_scores[left] <
                                                this->instability_scores[right];
                                    }
                                    return left > right;
                                });
#endif
                        this->adaptive_candidates.push_back({
                                row_a, row_b, std::move(support), primary_target,
                                best_actionability, utility, cycle_density});
                    }
                }

                std::vector<uint8_t> feasible_core(this->bit_count, 0);
                for (const auto &candidate: this->adaptive_candidates) {
                    if (candidate.utility > 0.0) feasible_core[candidate.primary_target] = 1;
                }
                for (int target: this->region_core_variables) {
                    this->adaptive_feasible_unstable_variable_count += feasible_core[target];
                }

                std::vector<uint8_t> candidate_used(this->adaptive_candidates.size(), 0);
                std::vector<uint8_t> globally_covered(this->bit_count, 0);
                for (int view = 0; view < adaptive_max_views; view++) {
                    std::vector<uint8_t> target_selected(this->bit_count, 0);
                    std::vector<int> selected_indices;
                    int weight_limit = 0;
                    while (int(selected_indices.size()) < int(this->region_core_variables.size())) {
                        int best_index = -1;
                        double best_marginal_utility = 0.0;
                        for (int candidate_index = 0;
                             candidate_index < int(this->adaptive_candidates.size()); candidate_index++) {
                            if (candidate_used[candidate_index]) continue;
                            const auto &candidate = this->adaptive_candidates[candidate_index];
                            double view_utility = candidate.utility;
                            if (target_selected[candidate.primary_target] || view_utility <= 0.0) continue;
                            double max_jaccard = 0.0;
                            for (int selected_index: selected_indices) {
                                const auto &selected = this->adaptive_candidates[selected_index];
                                int overlap = sorted_support_overlap(
                                        candidate.support, selected.support);
                                int union_size = int(candidate.support.size() +
                                        selected.support.size()) - overlap;
                                max_jaccard = std::max(
                                        max_jaccard, double(overlap) / double(union_size));
                            }
                            double marginal_utility = view_utility / (1.0 + max_jaccard);
                            if (marginal_utility > best_marginal_utility ||
                                (marginal_utility == best_marginal_utility &&
                                 (best_index == -1 || candidate.support <
                                  this->adaptive_candidates[best_index].support))) {
                                best_index = candidate_index;
                                best_marginal_utility = marginal_utility;
                            }
                        }
                        if (best_index == -1) break;
                        candidate_used[best_index] = 1;
                        selected_indices.push_back(best_index);
                        const auto &selected = this->adaptive_candidates[best_index];
                        target_selected[selected.primary_target] = 1;
                        globally_covered[selected.primary_target] = 1;
                        weight_limit = std::max(weight_limit, int(selected.support.size()));
                    }
                    if (selected_indices.empty()) break;
                    this->adaptive_views.push_back({std::move(selected_indices), weight_limit});
                }
                for (int target: this->region_core_variables) {
                    this->adaptive_covered_unstable_variable_count += globally_covered[target];
                }
            }

            int region_core_variable_count() const {
                return int(this->region_core_variables.size());
            }

            int region_conflict_check_count() const {
                return int(this->region_conflict_checks.size());
            }

            int region_anchor_variable_count() const {
                return int(this->region_anchor_variables.size());
            }

            double region_selected_row_actionability(int view_index, int row_index) const {
                int candidate_index = this->adaptive_views.at(view_index).candidate_indices.at(row_index);
                return this->adaptive_candidates[candidate_index].actionability;
            }

            double region_selected_row_cycle_density(int view_index, int row_index) const {
                int candidate_index = this->adaptive_views.at(view_index).candidate_indices.at(row_index);
                return this->adaptive_candidates[candidate_index].cycle_density;
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
                                e.check_to_bit_msg = std::clamp(
                                        e.check_to_bit_msg, -1.0 + 1e-15, 1.0 - 1e-15);
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
                            if (this->trajectory_sample_count > 0 &&
                                this->decoding[i] == 0) {
                                this->sign_flip_count[i] += 1;
                            }
                            this->decoding[i] = 1;
                            for (auto &e: this->pcm.iterate_column(i)) {
                                this->candidate_syndrome[e.row_index] ^= 1;
                            }
                        } else {
                            if (this->trajectory_sample_count > 0 &&
                                this->decoding[i] == 1) {
                                this->sign_flip_count[i] += 1;
                            }
                            this->decoding[i] = 0;
                        }
                    }

                    this->observe_check_residuals(syndrome);
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
                                e.check_to_bit_msg = std::clamp(
                                        e.check_to_bit_msg, -1.0 + 1e-15, 1.0 - 1e-15);
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
                            if (this->trajectory_sample_count > 0 &&
                                this->decoding[bit_index] == 0) {
                                this->sign_flip_count[bit_index] += 1;
                            }
                            this->decoding[bit_index] = 1;
                        } else {
                            if (this->trajectory_sample_count > 0 &&
                                this->decoding[bit_index] == 1) {
                                this->sign_flip_count[bit_index] += 1;
                            }
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
                    this->observe_check_residuals(syndrome);
                    this->iterations = it;
                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        this->converge = true;
                        return this->decoding;
                    }
                }
                return this->decoding;
            }

        };

        class RegionAdaptiveViewDecoderPool {
        public:
            RegionAdaptiveViewDecoderPool() = default;
            RegionAdaptiveViewDecoderPool(const RegionAdaptiveViewDecoderPool &) = delete;
            RegionAdaptiveViewDecoderPool &operator=(const RegionAdaptiveViewDecoderPool &) = delete;

            void run(
                    RegionAdaptiveBpDecoder &main_decoder,
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
                            main_decoder.random_schedule_seed + view_index + 1,
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

            RegionAdaptiveBpDecoder &decoder(int view_index) {
                return *view_decoders_.at(view_index);
            }

            BpSparse &pcm(int view_index) {
                return *view_pcms_.at(view_index);
            }

        private:
            std::vector<std::unique_ptr<BpSparse>> view_pcms_;
            std::vector<std::unique_ptr<RegionAdaptiveBpDecoder>> view_decoders_;
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

            void ensure_slots(RegionAdaptiveBpDecoder &main_decoder, int count) {
                while (int(view_decoders_.size()) < count) {
                    auto pcm = std::make_unique<BpSparse>(
                            main_decoder.check_count, main_decoder.bit_count, 0);
                    BpSparse *pcm_pointer = pcm.get();
                    view_pcms_.push_back(std::move(pcm));
                    view_decoders_.push_back(std::make_unique<RegionAdaptiveBpDecoder>(
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
}  // namespace ldpc_modified::region_adaptive_bp

#endif
