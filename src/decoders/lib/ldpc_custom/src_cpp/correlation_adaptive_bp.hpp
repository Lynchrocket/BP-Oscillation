#ifndef LDPC_MODIFIED_CORRELATION_ADAPTIVE_BP_HPP
#define LDPC_MODIFIED_CORRELATION_ADAPTIVE_BP_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace ldpc_modified {
namespace correlation_adaptive_bp {

class CorrelationAdaptiveRegionRefiner {
  public:
    std::vector<uint8_t> decoding;
    std::vector<int> selected_region_variables;
    std::vector<double> logical_log_masses;
    int structural_atom_count = 0;
    int eligible_atom_count = 0;
    int probe_region_count = 0;
    int candidate_evaluation_count = 0;
    int selected_atom_count = 0;
    uint64_t probe_enumerated_state_count = 0;
    uint64_t selected_region_state_count = 0;
    uint64_t selected_feasible_state_count = 0;
    double total_correlation_bits = 0.0;
    double marginal_kl_bits = 0.0;
    double factorization_gap_bits = 0.0;
    double selected_utility = 0.0;
    bool main_syndrome_valid = false;
    bool used_feasible_seed = false;
    bool refined = false;
    bool syndrome_valid = false;
    int selected_logical_label = 0;

    CorrelationAdaptiveRegionRefiner(
        const std::vector<std::vector<uint8_t>>& pcm,
        const std::vector<std::vector<uint8_t>>& logicals,
        const std::vector<double>& priors,
        const std::vector<int>& check_times,
        const std::vector<int>& check_groups,
        const std::vector<int>& check_layers,
        uint64_t max_region_states,
        int max_probe_regions,
        int region_mode,
        int score_mode,
        int region_seed,
        bool refine_zero_syndrome)
        : h_(pcm), l_(logicals), priors_(priors), times_(check_times),
          groups_(check_groups), layers_(check_layers),
          max_region_states_(max_region_states), max_probe_regions_(max_probe_regions),
          region_mode_(region_mode), score_mode_(score_mode), region_seed_(region_seed),
          refine_zero_syndrome_(refine_zero_syndrome) {
        m_ = static_cast<int>(h_.size());
        n_ = m_ == 0 ? static_cast<int>(priors_.size()) : static_cast<int>(h_[0].size());
        q_ = static_cast<int>(l_.size());
        if (n_ <= 0 || static_cast<int>(priors_.size()) != n_ ||
            max_region_states_ < 2 || max_region_states_ > (uint64_t{1} << 62) ||
            max_probe_regions_ <= 0) {
            throw std::invalid_argument("invalid correlation-adaptive dimensions or budgets");
        }
        if (region_mode_ < 0 || region_mode_ > 2 || score_mode_ < 0 || score_mode_ > 4) {
            throw std::invalid_argument("invalid correlation-adaptive mode");
        }
        if (q_ > 20) throw std::invalid_argument("at most 20 observables are supported");
        for (const auto& row : h_) {
            if (static_cast<int>(row.size()) != n_) throw std::invalid_argument("ragged pcm");
        }
        for (const auto& row : l_) {
            if (static_cast<int>(row.size()) != n_) throw std::invalid_argument("ragged logical matrix");
        }
        if (static_cast<int>(times_.size()) != m_ || static_cast<int>(groups_.size()) != m_ ||
            static_cast<int>(layers_.size()) != m_) {
            throw std::invalid_argument("invalid check metadata");
        }
        for (double p : priors_) {
            if (!std::isfinite(p) || p <= 0.0 || p >= 1.0) {
                throw std::invalid_argument(
                    "priors must be finite probabilities strictly between zero and one");
            }
        }
        build_graph();
        build_elimination_and_atoms();
    }

    std::vector<uint8_t>& refine(const std::vector<uint8_t>& syndrome,
                                 const std::vector<uint8_t>& main_decoding,
                                 const std::vector<double>& posterior_llr) {
        if (static_cast<int>(syndrome.size()) != m_ ||
            static_cast<int>(main_decoding.size()) != n_ ||
            static_cast<int>(posterior_llr.size()) != n_) {
            throw std::invalid_argument("invalid correlation-adaptive decode input dimensions");
        }
        for (double llr : posterior_llr) {
            if (std::isnan(llr)) {
                throw std::invalid_argument("posterior LLR contains NaN");
            }
        }
        reset_stats();
        main_syndrome_valid = valid_syndrome(main_decoding, syndrome);
        const std::vector<uint8_t> seed = make_valid_seed(syndrome, main_decoding);
        used_feasible_seed = seed != main_decoding;
        decoding = seed;

        const bool zero_syndrome = std::none_of(
            syndrome.begin(), syndrome.end(), [](uint8_t bit) { return bit != 0; });
        if (zero_syndrome && !refine_zero_syndrome_) return finish(main_decoding, syndrome);
        const uint64_t seed_label = logical_label(seed);

        const std::vector<int> shortlist = candidate_shortlist(
            syndrome, main_decoding, seed, posterior_llr);
        probe_region_count = static_cast<int>(shortlist.size());
        if (shortlist.empty()) return finish(main_decoding, syndrome);

        Evaluation best;
        int best_atom = -1;
        for (int atom_index : shortlist) {
            Evaluation evaluation = evaluate(
                atoms_[atom_index].variables, seed, syndrome, posterior_llr,
                atoms_[atom_index].structural_score, seed_label);
            if (!evaluation.selectable) continue;
            if (best_atom < 0 || better_evaluation(evaluation, best)) {
                best = std::move(evaluation);
                best_atom = atom_index;
            }
        }
        if (best_atom < 0 || !(best.utility > 0.0)) {
            return finish(main_decoding, syndrome);
        }

        std::vector<uint8_t> selected_atom(atoms_.size(), 0);
        selected_atom[best_atom] = 1;
        selected_atom_count = 1;
        bool grew = true;
        while (grew) {
            grew = false;
            Evaluation best_union;
            int best_union_atom = -1;
            for (int atom_index : shortlist) {
                if (selected_atom[atom_index]) continue;
                if (!regions_adjacent(best.variables, atoms_[atom_index].variables)) continue;
                std::vector<int> merged = sorted_union(
                    best.variables, atoms_[atom_index].variables);
                if (merged.size() == best.variables.size() ||
                    merged.size() >= 63 ||
                    (uint64_t{1} << merged.size()) > max_region_states_) {
                    continue;
                }
                Evaluation evaluation = evaluate(
                    merged, seed, syndrome, posterior_llr,
                    region_structural_score(merged), seed_label);
                if (!evaluation.selectable || !(evaluation.utility > best.utility)) continue;
                if (best_union_atom < 0 || better_evaluation(evaluation, best_union)) {
                    best_union = std::move(evaluation);
                    best_union_atom = atom_index;
                }
            }
            if (best_union_atom >= 0) {
                best = std::move(best_union);
                selected_atom[best_union_atom] = 1;
                ++selected_atom_count;
                grew = true;
            }
        }

        selected_region_variables = best.variables;
        selected_region_state_count = best.enumerated_states;
        selected_feasible_state_count = best.feasible_states;
        total_correlation_bits = best.total_correlation_bits;
        marginal_kl_bits = best.marginal_kl_bits;
        factorization_gap_bits = best.factorization_gap_bits;
        selected_utility = best.utility;
        selected_logical_label = static_cast<int>(best.selected_logical_label);
        logical_log_masses.assign(size_t{1} << q_, neg_inf());
        for (const auto& [label, mass] : best.logical_log_masses) {
            logical_log_masses[label] = mass;
        }
        decoding = seed;
        for (size_t i = 0; i < best.variables.size(); ++i) {
            decoding[best.variables[i]] = (best.selected_assignment >> i) & uint64_t{1};
        }
        return finish(main_decoding, syndrome);
    }

  private:
    struct Atom {
        std::vector<int> variables;
        double structural_score = 0.0;
    };

    struct Evaluation {
        std::vector<int> variables;
        uint64_t enumerated_states = 0;
        uint64_t feasible_states = 0;
        double structural_score = 0.0;
        double total_correlation_bits = 0.0;
        double marginal_kl_bits = 0.0;
        double factorization_gap_bits = 0.0;
        double selection_score = 0.0;
        double utility = 0.0;
        uint64_t selected_logical_label = 0;
        uint64_t selected_assignment = 0;
        std::map<uint64_t, double> logical_log_masses;
        bool selectable = false;
    };

    struct CandidateRank {
        int atom_index = -1;
        int residual_support = 0;
        int repaired_support = 0;
        double structural_score = 0.0;
        double ambiguity = 0.0;
        int negative_size = 0;
    };

    struct FeasibleState {
        uint64_t state = 0;
        uint64_t logical_label = 0;
        double log_weight = 0.0;
    };

    int m_ = 0;
    int n_ = 0;
    int q_ = 0;
    std::vector<std::vector<uint8_t>> h_;
    std::vector<std::vector<uint8_t>> l_;
    std::vector<double> priors_;
    std::vector<int> times_;
    std::vector<int> groups_;
    std::vector<int> layers_;
    uint64_t max_region_states_;
    int max_probe_regions_;
    int region_mode_;
    int score_mode_;
    int region_seed_;
    bool refine_zero_syndrome_;
    std::vector<std::vector<int>> check_vars_;
    std::vector<std::vector<int>> var_checks_;
    std::vector<int> cycle_count_;
    std::vector<int> structural_cycle_count_;
    std::vector<std::vector<uint8_t>> rref_;
    std::vector<std::vector<uint8_t>> row_transform_;
    std::vector<int> pivot_columns_;
    std::vector<Atom> atoms_;
    std::vector<uint64_t> variable_logical_mask_;

    static constexpr double log_two() { return 0.693147180559945309417232121458176568; }
    static double neg_inf() { return -std::numeric_limits<double>::infinity(); }

    static int trailing_zero_index(uint64_t value) {
        assert(value != 0);
        int index = 0;
        while ((value & uint64_t{1}) == 0) {
            value >>= 1;
            ++index;
        }
        return index;
    }

    static double log_add(double a, double b) {
        if (!std::isfinite(a)) return b;
        if (!std::isfinite(b)) return a;
        const double hi = std::max(a, b);
        return hi + std::log1p(std::exp(std::min(a, b) - hi));
    }

    static double binary_entropy_nats(double p) {
        assert(p >= 0.0 && p <= 1.0);
        double result = 0.0;
        if (p > 0.0) result -= p * std::log(p);
        if (p < 1.0) result -= (1.0 - p) * std::log1p(-p);
        return result;
    }

    static double log_bp_probability(double llr, uint8_t bit) {
        if (std::isnan(llr)) throw std::runtime_error("BP posterior LLR is NaN");
        if (bit) {
            if (llr == std::numeric_limits<double>::infinity()) return neg_inf();
            if (llr == -std::numeric_limits<double>::infinity()) return 0.0;
            if (llr >= 0.0) return -llr - std::log1p(std::exp(-llr));
            return -std::log1p(std::exp(llr));
        }
        if (llr == std::numeric_limits<double>::infinity()) return 0.0;
        if (llr == -std::numeric_limits<double>::infinity()) return neg_inf();
        if (llr >= 0.0) return -std::log1p(std::exp(-llr));
        return llr - std::log1p(std::exp(llr));
    }

    static double binary_kl_from_llr_bits(double p, double llr) {
        assert(p >= 0.0 && p <= 1.0);
        const double log_q0 = log_bp_probability(llr, 0);
        const double log_q1 = log_bp_probability(llr, 1);
        if ((p < 1.0 && !std::isfinite(log_q0)) || (p > 0.0 && !std::isfinite(log_q1))) {
            return std::numeric_limits<double>::infinity();
        }
        double cross_entropy = 0.0;
        if (p < 1.0) cross_entropy -= (1.0 - p) * log_q0;
        if (p > 0.0) cross_entropy -= p * log_q1;
        const double kl_nats = cross_entropy - binary_entropy_nats(p);
        if (kl_nats < -1e-10) throw std::runtime_error("negative marginal KL divergence");
        return std::max(0.0, kl_nats) / log_two();
    }

    void build_graph() {
        check_vars_.assign(m_, {});
        var_checks_.assign(n_, {});
        for (int c = 0; c < m_; ++c) {
            for (int v = 0; v < n_; ++v) {
                if (!(h_[c][v] & uint8_t{1})) continue;
                check_vars_[c].push_back(v);
                var_checks_[v].push_back(c);
            }
        }

        cycle_count_.assign(n_, 0);
        structural_cycle_count_.assign(n_, 0);
        for (int a = 0; a < m_; ++a) {
            for (int b = a + 1; b < m_; ++b) {
                std::vector<int> shared;
                std::set_intersection(
                    check_vars_[a].begin(), check_vars_[a].end(),
                    check_vars_[b].begin(), check_vars_[b].end(),
                    std::back_inserter(shared));
                if (shared.size() < 2) continue;
                const bool temporal = groups_[a] >= 0 && groups_[a] == groups_[b] &&
                    times_[a] >= 0 && times_[a] != times_[b];
                const bool spatial = times_[a] >= 0 && times_[a] == times_[b] &&
                    groups_[a] >= 0 && groups_[a] != groups_[b];
                const bool recursive = layers_[a] >= 0 && layers_[b] >= 0 &&
                    layers_[a] != layers_[b];
                for (int v : shared) {
                    ++cycle_count_[v];
                    if (temporal || spatial || recursive) ++structural_cycle_count_[v];
                }
            }
        }

        variable_logical_mask_.assign(n_, 0);
        for (int logical = 0; logical < q_; ++logical) {
            for (int v = 0; v < n_; ++v) {
                if (l_[logical][v]) variable_logical_mask_[v] |= uint64_t{1} << logical;
            }
        }
    }

    double region_structural_score(const std::vector<int>& variables) const {
        double score = 0.0;
        for (int variable : variables) {
            score += 1.0 + cycle_count_[variable] + structural_cycle_count_[variable];
        }
        return score;
    }

    void add_atom(std::set<std::vector<int>>& unique, std::vector<int> support) {
        std::sort(support.begin(), support.end());
        support.erase(std::unique(support.begin(), support.end()), support.end());
        if (support.empty() || !unique.insert(support).second) return;
        ++structural_atom_count;
        for (int c = 0; c < m_; ++c) {
            uint8_t parity = 0;
            for (int v : support) parity ^= h_[c][v];
            if (parity) throw std::runtime_error("structural atom is not in the PCM kernel");
        }
        if (support.size() >= 63 ||
            (uint64_t{1} << support.size()) > max_region_states_) return;
        const double structural_score = region_structural_score(support);
        atoms_.push_back({std::move(support), structural_score});
        ++eligible_atom_count;
    }

    void build_elimination_and_atoms() {
        rref_ = h_;
        row_transform_.assign(m_, std::vector<uint8_t>(m_, 0));
        for (int row = 0; row < m_; ++row) row_transform_[row][row] = 1;
        pivot_columns_.clear();
        int row = 0;
        for (int col = 0; col < n_ && row < m_; ++col) {
            int pivot = row;
            while (pivot < m_ && !rref_[pivot][col]) ++pivot;
            if (pivot == m_) continue;
            std::swap(rref_[row], rref_[pivot]);
            std::swap(row_transform_[row], row_transform_[pivot]);
            for (int other = 0; other < m_; ++other) {
                if (other == row || !rref_[other][col]) continue;
                for (int j = 0; j < n_; ++j) rref_[other][j] ^= rref_[row][j];
                for (int j = 0; j < m_; ++j) {
                    row_transform_[other][j] ^= row_transform_[row][j];
                }
            }
            pivot_columns_.push_back(col);
            ++row;
        }

        std::set<std::vector<int>> unique;
        std::map<std::vector<int>, std::vector<int>> equal_columns;
        for (int v = 0; v < n_; ++v) equal_columns[var_checks_[v]].push_back(v);
        for (const auto& [boundary, variables] : equal_columns) {
            static_cast<void>(boundary);
            for (size_t i = 1; i < variables.size(); ++i) {
                add_atom(unique, {variables.front(), variables[i]});
            }
        }

        std::vector<uint8_t> is_pivot(n_, 0);
        for (int pivot : pivot_columns_) is_pivot[pivot] = 1;
        for (int free_column = 0; free_column < n_; ++free_column) {
            if (is_pivot[free_column]) continue;
            std::vector<int> support{free_column};
            for (size_t pivot_row = 0; pivot_row < pivot_columns_.size(); ++pivot_row) {
                if (rref_[pivot_row][free_column]) {
                    support.push_back(pivot_columns_[pivot_row]);
                }
            }
            add_atom(unique, std::move(support));
        }
    }

    std::vector<uint8_t> make_valid_seed(const std::vector<uint8_t>& syndrome,
                                         const std::vector<uint8_t>& main) const {
        if (valid_syndrome(main, syndrome)) return main;
        std::vector<uint8_t> residual(m_, 0);
        for (int c = 0; c < m_; ++c) {
            residual[c] = syndrome[c];
            for (int v : check_vars_[c]) residual[c] ^= main[v];
        }
        std::vector<uint8_t> transformed(m_, 0);
        for (int row = 0; row < m_; ++row) {
            for (int source = 0; source < m_; ++source) {
                transformed[row] ^= row_transform_[row][source] & residual[source];
            }
        }
        for (int row = static_cast<int>(pivot_columns_.size()); row < m_; ++row) {
            if (transformed[row]) throw std::runtime_error("inconsistent syndrome");
        }
        std::vector<uint8_t> seed = main;
        for (size_t row = 0; row < pivot_columns_.size(); ++row) {
            seed[pivot_columns_[row]] ^= transformed[row];
        }
        if (!valid_syndrome(seed, syndrome)) {
            throw std::runtime_error("failed to construct a syndrome-valid seed");
        }
        return seed;
    }

    std::vector<int> candidate_shortlist(const std::vector<uint8_t>& syndrome,
                                         const std::vector<uint8_t>& main,
                                         const std::vector<uint8_t>& seed,
                                         const std::vector<double>& llr) const {
        std::vector<uint8_t> residual_check(m_, 0);
        for (int c = 0; c < m_; ++c) {
            residual_check[c] = syndrome[c];
            for (int v : check_vars_[c]) residual_check[c] ^= main[v];
        }
        std::vector<int> order;
        if (region_mode_ == 2) {
            order.resize(atoms_.size());
            std::iota(order.begin(), order.end(), 0);
            std::mt19937 rng(static_cast<uint32_t>(region_seed_));
            std::shuffle(order.begin(), order.end(), rng);
        } else {
            std::vector<int> variable_residual_support(n_, 0);
            std::vector<uint8_t> variable_repaired_support(n_, 0);
            std::vector<double> variable_ambiguity(n_, 0.0);
            for (int v = 0; v < n_; ++v) {
                variable_repaired_support[v] = main[v] != seed[v];
                variable_ambiguity[v] = 1.0 / (1.0 + std::abs(llr[v]));
                for (int c : var_checks_[v]) {
                    variable_residual_support[v] += residual_check[c];
                }
            }

            std::vector<CandidateRank> ranked;
            ranked.reserve(atoms_.size());
            for (int index = 0; index < static_cast<int>(atoms_.size()); ++index) {
                CandidateRank rank;
                rank.atom_index = index;
                rank.structural_score = region_mode_ == 1
                    ? atoms_[index].structural_score : 0.0;
                rank.negative_size = -static_cast<int>(atoms_[index].variables.size());
                for (int v : atoms_[index].variables) {
                    rank.residual_support += variable_residual_support[v];
                    rank.repaired_support += variable_repaired_support[v];
                    rank.ambiguity += variable_ambiguity[v];
                }
                ranked.push_back(rank);
            }

            const auto better_rank = [&](const CandidateRank& left, const CandidateRank& right) {
                if (left.residual_support != right.residual_support) {
                    return left.residual_support > right.residual_support;
                }
                if (left.repaired_support != right.repaired_support) {
                    return left.repaired_support > right.repaired_support;
                }
                if (left.structural_score != right.structural_score) {
                    return left.structural_score > right.structural_score;
                }
                if (left.ambiguity != right.ambiguity) {
                    return left.ambiguity > right.ambiguity;
                }
                if (left.negative_size != right.negative_size) {
                    return left.negative_size > right.negative_size;
                }
                return atoms_[left.atom_index].variables > atoms_[right.atom_index].variables;
            };

            const size_t keep = std::min(
                ranked.size(), static_cast<size_t>(max_probe_regions_));
            if (keep < ranked.size()) {
                std::partial_sort(ranked.begin(), ranked.begin() + keep, ranked.end(), better_rank);
                ranked.resize(keep);
            } else {
                std::sort(ranked.begin(), ranked.end(), better_rank);
            }
            order.clear();
            order.reserve(ranked.size());
            for (const CandidateRank& rank : ranked) order.push_back(rank.atom_index);
        }
        if (order.size() > static_cast<size_t>(max_probe_regions_)) {
            order.resize(max_probe_regions_);
        }
        return order;
    }

    Evaluation evaluate(const std::vector<int>& variables,
                        const std::vector<uint8_t>& seed,
                        const std::vector<uint8_t>& syndrome,
                        const std::vector<double>& llr,
                        double structural_score,
                        uint64_t seed_label) {
        static_cast<void>(syndrome);
        Evaluation result;
        result.variables = variables;
        result.structural_score = structural_score;
        const size_t k = variables.size();
        if (k == 0 || k >= 63) return result;
        result.enumerated_states = uint64_t{1} << k;
        if (result.enumerated_states > max_region_states_) return result;
        ++candidate_evaluation_count;
        probe_enumerated_state_count += result.enumerated_states;

        std::vector<int> region_checks;
        std::vector<uint8_t> seen_check(m_, 0);
        for (int v : variables) {
            for (int c : var_checks_[v]) {
                if (!seen_check[c]) {
                    seen_check[c] = 1;
                    region_checks.push_back(c);
                }
            }
        }
        const size_t constraint_word_count = (region_checks.size() + 63) / 64;
        std::vector<uint64_t> check_variable_masks(region_checks.size(), 0);
        std::vector<uint64_t> required_words(constraint_word_count, 0);
        std::vector<std::vector<uint64_t>> variable_check_words(
            k, std::vector<uint64_t>(constraint_word_count, 0));
        for (size_t check_index = 0; check_index < region_checks.size(); ++check_index) {
            const int c = region_checks[check_index];
            uint8_t required = 0;
            uint64_t variable_mask = 0;
            for (size_t i = 0; i < k; ++i) {
                if (!h_[c][variables[i]]) continue;
                variable_mask |= uint64_t{1} << i;
                required ^= seed[variables[i]];
            }
            assert(variable_mask != 0);
            check_variable_masks[check_index] = variable_mask;
            const size_t word = check_index / 64;
            const uint64_t check_bit = uint64_t{1} << (check_index % 64);
            if (required) required_words[word] |= check_bit;
        }
        for (size_t check_index = 0; check_index < region_checks.size(); ++check_index) {
            const size_t word = check_index / 64;
            const uint64_t check_bit = uint64_t{1} << (check_index % 64);
            uint64_t remaining = check_variable_masks[check_index];
            while (remaining) {
                const int local_variable = trailing_zero_index(remaining);
                variable_check_words[local_variable][word] |= check_bit;
                remaining &= remaining - 1;
            }
        }

        std::vector<uint64_t> mismatch_words = required_words;
        size_t nonzero_mismatch_word_count = 0;
        for (uint64_t word : mismatch_words) nonzero_mismatch_word_count += word != 0;

        std::vector<double> log_zero(k, 0.0);
        std::vector<double> log_one(k, 0.0);
        uint64_t current_logical_label = seed_label;
        for (size_t i = 0; i < k; ++i) {
            const int variable = variables[i];
            log_zero[i] = std::log1p(-priors_[variable]);
            log_one[i] = std::log(priors_[variable]);
            if (seed[variable]) current_logical_label ^= variable_logical_mask_[variable];
        }

        std::vector<FeasibleState> feasible_records;
        uint64_t gray_state = 0;
        for (uint64_t step = 0; step < result.enumerated_states; ++step) {
            if (step != 0) {
                const int changed = trailing_zero_index(step);
                gray_state = step ^ (step >> 1);
                current_logical_label ^= variable_logical_mask_[variables[changed]];
                for (size_t word = 0; word < constraint_word_count; ++word) {
                    const uint64_t variable_word = variable_check_words[changed][word];
                    if (!variable_word) continue;
                    const bool was_nonzero = mismatch_words[word] != 0;
                    mismatch_words[word] ^= variable_word;
                    const bool is_nonzero = mismatch_words[word] != 0;
                    if (was_nonzero && !is_nonzero) --nonzero_mismatch_word_count;
                    if (!was_nonzero && is_nonzero) ++nonzero_mismatch_word_count;
                }
            }
            if (nonzero_mismatch_word_count != 0) continue;
            double log_weight = 0.0;
            for (size_t i = 0; i < k; ++i) {
                log_weight += ((gray_state >> i) & uint64_t{1}) ? log_one[i] : log_zero[i];
            }
            if (!std::isfinite(log_weight)) continue;
            feasible_records.push_back({gray_state, current_logical_label, log_weight});
        }
        std::sort(
            feasible_records.begin(), feasible_records.end(),
            [](const FeasibleState& left, const FeasibleState& right) {
                return left.state < right.state;
            });

        std::vector<uint64_t> assignments;
        std::vector<double> log_weights;
        std::map<uint64_t, double> best_log_weight;
        std::map<uint64_t, uint64_t> best_assignment;
        assignments.reserve(feasible_records.size());
        log_weights.reserve(feasible_records.size());
        for (const FeasibleState& record : feasible_records) {
            assignments.push_back(record.state);
            log_weights.push_back(record.log_weight);
            auto mass = result.logical_log_masses.find(record.logical_label);
            if (mass == result.logical_log_masses.end()) {
                result.logical_log_masses[record.logical_label] = record.log_weight;
            } else {
                mass->second = log_add(mass->second, record.log_weight);
            }
            auto best_weight = best_log_weight.find(record.logical_label);
            if (best_weight == best_log_weight.end() ||
                record.log_weight > best_weight->second ||
                (record.log_weight == best_weight->second &&
                 record.state < best_assignment[record.logical_label])) {
                best_log_weight[record.logical_label] = record.log_weight;
                best_assignment[record.logical_label] = record.state;
            }
        }
        result.feasible_states = assignments.size();
        if (assignments.empty()) {
            throw std::runtime_error("selected region has no positive-probability feasible assignment");
        }

        double log_partition = neg_inf();
        for (double weight : log_weights) log_partition = log_add(log_partition, weight);
        std::vector<double> marginal_one(k, 0.0);
        double joint_entropy_nats = 0.0;
        for (size_t state_index = 0; state_index < assignments.size(); ++state_index) {
            const double probability = std::exp(log_weights[state_index] - log_partition);
            joint_entropy_nats -= probability * (log_weights[state_index] - log_partition);
            for (size_t i = 0; i < k; ++i) {
                if ((assignments[state_index] >> i) & uint64_t{1}) marginal_one[i] += probability;
            }
        }
        double marginal_entropy_nats = 0.0;
        double marginal_kl = 0.0;
        for (size_t i = 0; i < k; ++i) {
            marginal_one[i] = std::min(1.0, std::max(0.0, marginal_one[i]));
            marginal_entropy_nats += binary_entropy_nats(marginal_one[i]);
            marginal_kl += binary_kl_from_llr_bits(marginal_one[i], llr[variables[i]]);
        }
        const double correlation_nats = marginal_entropy_nats - joint_entropy_nats;
        if (correlation_nats < -1e-10) {
            throw std::runtime_error("negative total correlation");
        }
        result.total_correlation_bits = std::max(0.0, correlation_nats) / log_two();
        result.marginal_kl_bits = marginal_kl;
        result.factorization_gap_bits = result.total_correlation_bits + result.marginal_kl_bits;
        if (score_mode_ == 0) result.selection_score = result.factorization_gap_bits;
        if (score_mode_ == 1) result.selection_score = result.total_correlation_bits;
        if (score_mode_ == 2) result.selection_score = result.marginal_kl_bits;
        if (score_mode_ == 3) result.selection_score = result.structural_score;
        if (score_mode_ == 4) {
            uint64_t hash = static_cast<uint64_t>(static_cast<uint32_t>(region_seed_)) +
                0x9e3779b97f4a7c15ULL;
            for (int variable : variables) {
                hash ^= static_cast<uint64_t>(variable + 1) + 0x9e3779b97f4a7c15ULL +
                    (hash << 6) + (hash >> 2);
            }
            result.selection_score = static_cast<double>(hash % 1000003ULL + 1ULL) /
                1000004.0;
        }
        result.utility = result.selection_score / static_cast<double>(result.enumerated_states);
        result.selectable = result.feasible_states >= 2;

        auto best_label = result.logical_log_masses.begin();
        for (auto it = std::next(result.logical_log_masses.begin());
             it != result.logical_log_masses.end(); ++it) {
            if (it->second > best_label->second) best_label = it;
        }
        result.selected_logical_label = best_label->first;
        result.selected_assignment = best_assignment[result.selected_logical_label];
        return result;
    }

    static bool better_evaluation(const Evaluation& left, const Evaluation& right) {
        if (left.utility != right.utility) return left.utility > right.utility;
        if (left.selection_score != right.selection_score) {
            return left.selection_score > right.selection_score;
        }
        if (left.enumerated_states != right.enumerated_states) {
            return left.enumerated_states < right.enumerated_states;
        }
        return left.variables < right.variables;
    }

    bool regions_adjacent(const std::vector<int>& left,
                          const std::vector<int>& right) const {
        size_t i = 0;
        size_t j = 0;
        while (i < left.size() && j < right.size()) {
            if (left[i] == right[j]) return true;
            if (left[i] < right[j]) ++i;
            else ++j;
        }
        std::vector<uint8_t> left_check(m_, 0);
        for (int v : left) for (int c : var_checks_[v]) left_check[c] = 1;
        for (int v : right) {
            for (int c : var_checks_[v]) if (left_check[c]) return true;
        }
        return false;
    }

    static std::vector<int> sorted_union(const std::vector<int>& left,
                                         const std::vector<int>& right) {
        std::vector<int> result;
        result.reserve(left.size() + right.size());
        std::set_union(left.begin(), left.end(), right.begin(), right.end(),
                       std::back_inserter(result));
        return result;
    }

    uint64_t logical_label(const std::vector<uint8_t>& error) const {
        uint64_t label = 0;
        for (int v = 0; v < n_; ++v) if (error[v]) label ^= variable_logical_mask_[v];
        return label;
    }

    bool valid_syndrome(const std::vector<uint8_t>& error,
                        const std::vector<uint8_t>& syndrome) const {
        for (int c = 0; c < m_; ++c) {
            uint8_t parity = 0;
            for (int v : check_vars_[c]) parity ^= error[v];
            if (parity != syndrome[c]) return false;
        }
        return true;
    }

    std::vector<uint8_t>& finish(const std::vector<uint8_t>& main,
                                 const std::vector<uint8_t>& syndrome) {
        refined = decoding != main;
        syndrome_valid = valid_syndrome(decoding, syndrome);
        if (!syndrome_valid) {
            throw std::runtime_error("correlation-adaptive refiner returned an invalid syndrome");
        }
        return decoding;
    }

    void reset_stats() {
        decoding.clear();
        selected_region_variables.clear();
        logical_log_masses.clear();
        probe_region_count = candidate_evaluation_count = selected_atom_count = 0;
        probe_enumerated_state_count = selected_region_state_count = 0;
        selected_feasible_state_count = 0;
        total_correlation_bits = marginal_kl_bits = factorization_gap_bits = 0.0;
        selected_utility = 0.0;
        main_syndrome_valid = used_feasible_seed = refined = syndrome_valid = false;
        selected_logical_label = 0;
    }
};

}  // namespace correlation_adaptive_bp
}  // namespace ldpc_modified

#endif
