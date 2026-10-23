#ifndef LDPC_MODIFIED_RECURSIVE_FRONTIER_BP_HPP
#define LDPC_MODIFIED_RECURSIVE_FRONTIER_BP_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ldpc_modified {
namespace recursive_frontier_bp {

class RecursiveFrontierRefiner {
  public:
    std::vector<uint8_t> decoding;
    std::vector<double> logical_log_masses;
    uint64_t expanded_transition_count = 0;
    uint64_t merged_path_count = 0;
    uint64_t pruned_state_count = 0;
    uint64_t closed_check_rejection_count = 0;
    uint64_t terminal_state_count = 0;
    uint64_t peak_state_count = 0;
    int peak_frontier_check_count = 0;
    int recursive_bag_count = 0;
    int peak_recursive_bag_size = 0;
    int selected_logical_label = 0;
    bool truncated = false;
    bool used_main_fallback = false;
    bool refined = false;
    bool syndrome_valid = false;

    RecursiveFrontierRefiner(
        const std::vector<std::vector<uint8_t>>& pcm,
        const std::vector<std::vector<uint8_t>>& logicals,
        const std::vector<double>& priors,
        const std::vector<int>& variable_order,
        const std::vector<int>& bag_ends,
        uint64_t max_frontier_states)
        : h_(pcm), l_(logicals), priors_(priors), order_(variable_order),
          bag_ends_(bag_ends),
          max_frontier_states_(max_frontier_states) {
        if (h_.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
            priors_.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
            l_.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
            (!h_.empty() && h_[0].size() > static_cast<size_t>(std::numeric_limits<int>::max()))) {
            throw std::invalid_argument("frontier dimensions exceed integer indexing");
        }
        m_ = static_cast<int>(h_.size());
        n_ = m_ == 0 ? static_cast<int>(priors_.size()) : static_cast<int>(h_[0].size());
        q_ = static_cast<int>(l_.size());
        if (static_cast<int>(priors_.size()) != n_ || max_frontier_states_ == 0 ||
            max_frontier_states_ > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("invalid frontier dimensions or state budget");
        }
        if (q_ > 20) throw std::invalid_argument("at most 20 observables are supported");
        for (const auto& row : h_) {
            if (static_cast<int>(row.size()) != n_) throw std::invalid_argument("ragged pcm");
        }
        for (const auto& row : l_) {
            if (static_cast<int>(row.size()) != n_) throw std::invalid_argument("ragged logical matrix");
        }
        for (double p : priors_) {
            if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
                throw std::invalid_argument("priors must be finite probabilities");
            }
        }
        validate_order();
        validate_bags();
        build_schedule();
    }

    std::vector<uint8_t>& refine(const std::vector<uint8_t>& syndrome,
                                 const std::vector<uint8_t>& main_decoding,
                                 const std::vector<double>& posterior_llr) {
        if (static_cast<int>(syndrome.size()) != m_ ||
            static_cast<int>(main_decoding.size()) != n_ ||
            static_cast<int>(posterior_llr.size()) != n_) {
            throw std::invalid_argument("invalid frontier decode input dimensions");
        }
        reset_stats();
        const bool main_valid = valid_syndrome(main_decoding, syndrome);
        for (int c : empty_checks_) {
            if (syndrome[c]) {
                decoding = main_decoding;
                used_main_fallback = true;
                syndrome_valid = false;
                return decoding;
            }
        }

        State initial;
        initial.key.parity.clear();
        initial.log_mass = 0.0;
        initial.best_log_weight = 0.0;
        initial.bp_tie_score = 0.0;
        std::vector<State> states{initial};
        Key anchor_key = initial.key;
        std::vector<std::vector<int>> history_prev;
        std::vector<std::vector<uint8_t>> history_bit;
        history_prev.reserve(n_);
        history_bit.reserve(n_);

        for (int position = 0; position < n_; ++position) {
            const int variable = order_[position];
            std::map<Key, int> merged_indices;
            std::vector<State> candidates;
            candidates.reserve(std::min<uint64_t>(2 * states.size(), max_frontier_states_ * 2));
            for (int state_index = 0; state_index < static_cast<int>(states.size()); ++state_index) {
                for (uint8_t bit = 0; bit <= 1; ++bit) {
                    ++expanded_transition_count;
                    const double branch_weight = bit ? log_probability_one_[variable]
                                                     : log_probability_zero_[variable];
                    if (!std::isfinite(branch_weight)) continue;
                    State candidate = states[state_index];
                    candidate.prev_index = state_index;
                    candidate.bit = bit;
                    candidate.log_mass += branch_weight;
                    candidate.best_log_weight += branch_weight;
                    candidate.bp_tie_score += bp_preference(posterior_llr[variable], bit);
                    Key transitioned;
                    if (!transition_key(candidate.key, transitioned, position, variable, bit, syndrome)) {
                        ++closed_check_rejection_count;
                        continue;
                    }
                    candidate.key = std::move(transitioned);
                    auto [entry, inserted] = merged_indices.emplace(candidate.key, candidates.size());
                    if (inserted) {
                        candidates.push_back(std::move(candidate));
                    } else {
                        ++merged_path_count;
                        State& existing = candidates[entry->second];
                        existing.log_mass = log_add(existing.log_mass, candidate.log_mass);
                        if (better_map_path(candidate, existing)) {
                            existing.best_log_weight = candidate.best_log_weight;
                            existing.bp_tie_score = candidate.bp_tie_score;
                            existing.prev_index = candidate.prev_index;
                            existing.bit = candidate.bit;
                        }
                    }
                }
            }

            if (main_valid) advance_anchor(
                anchor_key, position, variable, main_decoding[variable], syndrome
            );
            std::vector<int> keep(candidates.size());
            std::iota(keep.begin(), keep.end(), 0);
            std::stable_sort(keep.begin(), keep.end(), [&](int left, int right) {
                return better_beam_state(
                    candidates[left], candidates[right], bag_end_at_[position]
                );
            });
            if (keep.size() > max_frontier_states_) {
                truncated = true;
                pruned_state_count += keep.size() - max_frontier_states_;
                keep.resize(static_cast<size_t>(max_frontier_states_));
                if (main_valid) preserve_anchor(keep, candidates, anchor_key);
            }

            std::vector<State> next_states;
            std::vector<int> stage_prev;
            std::vector<uint8_t> stage_bit;
            next_states.reserve(keep.size());
            stage_prev.reserve(keep.size());
            stage_bit.reserve(keep.size());
            for (int index : keep) {
                stage_prev.push_back(candidates[index].prev_index);
                stage_bit.push_back(candidates[index].bit);
                next_states.push_back(std::move(candidates[index]));
            }
            states = std::move(next_states);
            history_prev.push_back(std::move(stage_prev));
            history_bit.push_back(std::move(stage_bit));
            peak_state_count = std::max<uint64_t>(peak_state_count, states.size());
            if (states.empty()) break;
        }

        terminal_state_count = states.size();
        logical_log_masses.assign(size_t(1) << q_, -std::numeric_limits<double>::infinity());
        if (states.empty() || history_prev.size() != static_cast<size_t>(n_)) {
            decoding = main_decoding;
            used_main_fallback = true;
            syndrome_valid = main_valid;
            refined = false;
            return decoding;
        }

        int best_index = 0;
        for (int index = 0; index < static_cast<int>(states.size()); ++index) {
            const uint64_t label = states[index].key.logical_label;
            logical_log_masses[label] = states[index].log_mass;
            if (states[index].log_mass > states[best_index].log_mass ||
                (states[index].log_mass == states[best_index].log_mass &&
                 better_map_path(states[index], states[best_index]))) {
                best_index = index;
            }
        }
        selected_logical_label = static_cast<int>(states[best_index].key.logical_label);
        decoding.assign(n_, 0);
        int state_index = best_index;
        for (int position = n_ - 1; position >= 0; --position) {
            decoding[order_[position]] = history_bit[position][state_index];
            state_index = history_prev[position][state_index];
        }
        syndrome_valid = valid_syndrome(decoding, syndrome);
        if (!syndrome_valid) throw std::runtime_error("frontier traceback violated syndrome");
        refined = decoding != main_decoding;
        return decoding;
    }

  private:
    struct Key {
        std::vector<uint64_t> parity;
        uint64_t logical_label = 0;

        bool operator<(const Key& other) const {
            if (logical_label != other.logical_label) return logical_label < other.logical_label;
            return parity < other.parity;
        }
        bool operator==(const Key& other) const {
            return logical_label == other.logical_label && parity == other.parity;
        }
    };

    struct State {
        Key key;
        double log_mass = -std::numeric_limits<double>::infinity();
        double best_log_weight = -std::numeric_limits<double>::infinity();
        double bp_tie_score = -std::numeric_limits<double>::infinity();
        int prev_index = -1;
        uint8_t bit = 0;
    };

    int m_ = 0;
    int n_ = 0;
    int q_ = 0;
    std::vector<std::vector<uint8_t>> h_;
    std::vector<std::vector<uint8_t>> l_;
    std::vector<double> priors_;
    std::vector<int> order_;
    std::vector<int> bag_ends_;
    std::vector<uint8_t> bag_end_at_;
    uint64_t max_frontier_states_;
    std::vector<std::vector<int>> closes_at_;
    std::vector<std::vector<int>> active_after_;
    std::vector<std::vector<int>> previous_slots_;
    std::vector<int> first_position_, last_position_;
    std::vector<int> empty_checks_;
    std::vector<uint64_t> logical_masks_;
    std::vector<double> log_probability_zero_, log_probability_one_;

    static double log_add(double left, double right) {
        if (!std::isfinite(left)) return right;
        if (!std::isfinite(right)) return left;
        const double high = std::max(left, right);
        const double low = std::min(left, right);
        return high + std::log1p(std::exp(low - high));
    }

    static double bp_preference(double llr, uint8_t bit) {
        if (std::isnan(llr)) return 0.0;
        const double clipped = std::max(-100.0, std::min(100.0, llr));
        return bit ? -0.5 * clipped : 0.5 * clipped;
    }

    static bool better_map_path(const State& left, const State& right) {
        if (left.best_log_weight != right.best_log_weight) {
            return left.best_log_weight > right.best_log_weight;
        }
        if (left.bp_tie_score != right.bp_tie_score) return left.bp_tie_score > right.bp_tie_score;
        if (left.prev_index != right.prev_index) return left.prev_index < right.prev_index;
        return left.bit < right.bit;
    }

    static bool better_beam_state(const State& left, const State& right, bool bag_boundary) {
        if (bag_boundary && left.log_mass != right.log_mass) return left.log_mass > right.log_mass;
        if (left.best_log_weight != right.best_log_weight) return left.best_log_weight > right.best_log_weight;
        if (!bag_boundary && left.log_mass != right.log_mass) return left.log_mass > right.log_mass;
        if (left.bp_tie_score != right.bp_tie_score) return left.bp_tie_score > right.bp_tie_score;
        return left.key < right.key;
    }

    static bool bit(const std::vector<uint64_t>& words, int index) {
        return (words[index >> 6] >> (index & 63)) & uint64_t(1);
    }

    static void toggle(std::vector<uint64_t>& words, int index) {
        words[index >> 6] ^= uint64_t(1) << (index & 63);
    }

    void validate_order() const {
        if (static_cast<int>(order_.size()) != n_) throw std::invalid_argument("variable order has wrong size");
        std::vector<uint8_t> seen(n_, 0);
        for (int variable : order_) {
            if (variable < 0 || variable >= n_ || seen[variable]) {
                throw std::invalid_argument("variable order must be a permutation");
            }
            seen[variable] = 1;
        }
    }

    void validate_bags() {
        if (n_ == 0) {
            if (!bag_ends_.empty()) throw std::invalid_argument("empty order requires empty bags");
            return;
        }
        if (bag_ends_.empty() || bag_ends_.back() != n_ - 1) {
            throw std::invalid_argument("recursive bags must end at the final variable");
        }
        int previous = -1;
        for (int end : bag_ends_) {
            if (end <= previous || end >= n_) throw std::invalid_argument("invalid recursive bag boundaries");
            peak_recursive_bag_size = std::max(peak_recursive_bag_size, end - previous);
            previous = end;
        }
        recursive_bag_count = static_cast<int>(bag_ends_.size());
    }

    void build_schedule() {
        bag_end_at_.assign(n_, 0);
        for (int end : bag_ends_) bag_end_at_[end] = 1;
        first_position_.assign(m_, n_);
        last_position_.assign(m_, -1);
        std::vector<int> position_of(n_, -1);
        for (int position = 0; position < n_; ++position) position_of[order_[position]] = position;
        for (int check = 0; check < m_; ++check) {
            for (int variable = 0; variable < n_; ++variable) if (h_[check][variable] & 1U) {
                first_position_[check] = std::min(first_position_[check], position_of[variable]);
                last_position_[check] = std::max(last_position_[check], position_of[variable]);
            }
            if (last_position_[check] < 0) empty_checks_.push_back(check);
        }
        closes_at_.assign(n_, {});
        active_after_.assign(n_, {});
        previous_slots_.assign(n_, std::vector<int>(m_, -1));
        for (int check = 0; check < m_; ++check) if (last_position_[check] >= 0) {
            closes_at_[last_position_[check]].push_back(check);
            for (int position = first_position_[check]; position < last_position_[check]; ++position) {
                active_after_[position].push_back(check);
            }
        }
        for (int position = 0; position < n_; ++position) {
            std::sort(active_after_[position].begin(), active_after_[position].end());
            peak_frontier_check_count = std::max<int>(
                peak_frontier_check_count, active_after_[position].size()
            );
            if (position > 0) {
                for (int slot = 0; slot < static_cast<int>(active_after_[position - 1].size()); ++slot) {
                    previous_slots_[position][active_after_[position - 1][slot]] = slot;
                }
            }
        }
        logical_masks_.assign(n_, 0);
        for (int logical = 0; logical < q_; ++logical) {
            for (int variable = 0; variable < n_; ++variable) if (l_[logical][variable] & 1U) {
                logical_masks_[variable] |= uint64_t(1) << logical;
            }
        }
        log_probability_zero_.resize(n_);
        log_probability_one_.resize(n_);
        for (int variable = 0; variable < n_; ++variable) {
            log_probability_zero_[variable] = priors_[variable] == 1.0
                ? -std::numeric_limits<double>::infinity() : std::log1p(-priors_[variable]);
            log_probability_one_[variable] = priors_[variable] == 0.0
                ? -std::numeric_limits<double>::infinity() : std::log(priors_[variable]);
        }
    }

    bool previous_parity(const Key& previous, int position, int check) const {
        if (position == 0) return false;
        const int slot = previous_slots_[position][check];
        return slot >= 0 && bit(previous.parity, slot);
    }

    bool transition_key(const Key& previous, Key& next, int position, int variable,
                        uint8_t value, const std::vector<uint8_t>& syndrome) const {
        next.logical_label = previous.logical_label ^ (value ? logical_masks_[variable] : 0);
        next.parity.assign((active_after_[position].size() + 63) / 64, 0);
        for (int check : closes_at_[position]) {
            bool parity = previous_parity(previous, position, check);
            parity ^= static_cast<bool>((h_[check][variable] & 1U) && value);
            if (parity != static_cast<bool>(syndrome[check] & 1U)) return false;
        }
        for (int slot = 0; slot < static_cast<int>(active_after_[position].size()); ++slot) {
            const int check = active_after_[position][slot];
            bool parity = previous_parity(previous, position, check);
            parity ^= static_cast<bool>((h_[check][variable] & 1U) && value);
            if (parity) toggle(next.parity, slot);
        }
        return true;
    }

    void advance_anchor(Key& key, int position, int variable, uint8_t value,
                        const std::vector<uint8_t>& syndrome) const {
        Key next;
        if (!transition_key(key, next, position, variable, value, syndrome)) {
            throw std::runtime_error("syndrome-valid anchor failed during frontier closure");
        }
        key = std::move(next);
    }

    static void preserve_anchor(std::vector<int>& keep,
                                const std::vector<State>& candidates,
                                const Key& anchor_key) {
        int anchor_index = -1;
        for (int index = 0; index < static_cast<int>(candidates.size()); ++index) {
            if (candidates[index].key == anchor_key) {
                anchor_index = index;
                break;
            }
        }
        if (anchor_index < 0 || std::find(keep.begin(), keep.end(), anchor_index) != keep.end()) return;
        keep.back() = anchor_index;
    }

    bool valid_syndrome(const std::vector<uint8_t>& candidate,
                        const std::vector<uint8_t>& syndrome) const {
        for (int check = 0; check < m_; ++check) {
            uint8_t parity = 0;
            for (int variable = 0; variable < n_; ++variable) {
                parity ^= h_[check][variable] & candidate[variable];
            }
            if (parity != (syndrome[check] & 1U)) return false;
        }
        return true;
    }

    void reset_stats() {
        decoding.clear();
        logical_log_masses.clear();
        expanded_transition_count = 0;
        merged_path_count = 0;
        pruned_state_count = 0;
        closed_check_rejection_count = 0;
        terminal_state_count = 0;
        peak_state_count = 0;
        selected_logical_label = 0;
        truncated = false;
        used_main_fallback = false;
        refined = false;
        syndrome_valid = false;
    }
};

}  // namespace recursive_frontier_bp
}  // namespace ldpc_modified

#endif
