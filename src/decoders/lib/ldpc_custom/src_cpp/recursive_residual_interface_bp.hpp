#ifndef LDPC_MODIFIED_RECURSIVE_RESIDUAL_INTERFACE_BP_HPP
#define LDPC_MODIFIED_RECURSIVE_RESIDUAL_INTERFACE_BP_HPP

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
namespace recursive_residual_interface_bp {

class RecursiveResidualInterfaceRefiner {
  public:
    std::vector<uint8_t> decoding;
    std::vector<double> logical_log_masses;
    uint64_t expanded_pair_count = 0;
    uint64_t merged_path_count = 0;
    uint64_t pruned_state_count = 0;
    uint64_t closed_check_rejection_count = 0;
    uint64_t terminal_state_count = 0;
    uint64_t peak_state_count = 0;
    int region_count = 0;
    int compressed_region_count = 0;
    int peak_boundary_check_count = 0;
    int selected_logical_label = 0;
    double residual_spectral_energy = 0.0;
    double mean_effective_walsh_rank = 0.0;
    bool truncated = false;
    bool used_main_fallback = false;
    bool refined = false;
    bool syndrome_valid = false;

    RecursiveResidualInterfaceRefiner(
        const std::vector<std::vector<uint8_t>>& pcm,
        const std::vector<std::vector<uint8_t>>& logicals,
        const std::vector<double>& priors,
        const std::vector<int>& variable_order,
        uint64_t max_interface_states,
        bool residual_compression)
        : h_(pcm), l_(logicals), priors_(priors), order_(variable_order),
          max_interface_states_(max_interface_states),
          residual_compression_(residual_compression) {
        m_ = static_cast<int>(h_.size());
        n_ = m_ == 0 ? static_cast<int>(priors_.size()) : static_cast<int>(h_[0].size());
        q_ = static_cast<int>(l_.size());
        if (n_ <= 0 || static_cast<int>(priors_.size()) != n_ ||
            max_interface_states_ == 0 ||
            max_interface_states_ > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("invalid interface dimensions or state budget");
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
        check_degree_.assign(m_, 0);
        for (int c = 0; c < m_; ++c) {
            check_degree_[c] = static_cast<int>(std::count(h_[c].begin(), h_[c].end(), uint8_t{1}));
            if (check_degree_[c] == 0) empty_checks_.push_back(c);
        }
        log_probability_zero_.resize(n_);
        log_probability_one_.resize(n_);
        for (int v = 0; v < n_; ++v) {
            log_probability_zero_[v] = priors_[v] == 1.0 ? neg_inf() : std::log1p(-priors_[v]);
            log_probability_one_[v] = priors_[v] == 0.0 ? neg_inf() : std::log(priors_[v]);
        }
        root_ = build_region(0, n_);
        region_count = static_cast<int>(nodes_.size());
    }

    std::vector<uint8_t>& refine(const std::vector<uint8_t>& syndrome,
                                 const std::vector<uint8_t>& main_decoding,
                                 const std::vector<double>& posterior_llr) {
        if (static_cast<int>(syndrome.size()) != m_ ||
            static_cast<int>(main_decoding.size()) != n_ ||
            static_cast<int>(posterior_llr.size()) != n_) {
            throw std::invalid_argument("invalid interface decode input dimensions");
        }
        reset_stats();
        const bool main_valid = valid_syndrome(main_decoding, syndrome);
        for (int c : empty_checks_) {
            if (syndrome[c]) return fallback(main_decoding, false);
        }
        results_.assign(nodes_.size(), {});
        anchor_keys_.assign(nodes_.size(), {});
        solve_region(root_, syndrome, main_decoding, posterior_llr, main_valid);
        auto& root_states = results_[root_];
        terminal_state_count = root_states.size();
        logical_log_masses.assign(size_t{1} << q_, neg_inf());
        if (root_states.empty()) return fallback(main_decoding, main_valid);

        int best = 0;
        for (int i = 0; i < static_cast<int>(root_states.size()); ++i) {
            const uint64_t label = root_states[i].key.logical;
            logical_log_masses[label] = root_states[i].log_mass;
            if (better_path(root_states[i], root_states[best])) best = i;
        }
        selected_logical_label = static_cast<int>(root_states[best].key.logical);
        decoding.assign(n_, 0);
        traceback(root_, best);
        syndrome_valid = valid_syndrome(decoding, syndrome);
        if (!syndrome_valid && main_valid) return fallback(main_decoding, true);
        refined = syndrome_valid;
        mean_effective_walsh_rank = walsh_region_count_ == 0
            ? 0.0 : walsh_rank_sum_ / static_cast<double>(walsh_region_count_);
        return decoding;
    }

  private:
    struct Key {
        std::vector<uint64_t> parity;
        uint64_t logical = 0;
        bool operator<(const Key& other) const {
            if (logical != other.logical) return logical < other.logical;
            return parity < other.parity;
        }
        bool operator==(const Key& other) const {
            return logical == other.logical && parity == other.parity;
        }
    };

    struct State {
        Key key;
        double log_mass = neg_inf();
        double log_bp_mass = neg_inf();
        double best_log_weight = neg_inf();
        double best_bp_weight = neg_inf();
        double compression_score = 0.0;
        double group_priority = 0.0;
        uint64_t compression_group = 0;
        int left_state = -1;
        int right_state = -1;
        uint8_t leaf_bit = 0;
    };

    struct Region {
        int left = -1;
        int right = -1;
        int variable = -1;
        std::vector<int> inside_count;
        std::vector<int> boundary;
        std::vector<int> closed_here;
        std::vector<int> left_boundary_slot;
        std::vector<int> right_boundary_slot;
        std::vector<int> left_closed_slot;
        std::vector<int> right_closed_slot;
    };

    std::vector<std::vector<uint8_t>> h_;
    std::vector<std::vector<uint8_t>> l_;
    std::vector<double> priors_;
    std::vector<int> order_;
    uint64_t max_interface_states_;
    bool residual_compression_;
    int m_ = 0;
    int n_ = 0;
    int q_ = 0;
    int root_ = -1;
    std::vector<int> check_degree_;
    std::vector<int> empty_checks_;
    std::vector<double> log_probability_zero_;
    std::vector<double> log_probability_one_;
    std::vector<Region> nodes_;
    std::vector<std::vector<State>> results_;
    std::vector<Key> anchor_keys_;
    double walsh_rank_sum_ = 0.0;
    int walsh_region_count_ = 0;

    static double neg_inf() { return -std::numeric_limits<double>::infinity(); }

    static double log_add(double a, double b) {
        if (!std::isfinite(a)) return b;
        if (!std::isfinite(b)) return a;
        const double hi = std::max(a, b);
        return hi + std::log1p(std::exp(std::min(a, b) - hi));
    }

    static bool get_bit(const std::vector<uint64_t>& words, int bit) {
        return (words[bit >> 6] >> (bit & 63)) & uint64_t{1};
    }

    static void set_bit(std::vector<uint64_t>& words, int bit, bool value) {
        if (value) words[bit >> 6] |= uint64_t{1} << (bit & 63);
    }

    static double bp_log_probability(double llr, uint8_t bit) {
        if (std::isnan(llr)) return -std::log(2.0);
        if (bit) {
            if (llr >= 0.0) return -llr - std::log1p(std::exp(-llr));
            return -std::log1p(std::exp(llr));
        }
        if (llr >= 0.0) return -std::log1p(std::exp(-llr));
        return llr - std::log1p(std::exp(llr));
    }

    void validate_order() const {
        if (static_cast<int>(order_.size()) != n_) throw std::invalid_argument("order must contain every variable");
        std::vector<uint8_t> seen(n_, 0);
        for (int v : order_) {
            if (v < 0 || v >= n_ || seen[v]) throw std::invalid_argument("order must be a permutation");
            seen[v] = 1;
        }
    }

    int build_region(int begin, int end) {
        Region node;
        if (end - begin == 1) {
            node.variable = order_[begin];
            node.inside_count.assign(m_, 0);
            for (int c = 0; c < m_; ++c) node.inside_count[c] = h_[c][node.variable] ? 1 : 0;
        } else {
            const int middle = begin + (end - begin) / 2;
            node.left = build_region(begin, middle);
            node.right = build_region(middle, end);
            node.inside_count.resize(m_);
            for (int c = 0; c < m_; ++c) {
                node.inside_count[c] = nodes_[node.left].inside_count[c] + nodes_[node.right].inside_count[c];
            }
        }
        for (int c = 0; c < m_; ++c) {
            if (node.inside_count[c] > 0 && node.inside_count[c] < check_degree_[c]) node.boundary.push_back(c);
        }
        peak_boundary_check_count = std::max(peak_boundary_check_count, static_cast<int>(node.boundary.size()));
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back(std::move(node));
        if (nodes_[id].left >= 0) build_relations(id);
        return id;
    }

    static std::vector<int> slots(const std::vector<int>& source, const std::vector<int>& checks) {
        std::vector<int> result;
        result.reserve(checks.size());
        for (int check : checks) {
            auto it = std::lower_bound(source.begin(), source.end(), check);
            result.push_back(it != source.end() && *it == check ? static_cast<int>(it - source.begin()) : -1);
        }
        return result;
    }

    void build_relations(int id) {
        Region& node = nodes_[id];
        const Region& left = nodes_[node.left];
        const Region& right = nodes_[node.right];
        for (int c = 0; c < m_; ++c) {
            if (node.inside_count[c] == check_degree_[c] &&
                left.inside_count[c] < check_degree_[c] && right.inside_count[c] < check_degree_[c]) {
                node.closed_here.push_back(c);
            }
        }
        node.left_boundary_slot = slots(left.boundary, node.boundary);
        node.right_boundary_slot = slots(right.boundary, node.boundary);
        node.left_closed_slot = slots(left.boundary, node.closed_here);
        node.right_closed_slot = slots(right.boundary, node.closed_here);
    }

    Key leaf_key(int node_id, uint8_t bit) const {
        const Region& node = nodes_[node_id];
        Key key;
        key.parity.assign((node.boundary.size() + 63) / 64, 0);
        if (bit) {
            for (int i = 0; i < static_cast<int>(node.boundary.size()); ++i) set_bit(key.parity, i, true);
            for (int logical = 0; logical < q_; ++logical) {
                if (l_[logical][node.variable]) key.logical |= uint64_t{1} << logical;
            }
        }
        return key;
    }

    bool leaf_internal_valid(int node_id, uint8_t bit, const std::vector<uint8_t>& syndrome) const {
        const int variable = nodes_[node_id].variable;
        for (int c = 0; c < m_; ++c) {
            if (h_[c][variable] && check_degree_[c] == 1 && bit != syndrome[c]) return false;
        }
        return true;
    }

    Key combine_key(const Region& node, const State& left, const State& right) const {
        Key key;
        key.logical = left.key.logical ^ right.key.logical;
        key.parity.assign((node.boundary.size() + 63) / 64, 0);
        for (int i = 0; i < static_cast<int>(node.boundary.size()); ++i) {
            bool value = false;
            if (node.left_boundary_slot[i] >= 0) value ^= get_bit(left.key.parity, node.left_boundary_slot[i]);
            if (node.right_boundary_slot[i] >= 0) value ^= get_bit(right.key.parity, node.right_boundary_slot[i]);
            set_bit(key.parity, i, value);
        }
        return key;
    }

    bool closed_valid(const Region& node, const State& left, const State& right,
                      const std::vector<uint8_t>& syndrome) const {
        for (int i = 0; i < static_cast<int>(node.closed_here.size()); ++i) {
            bool value = false;
            if (node.left_closed_slot[i] >= 0) value ^= get_bit(left.key.parity, node.left_closed_slot[i]);
            if (node.right_closed_slot[i] >= 0) value ^= get_bit(right.key.parity, node.right_closed_slot[i]);
            if (value != static_cast<bool>(syndrome[node.closed_here[i]])) return false;
        }
        return true;
    }

    static bool better_path(const State& a, const State& b) {
        if (a.log_mass != b.log_mass) return a.log_mass > b.log_mass;
        if (a.best_log_weight != b.best_log_weight) return a.best_log_weight > b.best_log_weight;
        return a.best_bp_weight > b.best_bp_weight;
    }

    void merge_candidate(std::map<Key, int>& index, std::vector<State>& candidates, State candidate) {
        auto [entry, inserted] = index.emplace(candidate.key, static_cast<int>(candidates.size()));
        if (inserted) {
            candidates.push_back(std::move(candidate));
            return;
        }
        ++merged_path_count;
        State& current = candidates[entry->second];
        current.log_mass = log_add(current.log_mass, candidate.log_mass);
        current.log_bp_mass = log_add(current.log_bp_mass, candidate.log_bp_mass);
        if (candidate.best_log_weight > current.best_log_weight ||
            (candidate.best_log_weight == current.best_log_weight && candidate.best_bp_weight > current.best_bp_weight)) {
            const double channel_mass = current.log_mass;
            const double bp_mass = current.log_bp_mass;
            current = std::move(candidate);
            current.log_mass = channel_mass;
            current.log_bp_mass = bp_mass;
        }
    }

    Key make_anchor(int node_id, const std::vector<uint8_t>& main) {
        const Region& node = nodes_[node_id];
        if (node.variable >= 0) return leaf_key(node_id, main[node.variable]);
        State left;
        State right;
        left.key = anchor_keys_[node.left];
        right.key = anchor_keys_[node.right];
        return combine_key(node, left, right);
    }

    void solve_region(int node_id, const std::vector<uint8_t>& syndrome,
                      const std::vector<uint8_t>& main, const std::vector<double>& posterior_llr,
                      bool main_valid) {
        const Region& node = nodes_[node_id];
        std::vector<State> candidates;
        std::map<Key, int> merged;
        if (node.variable >= 0) {
            for (uint8_t bit = 0; bit <= 1; ++bit) {
                if (!leaf_internal_valid(node_id, bit, syndrome)) {
                    ++closed_check_rejection_count;
                    continue;
                }
                State state;
                state.key = leaf_key(node_id, bit);
                state.log_mass = bit ? log_probability_one_[node.variable] : log_probability_zero_[node.variable];
                state.log_bp_mass = bp_log_probability(posterior_llr[node.variable], bit);
                state.best_log_weight = state.log_mass;
                state.best_bp_weight = state.log_bp_mass;
                state.leaf_bit = bit;
                if (std::isfinite(state.log_mass)) merge_candidate(merged, candidates, std::move(state));
            }
        } else {
            solve_region(node.left, syndrome, main, posterior_llr, main_valid);
            solve_region(node.right, syndrome, main, posterior_llr, main_valid);
            const auto& left_states = results_[node.left];
            const auto& right_states = results_[node.right];
            for (int li = 0; li < static_cast<int>(left_states.size()); ++li) {
                for (int ri = 0; ri < static_cast<int>(right_states.size()); ++ri) {
                    ++expanded_pair_count;
                    if (!closed_valid(node, left_states[li], right_states[ri], syndrome)) {
                        ++closed_check_rejection_count;
                        continue;
                    }
                    State state;
                    state.key = combine_key(node, left_states[li], right_states[ri]);
                    state.log_mass = left_states[li].log_mass + right_states[ri].log_mass;
                    state.log_bp_mass = left_states[li].log_bp_mass + right_states[ri].log_bp_mass;
                    state.best_log_weight = left_states[li].best_log_weight + right_states[ri].best_log_weight;
                    state.best_bp_weight = left_states[li].best_bp_weight + right_states[ri].best_bp_weight;
                    state.left_state = li;
                    state.right_state = ri;
                    merge_candidate(merged, candidates, std::move(state));
                }
            }
        }
        if (main_valid) anchor_keys_[node_id] = make_anchor(node_id, main);
        compress(node_id, candidates, main_valid);
        peak_state_count = std::max<uint64_t>(peak_state_count, results_[node_id].size());
    }

    static void fwht(std::vector<double>& values) {
        for (size_t length = 1; length < values.size(); length <<= 1) {
            for (size_t start = 0; start < values.size(); start += length << 1) {
                for (size_t i = 0; i < length; ++i) {
                    const double a = values[start + i];
                    const double b = values[start + length + i];
                    values[start + i] = a + b;
                    values[start + length + i] = a - b;
                }
            }
        }
    }

    uint64_t projected_signature(const State& state, const std::vector<int>& dimensions,
                                 int boundary_size) const {
        uint64_t signature = 0;
        for (int i = 0; i < static_cast<int>(dimensions.size()); ++i) {
            const int dimension = dimensions[i];
            const bool bit = dimension < boundary_size
                ? get_bit(state.key.parity, dimension)
                : ((state.key.logical >> (dimension - boundary_size)) & uint64_t{1});
            if (bit) signature |= uint64_t{1} << i;
        }
        return signature;
    }

    void residual_scores(std::vector<State>& candidates, int boundary_size) {
        const double max_channel = std::max_element(candidates.begin(), candidates.end(),
            [](const State& a, const State& b) { return a.log_mass < b.log_mass; })->log_mass;
        const double max_bp = std::max_element(candidates.begin(), candidates.end(),
            [](const State& a, const State& b) { return a.log_bp_mass < b.log_bp_mass; })->log_bp_mass;
        std::vector<double> p(candidates.size()), q(candidates.size());
        double p_sum = 0.0, q_sum = 0.0;
        for (size_t i = 0; i < candidates.size(); ++i) {
            p[i] = std::exp(candidates[i].log_mass - max_channel);
            q[i] = std::isfinite(max_bp)
                ? std::exp(candidates[i].log_bp_mass - max_bp)
                : 1.0;
            p_sum += p[i];
            q_sum += q[i];
        }
        for (size_t i = 0; i < candidates.size(); ++i) { p[i] /= p_sum; q[i] /= q_sum; }

        int projected_bits = 0;
        while (projected_bits < 16 && projected_bits < boundary_size + q_ &&
               (size_t{1} << (projected_bits + 1)) <= max_interface_states_) ++projected_bits;
        std::vector<std::pair<double, int>> entropy_dimensions;
        for (int dimension = 0; dimension < boundary_size + q_; ++dimension) {
            double probability_one = 0.0;
            for (size_t i = 0; i < candidates.size(); ++i) {
                const bool bit = dimension < boundary_size
                    ? get_bit(candidates[i].key.parity, dimension)
                    : ((candidates[i].key.logical >> (dimension - boundary_size)) & uint64_t{1});
                if (bit) probability_one += q[i];
            }
            entropy_dimensions.emplace_back(4.0 * probability_one * (1.0 - probability_one), dimension);
        }
        std::stable_sort(entropy_dimensions.begin(), entropy_dimensions.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });
        std::vector<int> dimensions;
        for (int i = 0; i < projected_bits; ++i) dimensions.push_back(entropy_dimensions[i].second);
        const size_t spectrum_size = size_t{1} << projected_bits;
        std::vector<double> residual(spectrum_size, 0.0);
        for (size_t i = 0; i < candidates.size(); ++i) {
            const auto signature = projected_signature(candidates[i], dimensions, boundary_size);
            residual[signature] += p[i] - q[i];
        }
        fwht(residual);
        const double normalization = static_cast<double>(spectrum_size);
        std::vector<std::pair<double, int>> modes;
        double energy = 0.0;
        for (int i = 0; i < static_cast<int>(residual.size()); ++i) {
            residual[i] /= normalization;
            const double mode_energy = residual[i] * residual[i];
            modes.emplace_back(mode_energy, i);
            energy += mode_energy;
        }
        residual_spectral_energy += energy;
        std::stable_sort(modes.begin(), modes.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
        double retained = 0.0;
        int rank = 0;
        const double target = 0.99 * energy;
        while (rank < static_cast<int>(modes.size()) && retained < target) retained += modes[rank++].first;
        walsh_rank_sum_ += rank;
        ++walsh_region_count_;
        std::vector<double> retained_spectrum(spectrum_size, 0.0);
        for (int i = 0; i < rank; ++i) retained_spectrum[modes[i].second] = residual[modes[i].second];
        // Coefficients were normalized by N above, so one unnormalized FWHT is the inverse.
        fwht(retained_spectrum);
        for (size_t i = 0; i < candidates.size(); ++i) {
            const auto signature = projected_signature(candidates[i], dimensions, boundary_size);
            candidates[i].compression_score = p[i];
            candidates[i].group_priority = std::abs(retained_spectrum[signature]);
            candidates[i].compression_group = signature;
        }
    }

    void compress(int node_id, std::vector<State>& candidates, bool main_valid) {
        if (candidates.size() <= max_interface_states_) {
            results_[node_id] = std::move(candidates);
            return;
        }
        truncated = true;
        ++compressed_region_count;
        pruned_state_count += candidates.size() - max_interface_states_;
        if (residual_compression_) {
            residual_scores(candidates, static_cast<int>(nodes_[node_id].boundary.size()));
        } else {
            const double max_mass = std::max_element(candidates.begin(), candidates.end(),
                [](const State& a, const State& b) { return a.log_mass < b.log_mass; })->log_mass;
            for (State& state : candidates) state.compression_score = std::exp(state.log_mass - max_mass);
        }
        std::vector<int> ranked(candidates.size());
        std::iota(ranked.begin(), ranked.end(), 0);
        std::stable_sort(ranked.begin(), ranked.end(), [&](int a, int b) {
            if (candidates[a].compression_score != candidates[b].compression_score)
                return candidates[a].compression_score > candidates[b].compression_score;
            return better_path(candidates[a], candidates[b]);
        });
        std::map<uint64_t, int> logical_best;
        for (int index : ranked) logical_best.emplace(candidates[index].key.logical, index);
        std::vector<int> logical_representatives;
        for (const auto& [label, index] : logical_best) logical_representatives.push_back(index);
        std::stable_sort(logical_representatives.begin(), logical_representatives.end(), [&](int a, int b) {
            if (candidates[a].compression_score != candidates[b].compression_score)
                return candidates[a].compression_score > candidates[b].compression_score;
            return better_path(candidates[a], candidates[b]);
        });
        std::vector<int> order;
        order.reserve(static_cast<size_t>(max_interface_states_));
        for (int index : logical_representatives) {
            if (order.size() == max_interface_states_) break;
            order.push_back(index);
        }
        if (residual_compression_ && order.size() < max_interface_states_) {
            std::map<uint64_t, int> group_best;
            for (int index : ranked) group_best.emplace(candidates[index].compression_group, index);
            std::vector<int> group_representatives;
            for (const auto& [group, index] : group_best) group_representatives.push_back(index);
            std::stable_sort(group_representatives.begin(), group_representatives.end(), [&](int a, int b) {
                if (candidates[a].group_priority != candidates[b].group_priority)
                    return candidates[a].group_priority > candidates[b].group_priority;
                return better_path(candidates[a], candidates[b]);
            });
            const size_t intervention_budget = static_cast<size_t>(std::sqrt(max_interface_states_));
            size_t interventions = 0;
            for (int index : group_representatives) {
                if (order.size() == max_interface_states_ || interventions == intervention_budget) break;
                if (std::find(order.begin(), order.end(), index) == order.end()) {
                    order.push_back(index);
                    ++interventions;
                }
            }
        }
        for (int index : ranked) {
            if (order.size() == max_interface_states_) break;
            if (std::find(order.begin(), order.end(), index) == order.end()) order.push_back(index);
        }
        if (main_valid) {
            auto anchor = std::find_if(candidates.begin(), candidates.end(), [&](const State& state) {
                return state.key == anchor_keys_[node_id];
            });
            if (anchor != candidates.end()) {
                const int anchor_index = static_cast<int>(anchor - candidates.begin());
                if (std::find(order.begin(), order.end(), anchor_index) == order.end()) {
                    auto replace = std::find_if(order.begin(), order.end(), [&](int index) {
                        return candidates[index].key.logical == candidates[anchor_index].key.logical;
                    });
                    if (replace == order.end()) replace = std::prev(order.end());
                    *replace = anchor_index;
                }
            }
        }
        results_[node_id].reserve(order.size());
        for (int index : order) results_[node_id].push_back(std::move(candidates[index]));
    }

    void traceback(int node_id, int state_index) {
        const Region& node = nodes_[node_id];
        const State& state = results_[node_id][state_index];
        if (node.variable >= 0) {
            decoding[node.variable] = state.leaf_bit;
            return;
        }
        traceback(node.left, state.left_state);
        traceback(node.right, state.right_state);
    }

    bool valid_syndrome(const std::vector<uint8_t>& error, const std::vector<uint8_t>& syndrome) const {
        for (int c = 0; c < m_; ++c) {
            uint8_t parity = 0;
            for (int v = 0; v < n_; ++v) parity ^= h_[c][v] & error[v];
            if (parity != syndrome[c]) return false;
        }
        return true;
    }

    std::vector<uint8_t>& fallback(const std::vector<uint8_t>& main, bool main_valid) {
        decoding = main;
        used_main_fallback = true;
        syndrome_valid = main_valid;
        refined = false;
        mean_effective_walsh_rank = walsh_region_count_ == 0
            ? 0.0 : walsh_rank_sum_ / static_cast<double>(walsh_region_count_);
        return decoding;
    }

    void reset_stats() {
        logical_log_masses.clear();
        expanded_pair_count = merged_path_count = pruned_state_count = 0;
        closed_check_rejection_count = terminal_state_count = peak_state_count = 0;
        compressed_region_count = selected_logical_label = 0;
        residual_spectral_energy = mean_effective_walsh_rank = 0.0;
        walsh_rank_sum_ = 0.0;
        walsh_region_count_ = 0;
        truncated = used_main_fallback = refined = syndrome_valid = false;
    }
};

}  // namespace recursive_residual_interface_bp
}  // namespace ldpc_modified

#endif
