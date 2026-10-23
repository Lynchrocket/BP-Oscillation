#ifndef LDPC_MODIFIED_RECURSIVE_LOOP_CORE_BP_HPP
#define LDPC_MODIFIED_RECURSIVE_LOOP_CORE_BP_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace ldpc_modified {
namespace recursive_loop_core_bp {

class RecursiveLoopCoreRefiner {
  public:
    std::vector<uint8_t> decoding;
    std::vector<int> active_variables;
    std::vector<std::vector<int>> fundamental_circuits;
    std::vector<double> logical_log_masses;
    int tree_variable_count = 0;
    int loop_core_variable_count = 0;
    int loop_component_count = 0;
    int exact_component_count = 0;
    int truncated_component_count = 0;
    uint64_t enumerated_state_count = 0;
    uint64_t feasible_state_count = 0;
    bool refined = false;
    bool syndrome_valid = false;
    int selected_logical_label = 0;

    RecursiveLoopCoreRefiner(
        const std::vector<std::vector<uint8_t>>& pcm,
        const std::vector<std::vector<uint8_t>>& logicals,
        const std::vector<double>& priors,
        const std::vector<int>& check_times,
        const std::vector<int>& check_groups,
        const std::vector<int>& check_layers,
        uint64_t max_factor_states,
        int cluster_mode,
        bool sum_logical_mass)
        : h_(pcm), l_(logicals), priors_(priors), times_(check_times),
          groups_(check_groups), layers_(check_layers),
          max_factor_states_(max_factor_states), cluster_mode_(cluster_mode),
          sum_logical_mass_(sum_logical_mass) {
        m_ = static_cast<int>(h_.size());
        n_ = m_ == 0 ? static_cast<int>(priors_.size()) : static_cast<int>(h_[0].size());
        if (static_cast<int>(priors_.size()) != n_ || max_factor_states_ == 0) {
            throw std::invalid_argument("invalid refiner dimensions or state budget");
        }
        for (double p : priors_) if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
            throw std::invalid_argument("priors must be finite probabilities");
        }
        for (const auto& row : h_) if (static_cast<int>(row.size()) != n_) throw std::invalid_argument("ragged pcm");
        for (const auto& row : l_) if (static_cast<int>(row.size()) != n_) throw std::invalid_argument("ragged logical matrix");
        if (l_.size() > 20) throw std::invalid_argument("at most 20 observables are supported");
        if (static_cast<int>(times_.size()) != m_ || static_cast<int>(groups_.size()) != m_ ||
            static_cast<int>(layers_.size()) != m_) throw std::invalid_argument("invalid check metadata");
        build_graph();
    }

    std::vector<uint8_t>& refine(const std::vector<uint8_t>& syndrome,
                                 const std::vector<uint8_t>& main_decoding,
                                 const std::vector<double>& posterior_llr) {
        if (static_cast<int>(syndrome.size()) != m_ || static_cast<int>(main_decoding.size()) != n_ ||
            static_cast<int>(posterior_llr.size()) != n_) throw std::invalid_argument("invalid decode input dimensions");
        reset_stats();
        decoding = make_valid_seed(syndrome, main_decoding);
        select_active_variables(syndrome, main_decoding, decoding, posterior_llr);
        enumerate_conditioned_assignments(syndrome, decoding);
        refined = decoding != main_decoding;
        syndrome_valid = valid_syndrome(decoding, syndrome);
        if (!syndrome_valid) throw std::runtime_error("refiner returned an invalid syndrome");
        return decoding;
    }

  private:
    int m_ = 0;
    int n_ = 0;
    std::vector<std::vector<uint8_t>> h_;
    std::vector<std::vector<uint8_t>> l_;
    std::vector<double> priors_;
    std::vector<int> times_, groups_, layers_;
    uint64_t max_factor_states_;
    int cluster_mode_;
    bool sum_logical_mass_;
    std::vector<std::vector<int>> check_vars_, var_checks_;
    std::vector<int> core_degree_, cycle_count_, structural_cycle_count_, component_id_;
    std::vector<std::vector<int>> dependency_cells_;

    static double log_add(double a, double b) {
        if (!std::isfinite(a)) return b;
        if (!std::isfinite(b)) return a;
        const double hi = std::max(a, b), lo = std::min(a, b);
        return hi + std::log1p(std::exp(lo - hi));
    }

    void build_graph() {
        check_vars_.assign(m_, {});
        var_checks_.assign(n_, {});
        for (int c = 0; c < m_; ++c) for (int v = 0; v < n_; ++v) if (h_[c][v] & 1U) {
            check_vars_[c].push_back(v);
            var_checks_[v].push_back(c);
        }
        identify_two_core();
        identify_cycles_and_components();
        identify_column_dependency_cells();
    }

    void identify_two_core() {
        std::vector<int> cd(m_), vd(n_);
        std::vector<uint8_t> alive_c(m_, 1), alive_v(n_, 1);
        std::queue<std::pair<bool, int>> leaves;
        for (int c = 0; c < m_; ++c) { cd[c] = check_vars_[c].size(); if (cd[c] < 2) leaves.push({true, c}); }
        for (int v = 0; v < n_; ++v) { vd[v] = var_checks_[v].size(); if (vd[v] < 2) leaves.push({false, v}); }
        while (!leaves.empty()) {
            auto [is_check, id] = leaves.front(); leaves.pop();
            if (is_check) {
                if (!alive_c[id] || cd[id] >= 2) continue;
                alive_c[id] = 0;
                for (int v : check_vars_[id]) if (alive_v[v] && --vd[v] < 2) leaves.push({false, v});
            } else {
                if (!alive_v[id] || vd[id] >= 2) continue;
                alive_v[id] = 0;
                for (int c : var_checks_[id]) if (alive_c[c] && --cd[c] < 2) leaves.push({true, c});
            }
        }
        core_degree_.assign(n_, 0);
        for (int v = 0; v < n_; ++v) if (alive_v[v]) {
            for (int c : var_checks_[v]) core_degree_[v] += alive_c[c];
            if (core_degree_[v] >= 2) ++loop_core_variable_count; else core_degree_[v] = 0;
        }
        tree_variable_count = n_ - loop_core_variable_count;
    }

    bool structural_pair(int a, int b) const {
        const bool temporal = groups_[a] >= 0 && groups_[a] == groups_[b] && times_[a] >= 0 && times_[a] != times_[b];
        const bool spatial = times_[a] >= 0 && times_[a] == times_[b] && groups_[a] >= 0 && groups_[a] != groups_[b];
        const bool recursive = layers_[a] >= 0 && layers_[b] >= 0 && layers_[a] != layers_[b];
        return temporal || spatial || recursive;
    }

    void identify_cycles_and_components() {
        cycle_count_.assign(n_, 0);
        structural_cycle_count_.assign(n_, 0);
        std::vector<std::vector<int>> adjacency(n_);
        for (int a = 0; a < m_; ++a) for (int b = a + 1; b < m_; ++b) {
            std::vector<int> shared;
            std::set_intersection(check_vars_[a].begin(), check_vars_[a].end(), check_vars_[b].begin(),
                                  check_vars_[b].end(), std::back_inserter(shared));
            if (shared.size() < 2) continue;
            const bool structural = structural_pair(a, b);
            for (int v : shared) { ++cycle_count_[v]; if (structural) ++structural_cycle_count_[v]; }
            for (size_t i = 1; i < shared.size(); ++i) {
                adjacency[shared[0]].push_back(shared[i]); adjacency[shared[i]].push_back(shared[0]);
            }
        }
        component_id_.assign(n_, -1);
        for (int root = 0; root < n_; ++root) if (cycle_count_[root] && component_id_[root] < 0) {
            std::queue<int> q; q.push(root); component_id_[root] = loop_component_count;
            while (!q.empty()) { int v = q.front(); q.pop(); for (int u : adjacency[v]) if (component_id_[u] < 0) { component_id_[u] = loop_component_count; q.push(u); } }
            ++loop_component_count;
        }
    }

    void identify_column_dependency_cells() {
        std::set<std::vector<int>> unique_cells;
        std::map<std::vector<int>, std::vector<int>> buckets;
        for (int v = 0; v < n_; ++v) buckets[var_checks_[v]].push_back(v);
        dependency_cells_.clear();
        for (auto& [boundary, variables] : buckets) {
            for (size_t i = 1; i < variables.size(); ++i) {
                unique_cells.insert({variables[0], variables[i]});
            }
        }

        std::vector<std::vector<uint8_t>> rref = h_;
        std::vector<int> pivot_columns;
        int row = 0;
        for (int col = 0; col < n_ && row < m_; ++col) {
            int pivot = row;
            while (pivot < m_ && !rref[pivot][col]) ++pivot;
            if (pivot == m_) continue;
            std::swap(rref[row], rref[pivot]);
            for (int other = 0; other < m_; ++other) if (other != row && rref[other][col]) {
                for (int j = col; j < n_; ++j) rref[other][j] ^= rref[row][j];
            }
            pivot_columns.push_back(col);
            ++row;
        }
        std::vector<uint8_t> is_pivot(n_, 0);
        for (int pivot : pivot_columns) is_pivot[pivot] = 1;
        for (int free_column = 0; free_column < n_; ++free_column) if (!is_pivot[free_column]) {
            std::vector<int> support{free_column};
            for (size_t pivot_row = 0; pivot_row < pivot_columns.size(); ++pivot_row) {
                if (rref[pivot_row][free_column]) support.push_back(pivot_columns[pivot_row]);
            }
            std::sort(support.begin(), support.end());
            if (support.size() >= 2) unique_cells.insert(std::move(support));
        }
        dependency_cells_.assign(unique_cells.begin(), unique_cells.end());
        fundamental_circuits = dependency_cells_;
    }

    void reset_stats() {
        active_variables.clear(); logical_log_masses.clear(); enumerated_state_count = 0;
        feasible_state_count = 0; refined = false; syndrome_valid = false; selected_logical_label = 0;
        exact_component_count = 0; truncated_component_count = 0;
    }

    std::vector<uint8_t> make_valid_seed(const std::vector<uint8_t>& syndrome,
                                         const std::vector<uint8_t>& base) const {
        if (valid_syndrome(base, syndrome)) return base;
        std::vector<std::vector<uint8_t>> a(m_, std::vector<uint8_t>(n_ + 1, 0));
        for (int c = 0; c < m_; ++c) {
            uint8_t current = 0;
            for (int v = 0; v < n_; ++v) { a[c][v] = h_[c][v]; current ^= h_[c][v] & base[v]; }
            a[c][n_] = syndrome[c] ^ current;
        }
        std::vector<int> pivots;
        int row = 0;
        for (int col = 0; col < n_ && row < m_; ++col) {
            int pivot = row; while (pivot < m_ && !a[pivot][col]) ++pivot;
            if (pivot == m_) continue;
            std::swap(a[row], a[pivot]);
            for (int r = 0; r < m_; ++r) if (r != row && a[r][col]) for (int j = col; j <= n_; ++j) a[r][j] ^= a[row][j];
            pivots.push_back(col); ++row;
        }
        for (int r = row; r < m_; ++r) if (a[r][n_]) throw std::runtime_error("inconsistent syndrome");
        std::vector<uint8_t> seed = base;
        for (int r = 0; r < row; ++r) seed[pivots[r]] ^= a[r][n_];
        if (!valid_syndrome(seed, syndrome)) throw std::runtime_error("failed to construct syndrome-valid seed");
        return seed;
    }

    void select_active_variables(const std::vector<uint8_t>& syndrome,
                                 const std::vector<uint8_t>& main,
                                 const std::vector<uint8_t>& seed,
                                 const std::vector<double>& llr) {
        int budget_bits = 0;
        while (budget_bits < 62 && (uint64_t(1) << (budget_bits + 1)) <= max_factor_states_) ++budget_bits;
        budget_bits = std::min(budget_bits, n_);
        std::vector<uint8_t> touches_residual(n_, 0), repaired(n_, 0);
        for (int c = 0; c < m_; ++c) {
            uint8_t residual = syndrome[c];
            for (int v : check_vars_[c]) residual ^= main[v];
            if (residual) for (int v : check_vars_[c]) touches_residual[v] = 1;
        }
        for (int v = 0; v < n_; ++v) repaired[v] = main[v] != seed[v];
        std::vector<uint8_t> selected(n_, 0);
        if (cluster_mode_ >= 1) {
            std::vector<int> cell_order(dependency_cells_.size());
            std::iota(cell_order.begin(), cell_order.end(), 0);
            std::stable_sort(cell_order.begin(), cell_order.end(), [&](int a, int b) {
                auto key = [&](int cell) {
                    int seed_support = 0, residual_support = 0, structural = 0, cycles = 0;
                    double ambiguity = 0.0;
                    for (int v : dependency_cells_[cell]) {
                        seed_support += seed[v];
                        residual_support += touches_residual[v];
                        structural += structural_cycle_count_[v];
                        cycles += cycle_count_[v];
                        ambiguity -= std::abs(llr[v]);
                    }
                    if (cluster_mode_ < 2) structural = 0;
                    return std::tuple<int,int,int,int,int,double,int>(
                        seed_support, residual_support, structural, cycles,
                        -static_cast<int>(dependency_cells_[cell].size()), ambiguity, -cell);
                };
                return key(a) > key(b);
            });
            for (int cell : cell_order) {
                const auto& variables = dependency_cells_[cell];
                size_t new_count = 0;
                for (int v : variables) new_count += !selected[v];
                if (active_variables.size() + new_count > static_cast<size_t>(budget_bits)) continue;
                for (int v : variables) if (!selected[v]) {
                    active_variables.push_back(v);
                    selected[v] = 1;
                }
            }
        }
        std::vector<int> order(n_); std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            auto key = [&](int v) {
                const int structural = cluster_mode_ >= 2 ? structural_cycle_count_[v] : 0;
                const int cycles = cluster_mode_ >= 1 ? cycle_count_[v] : 0;
                const int core = cluster_mode_ >= 1 ? core_degree_[v] : 0;
                return std::tuple<int,int,int,int,int,double,int>(repaired[v], touches_residual[v], structural, cycles, core, -std::abs(llr[v]), -v);
            };
            return key(a) > key(b);
        });
        if (cluster_mode_ == 0) {
            for (int v : order) if (!selected[v] && active_variables.size() < static_cast<size_t>(budget_bits)) {
                active_variables.push_back(v);
                selected[v] = 1;
            }
        }
        std::sort(active_variables.begin(), active_variables.end());
        std::set<int> components;
        for (int v : active_variables) if (component_id_[v] >= 0) components.insert(component_id_[v]);
        exact_component_count = components.size();
        std::set<int> all_components;
        for (int id : component_id_) if (id >= 0) all_components.insert(id);
        truncated_component_count = static_cast<int>(all_components.size() - components.size());
    }

    int logical_label(const std::vector<uint8_t>& x) const {
        if (l_.size() > 20) throw std::invalid_argument("at most 20 observables are supported");
        int label = 0;
        for (size_t r = 0; r < l_.size(); ++r) {
            uint8_t bit = 0; for (int v = 0; v < n_; ++v) bit ^= l_[r][v] & x[v];
            label |= static_cast<int>(bit) << r;
        }
        return label;
    }

    bool valid_syndrome(const std::vector<uint8_t>& x, const std::vector<uint8_t>& syndrome) const {
        for (int c = 0; c < m_; ++c) { uint8_t bit = 0; for (int v : check_vars_[c]) bit ^= x[v]; if (bit != syndrome[c]) return false; }
        return true;
    }

    void enumerate_conditioned_assignments(const std::vector<uint8_t>& syndrome, const std::vector<uint8_t>& base) {
        const size_t k = active_variables.size();
        if (k >= 63) throw std::invalid_argument("active region exceeds enumeration width");
        const uint64_t states = uint64_t(1) << k;
        std::map<int, double> masses, best_weight;
        std::map<int, std::vector<uint8_t>> best_state;
        std::vector<uint8_t> candidate = base;
        for (uint64_t state = 0; state < states; ++state) {
            ++enumerated_state_count;
            double log_weight = 0.0;
            for (size_t i = 0; i < k; ++i) {
                const int v = active_variables[i]; candidate[v] = (state >> i) & 1U;
                const double p = priors_[v];
                if ((candidate[v] && p == 0.0) || (!candidate[v] && p == 1.0)) {
                    log_weight = -std::numeric_limits<double>::infinity();
                } else if (std::isfinite(log_weight)) {
                    log_weight += candidate[v] ? std::log(p) : std::log1p(-p);
                }
            }
            if (!valid_syndrome(candidate, syndrome)) continue;
            ++feasible_state_count;
            const int label = logical_label(candidate);
            masses[label] = masses.count(label) ? log_add(masses[label], log_weight) : log_weight;
            if (!best_weight.count(label) || log_weight > best_weight[label]) { best_weight[label] = log_weight; best_state[label] = candidate; }
        }
        const int label_count = 1 << l_.size();
        logical_log_masses.assign(label_count, -std::numeric_limits<double>::infinity());
        for (const auto& [label, mass] : masses) logical_log_masses[label] = mass;
        if (masses.empty()) return;
        auto score = [&](int label) { return sum_logical_mass_ ? masses[label] : best_weight[label]; };
        selected_logical_label = masses.begin()->first;
        for (const auto& [label, mass] : masses) if (score(label) > score(selected_logical_label)) selected_logical_label = label;
        decoding = best_state[selected_logical_label];
    }
};

}  // namespace recursive_loop_core_bp
}  // namespace ldpc_modified

#endif
