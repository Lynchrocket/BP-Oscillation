#ifndef LDPC_MODIFIED_MBBP_LD_H
#define LDPC_MODIFIED_MBBP_LD_H

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "bp.hpp"

namespace ldpc_modified {
namespace mbbp_ld {

using BpEntry = ldpc_modified::bp::BpEntry;
using BpSparse = ldpc_modified::bp::BpSparse;

class MbbpLdDecoder {
public:
    BpSparse &pcm;
    int check_count;
    int bit_count;
    bool converge;
    std::vector<uint8_t> decoding;
    int candidate_count;
    int unique_candidate_count;
    int selected_frequency;
    int selected_weight;

    MbbpLdDecoder(
            BpSparse &parity_check_matrix,
            const std::vector<double> &channel_probabilities,
            int maximum_iterations,
            int bp_method,
            int schedule,
            double min_sum_scaling_factor,
            const std::vector<std::vector<int>> &check_permutations) :
            pcm(parity_check_matrix),
            check_count(parity_check_matrix.m),
            bit_count(parity_check_matrix.n),
            converge(false),
            decoding(bit_count, 0),
            candidate_count(0),
            unique_candidate_count(0),
            selected_frequency(0),
            selected_weight(0) {
        if (channel_probabilities.size() != static_cast<std::size_t>(bit_count)) {
            throw std::runtime_error(
                    "Channel probabilities must have length equal to the number of bits");
        }
        if (maximum_iterations < 0) {
            throw std::runtime_error("maximum_iterations must be non-negative");
        }
        if (bp_method != 0 && bp_method != 1) {
            throw std::runtime_error("Unsupported BP method");
        }
        if (schedule != 0 && schedule != 1) {
            throw std::runtime_error("Unsupported BP schedule");
        }

        copy_original_graph();
        validate_permutations(check_permutations);
        for (const auto &permutation: check_permutations) {
            partitions.push_back(build_maximal_subtrees(permutation));
        }
        build_bases(
                channel_probabilities,
                maximum_iterations,
                static_cast<ldpc_modified::bp::BpMethod>(bp_method),
                static_cast<ldpc_modified::bp::BpSchedule>(schedule),
                min_sum_scaling_factor);
        start_workers();
    }

    ~MbbpLdDecoder() {
        shutdown_workers();
    }

    std::vector<uint8_t> decode(std::vector<uint8_t> &syndrome) {
        std::lock_guard<std::mutex> decode_lock(decode_mutex);
        if (syndrome.size() != static_cast<std::size_t>(check_count)) {
            throw std::runtime_error("Syndrome length must equal the original check count");
        }

        reset_decode_state();
        run_basis_decoders(syndrome);

        for (std::size_t basis_index = 0; basis_index < basis_decoders.size(); basis_index++) {
            if (basis_converged[basis_index]) {
                converged_candidates.push_back(basis_candidates[basis_index]);
            }
        }

        candidate_count = checked_size(converged_candidates.size());
        if (converged_candidates.empty()) {
            return decoding;
        }

        aggregate_candidates();
        converge = true;
        return decoding;
    }

    int partition_count() const {
        return checked_size(partitions.size());
    }

    int partition_subtree_count(int partition_index) const {
        return checked_size(partitions.at(partition_index).size());
    }

    int partition_subtree_size(int partition_index, int subtree_index) const {
        return checked_size(partitions.at(partition_index).at(subtree_index).size());
    }

    int partition_subtree_row(int partition_index, int subtree_index, int row_index) const {
        return partitions.at(partition_index).at(subtree_index).at(row_index);
    }

    int basis_count() const {
        return checked_size(basis_source_rows.size());
    }

    int basis_row_count(int basis_index) const {
        return checked_size(basis_source_rows.at(basis_index).size());
    }

    int basis_source_row(int basis_index, int row_index) const {
        return basis_source_rows.at(basis_index).at(row_index);
    }

    int basis_candidate_bit(int basis_index, int bit_index) const {
        return basis_candidates.at(basis_index).at(bit_index);
    }

    bool basis_converged_at(int basis_index) const {
        return basis_converged.at(basis_index) != 0;
    }

    int basis_iterations_at(int basis_index) const {
        return basis_decoders.at(basis_index)->iterations;
    }

    int basis_edge_count_at(int basis_index) const {
        return basis_pcms.at(basis_index)->entry_count();
    }

    int converged_candidate_bit(int candidate_index, int bit_index) const {
        return converged_candidates.at(candidate_index).at(bit_index);
    }

    int unique_candidate_frequency_at(int candidate_index) const {
        return unique_frequencies.at(candidate_index);
    }

    int unique_candidate_weight_at(int candidate_index) const {
        return hamming_weight(unique_candidates.at(candidate_index));
    }

    int unique_candidate_bit(int candidate_index, int bit_index) const {
        return unique_candidates.at(candidate_index).at(bit_index);
    }

    double channel_probability(int bit_index) const {
        return basis_decoders.empty()
               ? 0.0
               : basis_decoders.front()->channel_probabilities.at(bit_index);
    }

private:
    std::vector<std::vector<int>> original_rows;
    std::vector<std::vector<int>> variable_checks;
    std::vector<std::vector<std::vector<int>>> partitions;
    std::vector<std::vector<int>> basis_source_rows;
    std::vector<std::unique_ptr<BpSparse>> basis_pcms;
    std::vector<std::unique_ptr<ldpc_modified::bp::BpDecoder>> basis_decoders;
    std::vector<std::vector<uint8_t>> converged_candidates;
    std::vector<std::vector<uint8_t>> unique_candidates;
    std::vector<int> unique_frequencies;
    std::vector<std::jthread> worker_threads;
    std::mutex decode_mutex;
    std::mutex worker_mutex;
    std::condition_variable work_available;
    std::condition_variable work_complete;
    bool stop_workers = false;
    std::size_t work_generation = 0;
    std::size_t workers_remaining = 0;
    std::atomic<std::size_t> next_basis{0};
    const std::vector<uint8_t> *active_syndrome = nullptr;
    std::vector<std::vector<uint8_t>> basis_candidates;
    std::vector<uint8_t> basis_converged;
    std::exception_ptr worker_error;

    void copy_original_graph() {
        original_rows.resize(check_count);
        variable_checks.resize(bit_count);
        for (int row = 0; row < check_count; row++) {
            for (auto &entry: pcm.iterate_row(row)) {
                original_rows[row].push_back(entry.col_index);
                variable_checks[entry.col_index].push_back(row);
            }
        }
    }

    void validate_permutations(const std::vector<std::vector<int>> &check_permutations) const {
        if (check_permutations.empty()) {
            throw std::runtime_error("At least one check permutation is required");
        }
        for (const auto &permutation: check_permutations) {
            if (permutation.size() != static_cast<std::size_t>(check_count)) {
                throw std::runtime_error("Each check permutation must contain every check");
            }
            std::vector<uint8_t> seen(check_count, 0);
            for (int check: permutation) {
                if (check < 0 || check >= check_count || seen[check]) {
                    throw std::runtime_error("Invalid check permutation");
                }
                seen[check] = 1;
            }
        }
    }

    std::vector<std::vector<int>> build_maximal_subtrees(
            const std::vector<int> &permutation) const {
        std::vector<int> check_assignment(check_count, -1);
        std::vector<std::vector<int>> subtrees;

        for (int root: permutation) {
            if (check_assignment[root] >= 0) {
                continue;
            }

            std::vector<uint8_t> visited_variables(bit_count, 0);
            std::vector<uint8_t> tested_checks(check_count, 0);
            std::vector<uint8_t> queued_checks(check_count, 0);
            std::queue<int> pending;
            std::vector<int> subtree;
            pending.push(root);
            queued_checks[root] = 1;

            while (!pending.empty()) {
                int check = pending.front();
                pending.pop();
                queued_checks[check] = 0;
                if (tested_checks[check]) {
                    continue;
                }
                tested_checks[check] = 1;

                int visited_neighbor_count = 0;
                for (int variable: original_rows[check]) {
                    visited_neighbor_count += visited_variables[variable] != 0;
                }
                if (visited_neighbor_count > 1) {
                    continue;
                }

                check_assignment[check] = static_cast<int>(subtrees.size());
                subtree.push_back(check);
                for (int variable: original_rows[check]) {
                    visited_variables[variable] = 1;
                    for (int neighbor_check: variable_checks[variable]) {
                        if (check_assignment[neighbor_check] < 0
                            && !tested_checks[neighbor_check]
                            && !queued_checks[neighbor_check]) {
                            pending.push(neighbor_check);
                            queued_checks[neighbor_check] = 1;
                        }
                    }
                }
            }
            subtrees.push_back(std::move(subtree));
        }
        return subtrees;
    }

    void build_bases(
            const std::vector<double> &channel_probabilities,
            int maximum_iterations,
            ldpc_modified::bp::BpMethod bp_method,
            ldpc_modified::bp::BpSchedule schedule,
            double min_sum_scaling_factor) {
        for (const auto &partition: partitions) {
            for (const auto &subtree: partition) {
                int duplicate_entries = 0;
                for (int row: subtree) {
                    duplicate_entries += static_cast<int>(original_rows[row].size());
                }
                auto augmented_pcm = std::make_unique<BpSparse>(
                        check_count + static_cast<int>(subtree.size()),
                        bit_count,
                        pcm.entry_count() + duplicate_entries);

                for (int row = 0; row < check_count; row++) {
                    for (int variable: original_rows[row]) {
                        augmented_pcm->insert_entry(row, variable);
                    }
                }
                for (std::size_t duplicate_index = 0; duplicate_index < subtree.size(); duplicate_index++) {
                    for (int variable: original_rows[subtree[duplicate_index]]) {
                        augmented_pcm->insert_entry(
                                check_count + static_cast<int>(duplicate_index), variable);
                    }
                }

                BpSparse *basis_pcm = augmented_pcm.get();
                basis_source_rows.push_back(subtree);
                basis_pcms.push_back(std::move(augmented_pcm));
                basis_decoders.push_back(std::make_unique<ldpc_modified::bp::BpDecoder>(
                        *basis_pcm,
                        channel_probabilities,
                        maximum_iterations,
                        bp_method,
                        schedule,
                        min_sum_scaling_factor,
                        1,
                        ldpc_modified::bp::NULL_INT_VECTOR,
                        0,
                        false,
                        ldpc_modified::bp::SYNDROME));
            }
        }
    }

    void reset_decode_state() {
        converge = false;
        std::fill(decoding.begin(), decoding.end(), 0);
        candidate_count = 0;
        unique_candidate_count = 0;
        selected_frequency = 0;
        selected_weight = 0;
        converged_candidates.clear();
        unique_candidates.clear();
        unique_frequencies.clear();
    }

    void start_workers() {
        unsigned int hardware_threads = std::thread::hardware_concurrency();
        std::size_t worker_count = std::min<std::size_t>(
                basis_decoders.size(), hardware_threads == 0 ? 1 : hardware_threads);
        worker_threads.reserve(worker_count);
        try {
            for (std::size_t worker = 0; worker < worker_count; worker++) {
                worker_threads.emplace_back([this]() { worker_loop(); });
            }
        } catch (...) {
            shutdown_workers();
            throw;
        }
    }

    void shutdown_workers() {
        {
            std::lock_guard<std::mutex> lock(worker_mutex);
            stop_workers = true;
        }
        work_available.notify_all();
        for (auto &worker: worker_threads) {
            if (worker.joinable()) {
                worker.join();
            }
        }
        worker_threads.clear();
    }

    void worker_loop() {
        std::size_t observed_generation = 0;
        while (true) {
            const std::vector<uint8_t> *syndrome;
            {
                std::unique_lock<std::mutex> lock(worker_mutex);
                work_available.wait(lock, [&]() {
                    return stop_workers || work_generation != observed_generation;
                });
                if (stop_workers) {
                    return;
                }
                observed_generation = work_generation;
                syndrome = active_syndrome;
            }

            while (true) {
                std::size_t basis_index = next_basis.fetch_add(1, std::memory_order_relaxed);
                if (basis_index >= basis_decoders.size()) {
                    break;
                }
                try {
                    std::vector<uint8_t> augmented_syndrome = *syndrome;
                    for (int source_row: basis_source_rows[basis_index]) {
                        augmented_syndrome.push_back((*syndrome)[source_row]);
                    }
                    basis_candidates[basis_index] =
                            basis_decoders[basis_index]->decode(augmented_syndrome);
                    basis_converged[basis_index] = basis_decoders[basis_index]->converge ? 1 : 0;
                } catch (...) {
                    std::lock_guard<std::mutex> lock(worker_mutex);
                    if (!worker_error) {
                        worker_error = std::current_exception();
                    }
                }
            }

            {
                std::lock_guard<std::mutex> lock(worker_mutex);
                workers_remaining -= 1;
                if (workers_remaining == 0) {
                    work_complete.notify_one();
                }
            }
        }
    }

    void run_basis_decoders(const std::vector<uint8_t> &syndrome) {
        if (worker_threads.empty()) {
            return;
        }
        basis_candidates.assign(basis_decoders.size(), {});
        basis_converged.assign(basis_decoders.size(), 0);
        {
            std::lock_guard<std::mutex> lock(worker_mutex);
            next_basis.store(0, std::memory_order_relaxed);
            active_syndrome = &syndrome;
            worker_error = nullptr;
            workers_remaining = worker_threads.size();
            work_generation += 1;
        }
        work_available.notify_all();

        std::exception_ptr decode_error;
        {
            std::unique_lock<std::mutex> lock(worker_mutex);
            work_complete.wait(lock, [&]() { return workers_remaining == 0; });
            active_syndrome = nullptr;
            decode_error = worker_error;
        }
        if (decode_error) {
            std::rethrow_exception(decode_error);
        }
    }

    static int hamming_weight(const std::vector<uint8_t> &candidate) {
        return static_cast<int>(std::count(candidate.begin(), candidate.end(), uint8_t{1}));
    }

    static int checked_size(std::size_t size) {
        if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("MBBP-LD collection exceeds the supported integer range");
        }
        return static_cast<int>(size);
    }

    void aggregate_candidates() {
        for (const auto &candidate: converged_candidates) {
            auto found = std::find(unique_candidates.begin(), unique_candidates.end(), candidate);
            if (found == unique_candidates.end()) {
                unique_candidates.push_back(candidate);
                unique_frequencies.push_back(1);
            } else {
                std::size_t index = static_cast<std::size_t>(found - unique_candidates.begin());
                unique_frequencies[index] += 1;
            }
        }

        unique_candidate_count = checked_size(unique_candidates.size());
        int best_index = 0;
        int best_weight = hamming_weight(unique_candidates[0]);
        for (int index = 1; index < unique_candidate_count; index++) {
            int weight = hamming_weight(unique_candidates[index]);
            long long candidate_score = static_cast<long long>(unique_frequencies[index])
                                        * (static_cast<long long>(best_weight) + 1);
            long long best_score = static_cast<long long>(unique_frequencies[best_index])
                                   * (static_cast<long long>(weight) + 1);
            if (candidate_score > best_score) {
                best_index = index;
                best_weight = weight;
            }
        }

        decoding = unique_candidates[best_index];
        selected_frequency = unique_frequencies[best_index];
        selected_weight = best_weight;
    }
};

}  // namespace mbbp_ld
}  // namespace ldpc_modified

#endif
