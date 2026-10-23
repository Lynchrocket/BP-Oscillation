#ifndef QUOTIENT_SPLIT_H
#define QUOTIENT_SPLIT_H

#include <algorithm>
#include <array>
#include <map>
#include <numeric>
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
#include <cassert>
#include <cstdint>
#include <string>

#include "math.h"
#include "sparse_matrix_base.hpp"
#include "gf2sparse.hpp"
#include "gf2dense.hpp"
#include "rng.hpp"

namespace ldpc_modified {
    namespace quotient_split {

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

        using SparseRows = std::vector<std::vector<int>>;

        /**
         * GF(2) quotient decomposition used by the Python orchestration layer.
         *
         * For every legal syndrome s, the model represents every compatible
         * error as
         *
         *     x = T s xor D u xor Q lambda.
         *
         * K spans ker(H), D spans ker(H) intersect ker(L), and Q contains the
         * remaining syndrome-null logical directions.  Existing gf2dense
         * elimination/kernel helpers perform the expensive linear algebra;
         * this class only assembles their sparse row representations.
         */
        class QuotientModel {
        public:
            BpSparse &h_matrix;
            BpSparse &logical_matrix;

            int check_count;
            int bit_count;
            int logical_count;
            int rank_h = 0;
            int logical_rank = 0;
            int kernel_dimension = 0;
            int degeneracy_dimension = 0;

            // Matrices are CSR-like: matrix_rows[i] lists nonzero columns in
            // row i.  Their shapes are documented beside each field.
            SparseRows T;            // n x m
            SparseRows K;            // n x dim ker(H)
            SparseRows Q;            // n x logical_rank
            SparseRows D;            // n x degeneracy_dimension
            SparseRows logical_map;  // logical_count x logical_rank, equals LQ

            QuotientModel(BpSparse &h, BpSparse &logical)
                    : h_matrix(h),
                      logical_matrix(logical),
                      check_count(h.m),
                      bit_count(h.n),
                      logical_count(logical.m) {
                if (logical.n != bit_count) {
                    throw std::invalid_argument(
                            "H and L must have the same number of columns.");
                }
                build();
            }

            ~QuotientModel() = default;

            BpSparse &d_matrix() {
                return *d_sparse;
            }

            std::vector<uint8_t> particular_solution(std::vector<uint8_t> syndrome) const {
                if (static_cast<int>(syndrome.size()) != check_count) {
                    throw std::invalid_argument(
                            "The syndrome length must equal the number of rows in H.");
                }
                auto solution = mulvec(T, check_count, syndrome);
                if (mulvec(h_rows, bit_count, solution) != syndrome) {
                    throw std::invalid_argument("The syndrome is not in im(H).");
                }
                return solution;
            }

            /**
             * Return the root offset used by quotient splitting.
             *
             * Starting from any x0 with H x0 = s, canonical logical directions
             * with L Q = I remove its logical label:
             *
             *     x_s = x0 xor Q (L x0),
             *     H x_s = s,  L x_s = 0.
             *
             * The split bits can therefore be read directly as the observable
             * label of x_s xor Q ell.
             */
            std::vector<uint8_t> logical_neutral_solution(
                    std::vector<uint8_t> syndrome) const {
                if (logical_rank != logical_count
                    || logical_map != identity(logical_count)) {
                    throw std::invalid_argument(
                            "Quotient splitting requires full logical rank and L Q = I.");
                }

                auto x0 = particular_solution(std::move(syndrome));
                const auto base_label = mulvec(logical_rows, bit_count, x0);
                const auto logical_shift = mulvec(Q, logical_count, base_label);
                auto x_s = vector_xor(x0, logical_shift);
                const std::vector<uint8_t> zero_label(logical_count, 0);
                if (mulvec(logical_rows, bit_count, x_s) != zero_label) {
                    throw std::logic_error(
                            "Failed to construct a logical-neutral split root.");
                }
                return x_s;
            }

            std::vector<uint8_t> representative(
                    std::vector<uint8_t> syndrome,
                    std::vector<uint8_t> observable_label) const {
                if (static_cast<int>(observable_label.size()) != logical_count) {
                    throw std::invalid_argument(
                            "The logical label length must equal the number of rows in L.");
                }

                auto x_s = particular_solution(syndrome);
                auto base_label = mulvec(logical_rows, bit_count, x_s);
                auto delta = vector_xor(base_label, observable_label);
                auto logical_coordinate = solve_logical_coordinate(delta);
                auto logical_shift = mulvec(Q, logical_rank, logical_coordinate);
                auto offset = vector_xor(x_s, logical_shift);

                // The representative is the offset used by D-BP.  Its full
                // coset is {offset xor D u : u in GF(2)^dim(D)}.
                if (mulvec(h_rows, bit_count, offset) != syndrome
                    || mulvec(logical_rows, bit_count, offset) != observable_label) {
                    throw std::logic_error(
                            "Internal quotient representative invariant failed.");
                }
                return offset;
            }

        private:
            SparseRows h_rows;
            SparseRows logical_rows;
            std::unique_ptr<BpSparse> d_sparse;

            static void normalize_rows(SparseRows &rows) {
                for (auto &row: rows) {
                    std::sort(row.begin(), row.end());
                    row.erase(std::unique(row.begin(), row.end()), row.end());
                }
            }

            static int count_entries(const SparseRows &rows) {
                int count = 0;
                for (const auto &row: rows) {
                    count += static_cast<int>(row.size());
                }
                return count;
            }

            static ldpc_modified::gf2dense::CscMatrix rows_to_csc(
                    const SparseRows &rows,
                    int column_count) {
                ldpc_modified::gf2dense::CscMatrix columns(
                        column_count, std::vector<int>{});
                for (int row = 0; row < static_cast<int>(rows.size()); row++) {
                    for (int column: rows[row]) {
                        if (column < 0 || column >= column_count) {
                            throw std::logic_error("Sparse matrix column is out of range.");
                        }
                        columns[column].push_back(row);
                    }
                }
                return columns;
            }

            static std::vector<uint8_t> mulvec(
                    const SparseRows &rows,
                    int column_count,
                    const std::vector<uint8_t> &input) {
                if (static_cast<int>(input.size()) != column_count) {
                    throw std::invalid_argument("GF(2) matrix-vector dimensions do not match.");
                }
                std::vector<uint8_t> output(rows.size(), 0);
                for (int row = 0; row < static_cast<int>(rows.size()); row++) {
                    for (int column: rows[row]) {
                        output[row] ^= input[column] & 1u;
                    }
                }
                return output;
            }

            static SparseRows multiply(
                    const SparseRows &left,
                    int inner_dimension,
                    const SparseRows &right,
                    int right_column_count) {
                if (static_cast<int>(right.size()) != inner_dimension) {
                    throw std::invalid_argument("GF(2) matrix dimensions do not match.");
                }
                SparseRows output(left.size(), std::vector<int>{});
                std::vector<uint8_t> parity(right_column_count, 0);
                for (int row = 0; row < static_cast<int>(left.size()); row++) {
                    std::fill(parity.begin(), parity.end(), 0);
                    for (int inner: left[row]) {
                        if (inner < 0 || inner >= inner_dimension) {
                            throw std::logic_error("Sparse matrix inner index is out of range.");
                        }
                        for (int column: right[inner]) {
                            parity[column] ^= 1;
                        }
                    }
                    for (int column = 0; column < right_column_count; column++) {
                        if (parity[column]) {
                            output[row].push_back(column);
                        }
                    }
                }
                return output;
            }

            static SparseRows select_columns(
                    const SparseRows &matrix,
                    int column_count,
                    const std::vector<int> &selected_columns) {
                std::vector<int> new_index(column_count, -1);
                for (int index = 0;
                     index < static_cast<int>(selected_columns.size());
                     index++) {
                    const int old_index = selected_columns[index];
                    if (old_index < 0 || old_index >= column_count) {
                        throw std::logic_error("Selected GF(2) column is out of range.");
                    }
                    new_index[old_index] = index;
                }

                SparseRows output(matrix.size(), std::vector<int>{});
                for (int row = 0; row < static_cast<int>(matrix.size()); row++) {
                    for (int column: matrix[row]) {
                        if (new_index[column] >= 0) {
                            output[row].push_back(new_index[column]);
                        }
                    }
                }
                return output;
            }

            static SparseRows identity(int dimension) {
                SparseRows output(dimension, std::vector<int>{});
                for (int index = 0; index < dimension; index++) {
                    output[index].push_back(index);
                }
                return output;
            }

            static std::vector<uint8_t> vector_xor(
                    const std::vector<uint8_t> &left,
                    const std::vector<uint8_t> &right) {
                if (left.size() != right.size()) {
                    throw std::invalid_argument("GF(2) vector dimensions do not match.");
                }
                std::vector<uint8_t> output(left.size(), 0);
                for (int index = 0; index < static_cast<int>(left.size()); index++) {
                    output[index] = (left[index] ^ right[index]) & 1u;
                }
                return output;
            }

            static bool is_zero(const SparseRows &matrix) {
                for (const auto &row: matrix) {
                    if (!row.empty()) {
                        return false;
                    }
                }
                return true;
            }

            static SparseRows invert_square(const SparseRows &matrix, int dimension) {
                if (static_cast<int>(matrix.size()) != dimension) {
                    throw std::invalid_argument("GF(2) inverse requires a square matrix.");
                }
                auto csc = rows_to_csc(matrix, dimension);
                ldpc_modified::gf2dense::PluDecomposition decomposition(
                        dimension, dimension, csc);
                decomposition.rref(true, true);
                if (decomposition.matrix_rank != dimension) {
                    throw std::invalid_argument("Logical basis matrix is singular over GF(2).");
                }

                SparseRows inverse(dimension, std::vector<int>{});
                for (int column = 0; column < dimension; column++) {
                    std::vector<uint8_t> unit(dimension, 0);
                    unit[column] = 1;
                    auto solution = decomposition.lu_solve(unit);
                    for (int row = 0; row < dimension; row++) {
                        if (solution[row]) {
                            inverse[row].push_back(column);
                        }
                    }
                }
                return inverse;
            }

            void build() {
                h_rows = h_matrix.row_adjacency_list();
                logical_rows = logical_matrix.row_adjacency_list();
                normalize_rows(h_rows);
                normalize_rows(logical_rows);

                auto h_for_pivots = h_rows;
                const auto independent_h_rows = ldpc_modified::gf2dense::pivot_rows(
                        check_count, bit_count, h_for_pivots);
                rank_h = static_cast<int>(independent_h_rows.size());

                // Build a generalized right inverse T on im(H).  Independent
                // rows B of H have full row rank.  Solving B x_j = e_j gives
                // columns that map the corresponding independent syndrome
                // coordinates back to an error representative, hence H T H=H.
                T.assign(bit_count, std::vector<int>{});
                if (rank_h > 0) {
                    SparseRows row_basis(rank_h, std::vector<int>{});
                    for (int row = 0; row < rank_h; row++) {
                        row_basis[row] = h_rows[independent_h_rows[row]];
                    }
                    auto row_basis_csc = rows_to_csc(row_basis, bit_count);
                    ldpc_modified::gf2dense::PluDecomposition basis_decomposition(
                            rank_h, bit_count, row_basis_csc);
                    basis_decomposition.rref(true, true);
                    if (basis_decomposition.matrix_rank != rank_h) {
                        throw std::logic_error("Independent H row basis lost rank.");
                    }
                    for (int basis_column = 0; basis_column < rank_h; basis_column++) {
                        std::vector<uint8_t> unit(rank_h, 0);
                        unit[basis_column] = 1;
                        auto solution = basis_decomposition.lu_solve(unit);
                        for (int bit = 0; bit < bit_count; bit++) {
                            if (solution[bit]) {
                                T[bit].push_back(independent_h_rows[basis_column]);
                            }
                        }
                    }
                }

                auto h_for_kernel = h_rows;
                K = ldpc_modified::gf2dense::kernel(
                        check_count, bit_count, h_for_kernel);
                kernel_dimension = bit_count - rank_h;
                if (static_cast<int>(K.size()) != bit_count) {
                    throw std::logic_error("gf2dense::kernel returned an invalid K shape.");
                }

                // M maps a syndrome-null direction K a to its observable:
                // L(K a) = (L K) a = M a.
                auto M = multiply(
                        logical_rows, bit_count, K, kernel_dimension);
                auto m_csc = rows_to_csc(M, kernel_dimension);
                ldpc_modified::gf2dense::PluDecomposition m_decomposition(
                        logical_count, kernel_dimension, m_csc);
                m_decomposition.rref(false, false);
                logical_rank = m_decomposition.matrix_rank;

                auto m_for_kernel = M;
                auto kernel_m = ldpc_modified::gf2dense::kernel(
                        logical_count, kernel_dimension, m_for_kernel);
                degeneracy_dimension = kernel_dimension - logical_rank;
                D = multiply(
                        K, kernel_dimension, kernel_m, degeneracy_dimension);

                const auto logical_pivots = m_decomposition.pivot_cols;
                auto Q0 = select_columns(K, kernel_dimension, logical_pivots);
                auto R = multiply(
                        logical_rows, bit_count, Q0, logical_rank);
                if (logical_rank == logical_count) {
                    auto inverse_R = invert_square(R, logical_count);
                    Q = multiply(Q0, logical_rank, inverse_R, logical_count);
                    logical_map = identity(logical_count);
                } else {
                    Q = std::move(Q0);
                    logical_map = std::move(R);
                }

                validate_invariants();

                d_sparse = std::make_unique<BpSparse>(
                        bit_count, degeneracy_dimension, count_entries(D));
                for (int row = 0; row < bit_count; row++) {
                    for (int column: D[row]) {
                        d_sparse->insert_entry(row, column);
                    }
                }
            }

            void validate_invariants() const {
                // T is a generalized right inverse on im(H): H T H = H.
                const auto HT = multiply(h_rows, bit_count, T, check_count);
                const auto HTH = multiply(HT, check_count, h_rows, bit_count);
                const auto HK = multiply(h_rows, bit_count, K, kernel_dimension);
                const auto HD = multiply(h_rows, bit_count, D, degeneracy_dimension);
                const auto LD = multiply(
                        logical_rows, bit_count, D, degeneracy_dimension);
                const auto HQ = multiply(h_rows, bit_count, Q, logical_rank);
                const auto LQ = multiply(
                        logical_rows, bit_count, Q, logical_rank);

                if (HTH != h_rows || !is_zero(HK) || !is_zero(HD)
                    || !is_zero(LD) || !is_zero(HQ) || LQ != logical_map
                    || degeneracy_dimension + logical_rank != kernel_dimension) {
                    throw std::logic_error("T/K/Q/D quotient decomposition invariant failed.");
                }
            }

            std::vector<uint8_t> solve_logical_coordinate(
                    const std::vector<uint8_t> &delta) const {
                if (static_cast<int>(delta.size()) != logical_count) {
                    throw std::invalid_argument("Logical delta has an invalid length.");
                }
                if (logical_rank == logical_count) {
                    // The canonical basis construction above guarantees LQ=I.
                    return delta;
                }
                if (logical_rank == 0) {
                    if (std::any_of(delta.begin(), delta.end(), [](uint8_t bit) {
                        return bit != 0;
                    })) {
                        throw std::invalid_argument("The logical label is unreachable.");
                    }
                    return {};
                }

                auto map_csc = rows_to_csc(logical_map, logical_rank);
                ldpc_modified::gf2dense::PluDecomposition decomposition(
                        logical_count, logical_rank, map_csc);
                decomposition.rref(true, true);
                auto right_hand_side = delta;
                auto coordinate = decomposition.lu_solve(right_hand_side);
                if (mulvec(logical_map, logical_rank, coordinate) != delta) {
                    throw std::invalid_argument("The logical label is unreachable.");
                }
                return coordinate;
            }
        };

        /**
         * Sum-product BP on the quotient factor graph D.
         *
         * D has one factor per physical error mechanism and one variable per
         * degeneracy direction.  For offset x, factor i is
         *
         *   psi_i(u_N(i)) = P_i(x_i xor XOR_(a in N(i)) u_a).
         *
         * This is not ordinary syndrome BP, so it has a dedicated update and
         * returns a Bethe log-partition rather than a binary decoding.
         */
        class DPartitionBpDecoder {
        public:
            BpSparse &d_matrix;
            std::vector<double> channel_probabilities;
            int maximum_iterations;
            double damping;
            double tolerance;
            int iterations = 0;
            bool converge = false;
            double max_message_change = std::numeric_limits<double>::infinity();
            double log_partition = -std::numeric_limits<double>::infinity();

            DPartitionBpDecoder(
                    BpSparse &d,
                    std::vector<double> probabilities,
                    int maximum_iterations = 100,
                    double damping = 0.2,
                    double tolerance = 1e-10)
                    : d_matrix(d),
                      channel_probabilities(std::move(probabilities)),
                      maximum_iterations(maximum_iterations),
                      damping(damping),
                      tolerance(tolerance) {
                validate_parameters();
            }

            ~DPartitionBpDecoder() = default;

            void set_channel_probabilities(
                    const std::vector<double> &probabilities) {
                channel_probabilities = probabilities;
                validate_probabilities();
            }

            double decode(std::vector<uint8_t> &offset) {
                if (static_cast<int>(offset.size()) != d_matrix.m) {
                    throw std::invalid_argument(
                            "The D-BP offset must have one entry per row of D.");
                }
                validate_probabilities();
                initialise_messages();
                converge = false;
                iterations = 0;
                max_message_change = std::numeric_limits<double>::infinity();

                constexpr double epsilon = 1e-15;
                for (int iteration = 1; iteration <= maximum_iterations; iteration++) {
                    double iteration_change = 0.0;

                    for (int factor = 0; factor < d_matrix.m; factor++) {
                        const double coupling =
                                (1.0 - 2.0 * channel_probabilities[factor])
                                * (offset[factor] ? -1.0 : 1.0);
                        for (auto &edge: d_matrix.iterate_row(factor)) {
                            double product_other = 1.0;
                            for (auto &other: d_matrix.iterate_row(factor)) {
                                if (&other != &edge) {
                                    product_other *= std::tanh(
                                            other.bit_to_check_msg / 2.0);
                                }
                            }

                            // With c_i=(1-2p_i)(-1)^x_i, normalized factor
                            // message is
                            //   m_(i->a)(u) = [1 + (-1)^u q]/2,
                            //   q = c_i product_(b!=a) tanh(m_(b->i)/2).
                            const double q = std::clamp(
                                    coupling * product_other,
                                    -1.0 + epsilon,
                                    1.0 - epsilon);
                            const double proposed_probability_zero = (1.0 + q) / 2.0;
                            const double old_probability_zero = probability_zero(
                                    edge.check_to_bit_msg);
                            const double next_probability_zero =
                                    (1.0 - damping) * proposed_probability_zero
                                    + damping * old_probability_zero;
                            iteration_change = std::max(
                                    iteration_change,
                                    std::abs(next_probability_zero
                                             - old_probability_zero));
                            edge.check_to_bit_msg = llr_from_probability_zero(
                                    next_probability_zero);
                        }
                    }

                    for (int variable = 0; variable < d_matrix.n; variable++) {
                        for (auto &edge: d_matrix.iterate_column(variable)) {
                            // Degeneracy variables have a uniform prior, hence
                            // m_(a->i) = sum_(j in N(a)\\i) m_(j->a) in LLR form.
                            double proposed_llr = 0.0;
                            for (auto &other: d_matrix.iterate_column(variable)) {
                                if (&other != &edge) {
                                    proposed_llr += other.check_to_bit_msg;
                                }
                            }
                            const double proposed_probability_zero = probability_zero(
                                    proposed_llr);
                            const double old_probability_zero = probability_zero(
                                    edge.bit_to_check_msg);
                            const double next_probability_zero =
                                    (1.0 - damping) * proposed_probability_zero
                                    + damping * old_probability_zero;
                            iteration_change = std::max(
                                    iteration_change,
                                    std::abs(next_probability_zero
                                             - old_probability_zero));
                            edge.bit_to_check_msg = llr_from_probability_zero(
                                    next_probability_zero);
                        }
                    }

                    iterations = iteration;
                    max_message_change = iteration_change;
                    if (iteration_change < tolerance) {
                        converge = true;
                        break;
                    }
                }

                log_partition = compute_bethe_log_partition(offset);
                return log_partition;
            }

        private:
            static double probability_zero(double llr) {
                if (llr >= 0.0) {
                    return 1.0 / (1.0 + std::exp(-llr));
                }
                const double exponential = std::exp(llr);
                return exponential / (1.0 + exponential);
            }

            static double llr_from_probability_zero(double probability) {
                constexpr double epsilon = 1e-15;
                const double clipped = std::clamp(
                        probability, epsilon, 1.0 - epsilon);
                return std::log(clipped / (1.0 - clipped));
            }

            static double log_sum_exp(double left, double right) {
                const double maximum = std::max(left, right);
                return maximum + std::log(
                        std::exp(left - maximum) + std::exp(right - maximum));
            }

            void validate_parameters() const {
                if (maximum_iterations < 1) {
                    throw std::invalid_argument("D-BP maximum_iterations must be positive.");
                }
                if (!(damping >= 0.0 && damping < 1.0)) {
                    throw std::invalid_argument("D-BP damping must lie in [0, 1).");
                }
                if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
                    throw std::invalid_argument("D-BP tolerance must be finite and positive.");
                }
                validate_probabilities();
            }

            void validate_probabilities() const {
                if (static_cast<int>(channel_probabilities.size()) != d_matrix.m) {
                    throw std::invalid_argument(
                            "D-BP needs one channel probability per row of D.");
                }
                for (double probability: channel_probabilities) {
                    if (!(probability > 0.0 && probability < 1.0)
                        || !std::isfinite(probability)) {
                        throw std::invalid_argument(
                                "D-BP probabilities must be finite and lie in (0, 1).");
                    }
                }
            }

            void initialise_messages() {
                // u variables have uniform priors, so both directions start
                // with LLR log(P(u=0)/P(u=1)) = 0.
                for (int factor = 0; factor < d_matrix.m; factor++) {
                    for (auto &edge: d_matrix.iterate_row(factor)) {
                        edge.bit_to_check_msg = 0.0;
                        edge.check_to_bit_msg = 0.0;
                    }
                }
            }

            double compute_bethe_log_partition(
                    const std::vector<uint8_t> &offset) {
                // Bethe approximation:
                //   log Z_B = sum_i log Z_i + sum_a log Z_a
                //             - sum_(i,a) log Z_(i,a).
                double result = 0.0;

                for (int factor = 0; factor < d_matrix.m; factor++) {
                    double factor_normalization;
                    if (d_matrix.get_row_degree(factor) == 0) {
                        factor_normalization = offset[factor]
                                               ? channel_probabilities[factor]
                                               : 1.0 - channel_probabilities[factor];
                    } else {
                        const double coupling =
                                (1.0 - 2.0 * channel_probabilities[factor])
                                * (offset[factor] ? -1.0 : 1.0);
                        double product_magnetization = 1.0;
                        for (auto &edge: d_matrix.iterate_row(factor)) {
                            product_magnetization *= std::tanh(
                                    edge.bit_to_check_msg / 2.0);
                        }
                        factor_normalization = 0.5
                                               * (1.0 + coupling
                                                        * product_magnetization);
                    }
                    if (!(factor_normalization > 0.0)) {
                        throw std::runtime_error(
                                "D-BP factor normalization is non-positive.");
                    }
                    result += std::log(factor_normalization);
                }

                for (int variable = 0; variable < d_matrix.n; variable++) {
                    double log_zero = 0.0;
                    double log_one = 0.0;
                    for (auto &edge: d_matrix.iterate_column(variable)) {
                        const double probability = probability_zero(
                                edge.check_to_bit_msg);
                        log_zero += std::log(probability);
                        log_one += std::log(1.0 - probability);
                    }
                    result += log_sum_exp(log_zero, log_one);
                }

                for (int factor = 0; factor < d_matrix.m; factor++) {
                    for (auto &edge: d_matrix.iterate_row(factor)) {
                        const double variable_zero = probability_zero(
                                edge.bit_to_check_msg);
                        const double factor_zero = probability_zero(
                                edge.check_to_bit_msg);
                        const double edge_normalization =
                                variable_zero * factor_zero
                                + (1.0 - variable_zero) * (1.0 - factor_zero);
                        if (!(edge_normalization > 0.0)) {
                            throw std::runtime_error(
                                    "D-BP edge normalization is non-positive.");
                        }
                        result -= std::log(edge_normalization);
                    }
                }

                if (!std::isfinite(result)) {
                    throw std::runtime_error("D-BP Bethe log partition is non-finite.");
                }
                return result;
            }
        };

        class BpDecoder {
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

            // Optional quotient-BP trajectory tracking.  The ordinary BP path
            // leaves logical_matrix unset, so its behavior is unchanged.
            BpSparse *logical_matrix = nullptr;
            bool collect_full_trajectory = false;
            std::vector<std::vector<uint8_t>> trajectory_logical_cosets;
            std::vector<std::vector<uint8_t>> logical_cosets;
            std::set<std::vector<uint8_t>> logical_coset_set;

            BpDecoder(
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


                if (this->channel_probabilities.size() != this->bit_count) {
                    throw std::runtime_error(
                            "Channel probabilities vector must have length equal to the number of bits");
                }
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

            ~BpDecoder() = default;

            void set_omp_thread_count(int count) {
                this->omp_thread_count = count;
                // omp_set_num_threads(this->omp_thread_count);
                // NotImplemented
            }

            void set_random_schedule_seed(int seed) {
                this->random_schedule_seed = seed;
                this->rng_list_shuffle.seed(seed);
            }

            void enable_logical_tracking(BpSparse &logical, bool collect_full = true) {
                if (logical.n != this->bit_count) {
                    throw std::invalid_argument(
                            "The logical matrix and H must have the same number of columns.");
                }
                this->logical_matrix = &logical;
                this->collect_full_trajectory = collect_full;
            }

            void reset_logical_tracking() {
                this->trajectory_logical_cosets.clear();
                this->logical_cosets.clear();
                this->logical_coset_set.clear();
            }

            void record_logical_coset() {
                if (this->logical_matrix == nullptr) {
                    return;
                }

                // At H-BP iteration t, belief b_i^(t) produces
                // x_hat_i^(t) = 1[b_i^(t) <= 0].  Its logical coset label is
                //
                //     ell^(t) = L x_hat^(t)  (mod 2).
                //
                // The set removes repeated visits while logical_cosets keeps
                // deterministic first-appearance order for D-BP reranking.
                auto label = this->logical_matrix->mulvec(this->decoding);
                this->trajectory_logical_cosets.push_back(label);
                if (this->logical_coset_set.insert(label).second) {
                    this->logical_cosets.push_back(std::move(label));
                }
            }

            void initialise_log_domain_bp() {
                // initialise BP
                for (int i = 0; i < this->bit_count; i++) {
                    this->initial_log_prob_ratios[i] = std::log(
                            (1 - this->channel_probabilities[i]) / this->channel_probabilities[i]);

                    for (auto &e: this->pcm.iterate_column(i)) {
                        e.bit_to_check_msg = this->initial_log_prob_ratios[i];
                    }
                }
            }

            std::vector<uint8_t> decode(std::vector<uint8_t> &input_vector) {

                this->reset_logical_tracking();


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
                        // if(isnan(log_prob_ratios[i])) log_prob_ratios[i] = initial_log_prob_ratios[i];
                        if (temp <= 0) {
                            this->decoding[i] = 1;
                            for (auto &e: this->pcm.iterate_column(i)) {
                                this->candidate_syndrome[e.row_index] ^= 1;
                            }
                        } else {
                            this->decoding[i] = 0;
                        }
                    }

                    this->record_logical_coset();

                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        this->converge = true;
                    }

                    this->iterations = it;

                    if (this->converge && !this->collect_full_trajectory) {
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
                                      const double reliability1 = it != 1
                                                                  ? this->log_prob_ratios[bit1]
                                                                  : std::log((1 - channel_probabilities[bit1]) /
                                                                             channel_probabilities[bit1]);
                                      const double reliability2 = it != 1
                                                                  ? this->log_prob_ratios[bit2]
                                                                  : std::log((1 - channel_probabilities[bit2]) /
                                                                             channel_probabilities[bit2]);
                                      if (reliability1 == reliability2) {
                                          // A fixed tie-break prevents the previous decode's
                                          // in-place order from changing the next trajectory.
                                          return bit1 < bit2;
                                      }
                                      return reliability1 > reliability2;
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
                        if (this->log_prob_ratios[bit_index] <= 0) {
                            this->decoding[bit_index] = 1;
                        } else {
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
                    this->record_logical_coset();
                    this->iterations = it;
                    if (std::equal(candidate_syndrome.begin(), candidate_syndrome.end(), syndrome.begin())) {
                        this->converge = true;
                        if (!this->collect_full_trajectory) {
                            return this->decoding;
                        }
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

        /** A pair of normalized binary messages stored for one factor edge. */
        struct GeneratorMessagePair {
            std::array<double, 2> variable_to_factor{0.5, 0.5};
            std::array<double, 2> factor_to_variable{0.5, 0.5};
        };

        using GeneratorMessageState =
                std::map<std::pair<int, int>, GeneratorMessagePair>;

        /** Result of sum-product BP on an affine generator graph. */
        struct GeneratorBpResult {
            bool converged = false;
            int iterations = 0;
            std::vector<double> variable_probability_one;
            double log_partition = -std::numeric_limits<double>::infinity();
            double max_message_change = std::numeric_limits<double>::infinity();
            GeneratorMessageState message_state;
        };

        /**
         * Sum-product BP for an affine family x = offset xor G w.
         *
         * Each physical error mechanism i is a factor with potential
         *
         *   psi_i(w_N(i)) = P_i(offset_i xor XOR_(a in N(i)) w_a).
         *
         * coordinate_ids give every column a stable global identity.  Message
         * states can therefore warm-start a child graph after one Q column has
         * been removed: messages on all surviving (factor, coordinate) edges
         * retain their previous values.
         */
        class GeneratorGraphBp {
        public:
            GeneratorGraphBp(
                    SparseRows generator_rows,
                    int variable_count,
                    std::vector<uint8_t> offset,
                    std::vector<double> probabilities,
                    std::vector<int> coordinate_ids)
                    : generator_rows(std::move(generator_rows)),
                      variable_count(variable_count),
                      offset(std::move(offset)),
                      probabilities(std::move(probabilities)),
                      coordinate_ids(std::move(coordinate_ids)) {
                validate_input();
                build_edges();
            }

            GeneratorBpResult run(
                    int maximum_iterations,
                    double damping,
                    double tolerance,
                    const GeneratorMessageState *initial_state = nullptr) const {
                if (maximum_iterations < 0) {
                    throw std::invalid_argument(
                            "Generator BP iterations must be non-negative.");
                }
                if (!(damping >= 0.0 && damping < 1.0)) {
                    throw std::invalid_argument(
                            "Generator BP damping must lie in [0, 1).");
                }
                if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
                    throw std::invalid_argument(
                            "Generator BP tolerance must be finite and positive.");
                }

                std::vector<BinaryMessage> variable_to_factor(
                        edges.size(), BinaryMessage{0.5, 0.5});
                std::vector<BinaryMessage> factor_to_variable(
                        edges.size(), BinaryMessage{0.5, 0.5});
                if (initial_state != nullptr) {
                    for (int edge = 0; edge < static_cast<int>(edges.size()); edge++) {
                        const auto stored = initial_state->find(edge_key(edge));
                        if (stored != initial_state->end()) {
                            variable_to_factor[edge] = normalize(
                                    stored->second.variable_to_factor);
                            factor_to_variable[edge] = normalize(
                                    stored->second.factor_to_variable);
                        }
                    }
                }

                GeneratorBpResult result;
                for (int iteration = 1; iteration <= maximum_iterations; iteration++) {
                    const auto proposed_factor = factor_update(variable_to_factor);
                    auto next_factor = damp_messages(
                            factor_to_variable, proposed_factor, damping);
                    const auto proposed_variable = variable_update(next_factor);
                    auto next_variable = damp_messages(
                            variable_to_factor, proposed_variable, damping);

                    result.max_message_change = std::max(
                            max_difference(factor_to_variable, next_factor),
                            max_difference(variable_to_factor, next_variable));
                    factor_to_variable = std::move(next_factor);
                    variable_to_factor = std::move(next_variable);
                    result.iterations = iteration;
                    if (result.max_message_change < tolerance) {
                        result.converged = true;
                        break;
                    }
                }

                result.variable_probability_one = variable_beliefs(
                        factor_to_variable);
                result.log_partition = bethe_log_partition(
                        variable_to_factor, factor_to_variable);
                for (int edge = 0; edge < static_cast<int>(edges.size()); edge++) {
                    result.message_state.emplace(
                            edge_key(edge),
                            GeneratorMessagePair{
                                    variable_to_factor[edge],
                                    factor_to_variable[edge]});
                }
                return result;
            }

        private:
            using BinaryMessage = std::array<double, 2>;

            struct Edge {
                int factor;
                int local_variable;
            };

            SparseRows generator_rows;
            int variable_count;
            std::vector<uint8_t> offset;
            std::vector<double> probabilities;
            std::vector<int> coordinate_ids;
            std::vector<Edge> edges;
            std::vector<std::vector<int>> factor_edges;
            std::vector<std::vector<int>> variable_edges;

            static BinaryMessage normalize(BinaryMessage values) {
                constexpr double epsilon = 1e-15;
                values[0] = std::max(values[0], epsilon);
                values[1] = std::max(values[1], epsilon);
                const double total = values[0] + values[1];
                if (!(total > 0.0) || !std::isfinite(total)) {
                    throw std::runtime_error(
                            "Generator BP produced an invalid binary message.");
                }
                values[0] /= total;
                values[1] /= total;
                return values;
            }

            void validate_input() const {
                const int factor_count = static_cast<int>(generator_rows.size());
                if (static_cast<int>(offset.size()) != factor_count
                    || static_cast<int>(probabilities.size()) != factor_count) {
                    throw std::invalid_argument(
                            "Generator, offset, and probabilities have incompatible shapes.");
                }
                if (variable_count < 0
                    || static_cast<int>(coordinate_ids.size()) != variable_count) {
                    throw std::invalid_argument(
                            "Generator coordinate ids do not match its column count.");
                }
                auto sorted_ids = coordinate_ids;
                std::sort(sorted_ids.begin(), sorted_ids.end());
                if (std::adjacent_find(sorted_ids.begin(), sorted_ids.end())
                    != sorted_ids.end()) {
                    throw std::invalid_argument(
                            "Generator coordinate ids must be unique.");
                }
                for (double probability: probabilities) {
                    if (!std::isfinite(probability)
                        || !(probability > 0.0 && probability < 1.0)) {
                        throw std::invalid_argument(
                                "Every channel probability must lie strictly in (0, 1).");
                    }
                }
                for (const auto &row: generator_rows) {
                    for (int variable: row) {
                        if (variable < 0 || variable >= variable_count) {
                            throw std::invalid_argument(
                                    "Generator entry is outside its declared shape.");
                        }
                    }
                }
            }

            void build_edges() {
                factor_edges.assign(generator_rows.size(), std::vector<int>{});
                variable_edges.assign(variable_count, std::vector<int>{});
                for (int factor = 0;
                     factor < static_cast<int>(generator_rows.size());
                     factor++) {
                    for (int local_variable: generator_rows[factor]) {
                        const int edge = static_cast<int>(edges.size());
                        edges.push_back(Edge{factor, local_variable});
                        factor_edges[factor].push_back(edge);
                        variable_edges[local_variable].push_back(edge);
                    }
                }
            }

            std::pair<int, int> edge_key(int edge) const {
                return {
                        edges[edge].factor,
                        coordinate_ids[edges[edge].local_variable]};
            }

            std::vector<BinaryMessage> factor_update(
                    const std::vector<BinaryMessage> &variable_to_factor) const {
                std::vector<BinaryMessage> output(
                        edges.size(), BinaryMessage{0.5, 0.5});
                for (int factor = 0;
                     factor < static_cast<int>(factor_edges.size());
                     factor++) {
                    const auto &neighbors = factor_edges[factor];
                    const int degree = static_cast<int>(neighbors.size());
                    std::vector<double> biases(degree, 0.0);
                    std::vector<double> prefix(degree + 1, 1.0);
                    std::vector<double> suffix(degree + 1, 1.0);
                    for (int local = 0; local < degree; local++) {
                        const int edge = neighbors[local];
                        biases[local] = variable_to_factor[edge][0]
                                - variable_to_factor[edge][1];
                        prefix[local + 1] = prefix[local] * biases[local];
                    }
                    for (int local = degree - 1; local >= 0; local--) {
                        suffix[local] = suffix[local + 1] * biases[local];
                    }

                    const BinaryMessage likelihood{
                            1.0 - probabilities[factor], probabilities[factor]};
                    for (int local = 0; local < degree; local++) {
                        const int target_edge = neighbors[local];
                        // E[(-1)^parity_except_target] is the product of the
                        // incoming biases. Prefix/suffix products evaluate it
                        // for every target edge in O(degree), rather than
                        // recomputing each product in O(degree^2).
                        const double parity_bias =
                                prefix[local] * suffix[local + 1];
                        const double probability_even = 0.5 * (1.0 + parity_bias);
                        const double probability_odd = 0.5 * (1.0 - parity_bias);
                        BinaryMessage message{};
                        for (int bit = 0; bit < 2; bit++) {
                            const int even_error = offset[factor] ^ bit;
                            message[bit] =
                                    probability_even * likelihood[even_error]
                                    + probability_odd * likelihood[even_error ^ 1];
                        }
                        output[target_edge] = normalize(message);
                    }
                }
                return output;
            }

            std::vector<BinaryMessage> variable_update(
                    const std::vector<BinaryMessage> &factor_to_variable) const {
                constexpr double epsilon = 1e-15;
                std::vector<BinaryMessage> output(
                        edges.size(), BinaryMessage{0.5, 0.5});
                for (const auto &neighbors: variable_edges) {
                    BinaryMessage total_log_message{0.0, 0.0};
                    for (int edge: neighbors) {
                        total_log_message[0] += std::log(std::max(
                                factor_to_variable[edge][0], epsilon));
                        total_log_message[1] += std::log(std::max(
                                factor_to_variable[edge][1], epsilon));
                    }
                    for (int target_edge: neighbors) {
                        // m_{v->f}(b) is proportional to the product of all
                        // incoming m_{g->v}(b), g != f.  Subtracting the target
                        // log-message from the precomputed total evaluates all
                        // outgoing edges of this variable in O(degree).
                        BinaryMessage log_message{
                                total_log_message[0] - std::log(std::max(
                                        factor_to_variable[target_edge][0], epsilon)),
                                total_log_message[1] - std::log(std::max(
                                        factor_to_variable[target_edge][1], epsilon))};
                        const double maximum = std::max(
                                log_message[0], log_message[1]);
                        output[target_edge] = normalize(BinaryMessage{
                                std::exp(log_message[0] - maximum),
                                std::exp(log_message[1] - maximum)});
                    }
                }
                return output;
            }

            static std::vector<BinaryMessage> damp_messages(
                    const std::vector<BinaryMessage> &old_messages,
                    const std::vector<BinaryMessage> &proposed_messages,
                    double damping) {
                std::vector<BinaryMessage> output(old_messages.size());
                for (int edge = 0; edge < static_cast<int>(output.size()); edge++) {
                    output[edge] = normalize(BinaryMessage{
                            damping * old_messages[edge][0]
                            + (1.0 - damping) * proposed_messages[edge][0],
                            damping * old_messages[edge][1]
                            + (1.0 - damping) * proposed_messages[edge][1]});
                }
                return output;
            }

            static double max_difference(
                    const std::vector<BinaryMessage> &left,
                    const std::vector<BinaryMessage> &right) {
                double result = 0.0;
                for (int edge = 0; edge < static_cast<int>(left.size()); edge++) {
                    result = std::max(result, std::abs(left[edge][0] - right[edge][0]));
                    result = std::max(result, std::abs(left[edge][1] - right[edge][1]));
                }
                return result;
            }

            std::vector<double> variable_beliefs(
                    const std::vector<BinaryMessage> &factor_to_variable) const {
                std::vector<double> output(variable_count, 0.5);
                for (int variable = 0; variable < variable_count; variable++) {
                    BinaryMessage log_belief{0.0, 0.0};
                    for (int edge: variable_edges[variable]) {
                        log_belief[0] += std::log(factor_to_variable[edge][0]);
                        log_belief[1] += std::log(factor_to_variable[edge][1]);
                    }
                    const double maximum = std::max(log_belief[0], log_belief[1]);
                    const auto belief = normalize(BinaryMessage{
                            std::exp(log_belief[0] - maximum),
                            std::exp(log_belief[1] - maximum)});
                    output[variable] = belief[1];
                }
                return output;
            }

            double bethe_log_partition(
                    const std::vector<BinaryMessage> &variable_to_factor,
                    const std::vector<BinaryMessage> &factor_to_variable) const {
                constexpr double epsilon = 1e-15;
                double result = 0.0;
                for (int factor = 0;
                     factor < static_cast<int>(factor_edges.size());
                     factor++) {
                    double parity_bias = 1.0;
                    for (int edge: factor_edges[factor]) {
                        parity_bias *=
                                variable_to_factor[edge][0]
                                - variable_to_factor[edge][1];
                    }
                    const double probability_even = 0.5 * (1.0 + parity_bias);
                    const double probability_odd = 0.5 * (1.0 - parity_bias);
                    const double factor_normalization =
                            probability_even
                            * (offset[factor]
                               ? probabilities[factor]
                               : 1.0 - probabilities[factor])
                            + probability_odd
                            * (offset[factor]
                               ? 1.0 - probabilities[factor]
                               : probabilities[factor]);
                    result += std::log(std::max(factor_normalization, epsilon));
                }

                for (const auto &neighbors: variable_edges) {
                    double product_zero = 1.0;
                    double product_one = 1.0;
                    for (int edge: neighbors) {
                        product_zero *= factor_to_variable[edge][0];
                        product_one *= factor_to_variable[edge][1];
                    }
                    result += std::log(std::max(
                            product_zero + product_one, epsilon));
                }

                for (int edge = 0; edge < static_cast<int>(edges.size()); edge++) {
                    const double edge_normalization =
                            variable_to_factor[edge][0] * factor_to_variable[edge][0]
                            + variable_to_factor[edge][1] * factor_to_variable[edge][1];
                    result -= std::log(std::max(edge_normalization, epsilon));
                }
                if (!std::isfinite(result)) {
                    throw std::runtime_error(
                            "Generator BP Bethe log partition is non-finite.");
                }
                return result;
            }
        };

        struct QuotientSplitState {
            std::vector<uint8_t> logical_bits;
            std::vector<int> remaining;
            std::vector<uint8_t> offset;
            double path_score = 0.0;
            std::shared_ptr<GeneratorBpResult> bp_result;
        };

        struct QuotientSplitLeafResult {
            std::vector<uint8_t> logical_label;
            std::vector<uint8_t> correction;
            double path_score = 0.0;
            double log_partition = -std::numeric_limits<double>::infinity();
            bool bp_converged = false;
            int bp_iterations = 0;
        };

        struct QuotientSplitDecodeResult {
            std::vector<uint8_t> best_logical_label;
            std::vector<uint8_t> correction;
            std::vector<QuotientSplitLeafResult> leaves;
            std::vector<int> surviving_prefix_counts;
        };

        /**
         * Quotient successive-cancellation list BP.
         *
         * Internal nodes sum over u and all remaining logical bits using the
         * generator [D, Q_remaining].  Splitting coordinate j fixes ell_j to b
         * and updates the affine offset by b Q_j.  The beam score is
         *
         *   S_child = S_parent + log P(ell_j=b | current affine block).
         *
         * Once every logical bit is fixed, leaves are ranked only by a fresh
         * D-only Bethe partition estimate, exactly as in the Python reference.
         */
        class QuotientSplitDecoder {
        public:
            QuotientModel &model;
            std::vector<double> channel_probabilities;
            int list_width;
            std::vector<int> split_order;
            bool adaptive_split;
            int internal_iterations;
            int leaf_iterations;
            double damping;
            double tolerance;
            double probability_floor;

            QuotientSplitDecoder(
                    QuotientModel &model,
                    std::vector<double> probabilities,
                    int list_width = 2,
                    std::vector<int> split_order = {},
                    bool adaptive_split = false,
                    int internal_iterations = 5,
                    int leaf_iterations = 50,
                    double damping = 0.2,
                    double tolerance = 1e-9,
                    double probability_floor = 1e-300)
                    : model(model),
                      channel_probabilities(std::move(probabilities)),
                      list_width(list_width),
                      split_order(std::move(split_order)),
                      adaptive_split(adaptive_split),
                      internal_iterations(internal_iterations),
                      leaf_iterations(leaf_iterations),
                      damping(damping),
                      tolerance(tolerance),
                      probability_floor(probability_floor) {
                if (this->split_order.empty() && model.logical_count > 0) {
                    this->split_order.resize(model.logical_count);
                    std::iota(this->split_order.begin(), this->split_order.end(), 0);
                }
                validate_parameters();
            }

            void set_channel_probabilities(
                    const std::vector<double> &probabilities) {
                validate_probabilities(probabilities);
                channel_probabilities = probabilities;
            }

            QuotientSplitDecodeResult decode(std::vector<uint8_t> syndrome) const {
                QuotientSplitState root;
                root.logical_bits.assign(model.logical_count, unfixed_value);
                root.remaining.resize(model.logical_count);
                std::iota(root.remaining.begin(), root.remaining.end(), 0);
                root.offset = model.logical_neutral_solution(std::move(syndrome));

                std::vector<QuotientSplitState> active{
                        run_internal_bp(std::move(root))};
                std::vector<int> surviving_counts{1};

                for (int depth = 0; depth < model.logical_count; depth++) {
                    const int logical_index = choose_split_coordinate(active, depth);
                    std::vector<QuotientSplitState> candidates;
                    candidates.reserve(active.size() * 2);
                    for (const auto &state: active) {
                        auto children = branch(state, logical_index);
                        candidates.push_back(std::move(children[0]));
                        candidates.push_back(std::move(children[1]));
                    }

                    // stable_sort preserves bit=0 before bit=1 and the previous
                    // beam order when path scores tie.
                    std::stable_sort(
                            candidates.begin(),
                            candidates.end(),
                            [](const auto &left, const auto &right) {
                                return left.path_score > right.path_score;
                            });
                    if (static_cast<int>(candidates.size()) > list_width) {
                        candidates.resize(list_width);
                    }
                    active = std::move(candidates);
                    if (depth < model.logical_count - 1) {
                        for (auto &state: active) {
                            state = run_internal_bp(std::move(state));
                        }
                    }
                    surviving_counts.push_back(static_cast<int>(active.size()));
                }

                std::vector<QuotientSplitLeafResult> leaves;
                leaves.reserve(active.size());
                for (const auto &state: active) {
                    leaves.push_back(evaluate_leaf(state));
                }
                std::stable_sort(
                        leaves.begin(),
                        leaves.end(),
                        [](const auto &left, const auto &right) {
                            return left.log_partition > right.log_partition;
                        });
                if (leaves.empty()) {
                    throw std::logic_error("Quotient split produced no leaves.");
                }

                QuotientSplitDecodeResult result;
                result.best_logical_label = leaves[0].logical_label;
                result.correction = leaves[0].correction;
                result.leaves = std::move(leaves);
                result.surviving_prefix_counts = std::move(surviving_counts);
                return result;
            }

        private:
            static constexpr uint8_t unfixed_value = 255;

            void validate_probabilities(
                    const std::vector<double> &probabilities) const {
                if (static_cast<int>(probabilities.size()) != model.bit_count) {
                    throw std::invalid_argument(
                            "The error channel must have one entry per physical bit.");
                }
                for (double probability: probabilities) {
                    if (!std::isfinite(probability)
                        || !(probability > 0.0 && probability < 1.0)) {
                        throw std::invalid_argument(
                                "Every error-channel entry must lie strictly in (0, 1).");
                    }
                }
            }

            void validate_parameters() const {
                validate_probabilities(channel_probabilities);
                if (model.logical_rank != model.logical_count) {
                    throw std::invalid_argument(
                            "Quotient splitting requires full logical rank.");
                }
                if (list_width < 1) {
                    throw std::invalid_argument("list_width must be at least 1.");
                }
                if (internal_iterations < 0 || leaf_iterations < 0) {
                    throw std::invalid_argument(
                            "BP iteration counts must be non-negative.");
                }
                if (!(damping >= 0.0 && damping < 1.0)) {
                    throw std::invalid_argument("damping must lie in [0, 1).");
                }
                if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
                    throw std::invalid_argument(
                            "tolerance must be finite and positive.");
                }
                if (!(probability_floor > 0.0 && probability_floor < 0.5)
                    || !std::isfinite(probability_floor)) {
                    throw std::invalid_argument(
                            "probability_floor must be finite and lie in (0, 0.5).");
                }
                if (static_cast<int>(split_order.size()) != model.logical_count) {
                    throw std::invalid_argument(
                            "split_order must be a permutation of logical coordinates.");
                }
                auto sorted_order = split_order;
                std::sort(sorted_order.begin(), sorted_order.end());
                for (int index = 0; index < model.logical_count; index++) {
                    if (sorted_order[index] != index) {
                        throw std::invalid_argument(
                                "split_order must be a permutation of logical coordinates.");
                    }
                }
            }

            SparseRows build_generator(const std::vector<int> &remaining) const {
                SparseRows generator(model.bit_count, std::vector<int>{});
                for (int row = 0; row < model.bit_count; row++) {
                    generator[row] = model.D[row];
                    for (int local = 0;
                         local < static_cast<int>(remaining.size());
                         local++) {
                        if (std::binary_search(
                                model.Q[row].begin(),
                                model.Q[row].end(),
                                remaining[local])) {
                            generator[row].push_back(model.degeneracy_dimension + local);
                        }
                    }
                }
                return generator;
            }

            std::vector<int> coordinate_ids(
                    const std::vector<int> &remaining) const {
                std::vector<int> ids(model.degeneracy_dimension + remaining.size());
                std::iota(ids.begin(), ids.begin() + model.degeneracy_dimension, 0);
                for (int local = 0;
                     local < static_cast<int>(remaining.size());
                     local++) {
                    ids[model.degeneracy_dimension + local] =
                            model.degeneracy_dimension + remaining[local];
                }
                return ids;
            }

            QuotientSplitState run_internal_bp(QuotientSplitState state) const {
                const auto generator = build_generator(state.remaining);
                const auto ids = coordinate_ids(state.remaining);
                GeneratorGraphBp graph(
                        generator,
                        model.degeneracy_dimension
                                + static_cast<int>(state.remaining.size()),
                        state.offset,
                        channel_probabilities,
                        ids);
                const GeneratorMessageState *initial = state.bp_result
                        ? &state.bp_result->message_state
                        : nullptr;
                state.bp_result = std::make_shared<GeneratorBpResult>(graph.run(
                        internal_iterations, damping, tolerance, initial));
                return state;
            }

            int choose_split_coordinate(
                    const std::vector<QuotientSplitState> &active,
                    int depth) const {
                if (!adaptive_split) {
                    return split_order[depth];
                }
                if (active.empty()) {
                    throw std::logic_error("Cannot split an empty beam.");
                }

                std::vector<int> common = active[0].remaining;
                for (int state_index = 1;
                     state_index < static_cast<int>(active.size());
                     state_index++) {
                    std::vector<int> intersection;
                    std::set_intersection(
                            common.begin(), common.end(),
                            active[state_index].remaining.begin(),
                            active[state_index].remaining.end(),
                            std::back_inserter(intersection));
                    common = std::move(intersection);
                }
                if (common.empty()) {
                    throw std::logic_error(
                            "No common remaining logical coordinate.");
                }

                double maximum_score = -std::numeric_limits<double>::infinity();
                for (const auto &state: active) {
                    maximum_score = std::max(maximum_score, state.path_score);
                }
                std::vector<double> weights(active.size(), 0.0);
                double total_weight = 0.0;
                for (int index = 0; index < static_cast<int>(active.size()); index++) {
                    weights[index] = std::exp(active[index].path_score - maximum_score);
                    total_weight += weights[index];
                }
                for (double &weight: weights) {
                    weight /= total_weight;
                }

                int best_coordinate = -1;
                double best_entropy = std::numeric_limits<double>::infinity();
                for (int logical_index: common) {
                    double average_entropy = 0.0;
                    for (int state_index = 0;
                         state_index < static_cast<int>(active.size());
                         state_index++) {
                        const auto &state = active[state_index];
                        if (!state.bp_result) {
                            throw std::logic_error("Missing internal BP result.");
                        }
                        const auto position = std::find(
                                state.remaining.begin(),
                                state.remaining.end(),
                                logical_index);
                        const int local = model.degeneracy_dimension
                                + static_cast<int>(position - state.remaining.begin());
                        const double probability_one = clip_probability(
                                state.bp_result->variable_probability_one[local]);
                        const double probability_zero = 1.0 - probability_one;
                        average_entropy += weights[state_index] * (
                                -probability_zero
                                * std::log(std::max(probability_zero, probability_floor))
                                -probability_one
                                * std::log(std::max(probability_one, probability_floor)));
                    }
                    if (average_entropy < best_entropy) {
                        best_entropy = average_entropy;
                        best_coordinate = logical_index;
                    }
                }
                return best_coordinate;
            }

            std::array<QuotientSplitState, 2> branch(
                    const QuotientSplitState &parent,
                    int logical_index) const {
                if (!parent.bp_result) {
                    throw std::logic_error("Cannot branch without BP results.");
                }
                const auto position = std::find(
                        parent.remaining.begin(),
                        parent.remaining.end(),
                        logical_index);
                if (position == parent.remaining.end()) {
                    throw std::logic_error("Split coordinate is not active.");
                }
                const int local = model.degeneracy_dimension
                        + static_cast<int>(position - parent.remaining.begin());
                const double probability_one = clip_probability(
                        parent.bp_result->variable_probability_one[local]);
                const std::array<double, 2> probabilities{
                        1.0 - probability_one, probability_one};

                std::vector<int> remaining;
                remaining.reserve(parent.remaining.size() - 1);
                for (int coordinate: parent.remaining) {
                    if (coordinate != logical_index) {
                        remaining.push_back(coordinate);
                    }
                }

                std::array<QuotientSplitState, 2> children;
                for (int bit = 0; bit < 2; bit++) {
                    auto &child = children[bit];
                    child.logical_bits = parent.logical_bits;
                    child.logical_bits[logical_index] = static_cast<uint8_t>(bit);
                    child.remaining = remaining;
                    child.offset = parent.offset;
                    if (bit != 0) {
                        for (int row = 0; row < model.bit_count; row++) {
                            if (std::binary_search(
                                    model.Q[row].begin(),
                                    model.Q[row].end(),
                                    logical_index)) {
                                child.offset[row] ^= 1;
                            }
                        }
                    }
                    child.path_score = parent.path_score + std::log(std::max(
                            probabilities[bit], probability_floor));
                    child.bp_result = parent.bp_result;
                }
                return children;
            }

            QuotientSplitLeafResult evaluate_leaf(
                    const QuotientSplitState &state) const {
                const std::vector<int> no_logical_coordinates;
                GeneratorGraphBp graph(
                        model.D,
                        model.degeneracy_dimension,
                        state.offset,
                        channel_probabilities,
                        coordinate_ids(no_logical_coordinates));
                const GeneratorMessageState *initial = state.bp_result
                        ? &state.bp_result->message_state
                        : nullptr;
                const auto bp_result = graph.run(
                        leaf_iterations, damping, tolerance, initial);

                QuotientSplitLeafResult leaf;
                leaf.logical_label = state.logical_bits;
                leaf.correction = state.offset;
                leaf.path_score = state.path_score;
                leaf.log_partition = bp_result.log_partition;
                leaf.bp_converged = bp_result.converged;
                leaf.bp_iterations = bp_result.iterations;
                return leaf;
            }

            double clip_probability(double probability) const {
                if (!std::isfinite(probability)) {
                    throw std::runtime_error(
                            "Internal BP returned a non-finite marginal.");
                }
                return std::max(
                        probability_floor,
                        std::min(1.0 - probability_floor, probability));
            }
        };
    }
}  // namespace ldpc_modified::quotient_split

#endif
