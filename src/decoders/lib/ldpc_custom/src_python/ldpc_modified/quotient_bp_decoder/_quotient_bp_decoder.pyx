# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++

from typing import List, Optional, Union

import numpy as np
cimport numpy as np
import scipy.sparse

from libcpp.vector cimport vector

import ldpc_modified.helpers.scipy_helpers


cdef BpSparse* _python_to_bp_sparse(matrix):
    """Convert a dense/SciPy matrix into the existing C++ GF2Sparse type."""
    if not isinstance(matrix, np.ndarray) and not scipy.sparse.issparse(matrix):
        raise TypeError(
            "The input matrix must be a numpy.ndarray or scipy sparse matrix, "
            f"not {type(matrix)}."
        )

    matrix = ldpc_modified.helpers.scipy_helpers.convert_to_binary_sparse(matrix)
    cdef int row_count = matrix.shape[0]
    cdef int column_count = matrix.shape[1]
    cdef int nonzero_count
    cdef int row
    cdef int column

    if isinstance(matrix, np.ndarray):
        nonzero_count = int(np.count_nonzero(matrix))
    else:
        nonzero_count = int(matrix.nnz)

    cdef BpSparse* output = new BpSparse(
        row_count, column_count, nonzero_count
    )
    rows, columns = matrix.nonzero()
    for index in range(len(rows)):
        row = rows[index]
        column = columns[index]
        output.insert_entry(row, column)
    return output


cdef np.ndarray _uint8_vector_to_numpy(const vector[uint8_t]& values):
    cdef int size = values.size()
    cdef np.ndarray[np.uint8_t, ndim=1] output = np.zeros(size, dtype=np.uint8)
    cdef int index
    for index in range(size):
        output[index] = values[index]
    return output


cdef np.ndarray _double_vector_to_numpy(const vector[double]& values):
    cdef int size = values.size()
    cdef np.ndarray[np.float64_t, ndim=1] output = np.zeros(size, dtype=np.float64)
    cdef int index
    for index in range(size):
        output[index] = values[index]
    return output


cdef np.ndarray _int_vector_to_numpy(const vector[int]& values):
    cdef int size = values.size()
    cdef np.ndarray[np.int64_t, ndim=1] output = np.zeros(size, dtype=np.int64)
    cdef int index
    for index in range(size):
        output[index] = values[index]
    return output


cdef np.ndarray _uint8_matrix_to_numpy(
    const vector[vector[uint8_t]]& values,
    int column_count,
):
    cdef int row_count = values.size()
    cdef np.ndarray[np.uint8_t, ndim=2] output = np.zeros(
        (row_count, column_count), dtype=np.uint8
    )
    cdef int row
    cdef int column
    for row in range(row_count):
        if values[row].size() != column_count:
            raise RuntimeError("C++ returned a ragged binary matrix.")
        for column in range(column_count):
            output[row, column] = values[row][column]
    return output


cdef object _sparse_rows_to_scipy(
    const vector[vector[int]]& rows,
    int column_count,
):
    cdef int row_count = rows.size()
    cdef int entry_count = 0
    cdef int row
    cdef int index
    cdef int position = 0
    for row in range(row_count):
        entry_count += rows[row].size()

    cdef np.ndarray[np.int32_t, ndim=1] row_indices = np.zeros(
        entry_count, dtype=np.int32
    )
    cdef np.ndarray[np.int32_t, ndim=1] column_indices = np.zeros(
        entry_count, dtype=np.int32
    )
    cdef np.ndarray[np.uint8_t, ndim=1] data = np.ones(
        entry_count, dtype=np.uint8
    )
    for row in range(row_count):
        for index in range(rows[row].size()):
            row_indices[position] = row
            column_indices[position] = rows[row][index]
            position += 1
    return scipy.sparse.csr_matrix(
        (data, (row_indices, column_indices)),
        shape=(row_count, column_count),
        dtype=np.uint8,
    )


cdef BpMethod _parse_bp_method(object value):
    if isinstance(value, str):
        normalized = value.lower().replace("-", "_")
        if normalized in {"product_sum", "sum_product", "ps", "0"}:
            return PRODUCT_SUM
        if normalized in {"minimum_sum", "min_sum", "ms", "1"}:
            return MINIMUM_SUM
    elif isinstance(value, (int, np.integer)):
        if value == 0:
            return PRODUCT_SUM
        if value == 1:
            return MINIMUM_SUM
    raise ValueError(
        "bp_method must be product_sum/sum_product or minimum_sum/min_sum."
    )


cdef BpSchedule _parse_schedule(object value):
    if isinstance(value, str):
        normalized = value.lower().replace("-", "_")
        # Match the public BpDecoder API: its Python integer aliases describe
        # the user-facing order, which differs from the C++ enum values.
        if normalized in {"parallel", "p", "0"}:
            return PARALLEL
        if normalized in {"serial", "s", "1"}:
            return SERIAL
        if normalized in {"serial_relative", "relative", "sr", "2"}:
            return SERIAL_RELATIVE
    elif isinstance(value, (int, np.integer)):
        if value == 0:
            return PARALLEL
        if value == 1:
            return SERIAL
        if value == 2:
            return SERIAL_RELATIVE
    raise ValueError("schedule must be parallel, serial, or serial_relative.")


cdef vector[double] _build_error_channel(
    int bit_count,
    object error_rate,
    object error_channel,
):
    if error_rate is not None and error_channel is not None:
        raise ValueError("Specify either error_rate or error_channel, not both.")
    if error_rate is None and error_channel is None:
        raise ValueError("Specify either error_rate or error_channel.")

    cdef vector[double] probabilities
    probabilities.resize(bit_count)
    cdef int bit
    cdef double probability

    if error_channel is not None:
        channel = np.asarray(error_channel, dtype=np.float64)
        if channel.ndim != 1 or channel.size != bit_count:
            raise ValueError(
                f"error_channel must be one-dimensional with length {bit_count}."
            )
        for bit in range(bit_count):
            probability = channel[bit]
            if not np.isfinite(probability) or not 0.0 < probability < 1.0:
                raise ValueError(
                    "Every error_channel entry must be finite and lie in (0, 1)."
                )
            probabilities[bit] = probability
        return probabilities

    if error_rate is True or error_rate is False or isinstance(error_rate, np.bool_):
        raise TypeError("error_rate must be a real scalar, not bool.")
    if not isinstance(error_rate, (int, float, np.integer, np.floating)):
        raise TypeError("error_rate must be a real scalar.")
    probability = error_rate
    if not np.isfinite(probability) or not 0.0 < probability < 1.0:
        raise ValueError("error_rate must be finite and lie in (0, 1).")
    for bit in range(bit_count):
        probabilities[bit] = probability
    return probabilities


cdef vector[int] _build_serial_schedule(int bit_count, object value):
    cdef vector[int] output
    cdef int index
    cdef int bit
    if value is None:
        return NULL_INT_VECTOR
    schedule = np.asarray(value)
    if schedule.ndim != 1 or schedule.size != bit_count:
        raise ValueError(
            f"serial_schedule_order must have length {bit_count}."
        )
    if set(int(item) for item in schedule.tolist()) != set(range(bit_count)):
        raise ValueError(
            "serial_schedule_order must be a permutation of range(bit_count)."
        )
    output.resize(bit_count)
    for index in range(bit_count):
        bit = schedule[index]
        output[index] = bit
    return output


cdef class QuotientBpDecoder:
    """Trajectory-screened quotient belief-propagation decoder.

    C++ computes the quotient decomposition ``T/K/Q/D``.  The first BP run on
    ``H`` stores the unique logical labels ``L @ hard_decision`` visited by its
    full iteration trajectory.  Only those labels are reranked by sum-product
    BP on ``D`` using their Bethe log-partition estimates.
    """

    def __cinit__(
        self,
        H,
        L,
        error_rate: Optional[float] = None,
        error_channel: Optional[Union[np.ndarray, List[float]]] = None,
        max_iter: int = 0,
        bp_method: Union[str, int] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        schedule: Union[str, int] = "parallel",
        omp_thread_count: int = 1,
        random_schedule_seed: int = 0,
        serial_schedule_order=None,
        random_serial_schedule=False,
        d_max_iter: int = 100,
        d_damping: float = 0.2,
        d_tolerance: float = 1e-10,
    ):
        cdef vector[double] probabilities
        cdef BpMethod parsed_method
        cdef BpSchedule parsed_schedule
        cdef int effective_max_iter

        self.h_pcm = NULL
        self.logical_pcm = NULL
        self.model = NULL
        self.h_bpd = NULL
        self.d_bpd = NULL

        if not hasattr(H, "shape") or len(H.shape) != 2:
            raise ValueError("H must be a two-dimensional matrix.")
        if not hasattr(L, "shape") or len(L.shape) != 2:
            raise ValueError("L must be a two-dimensional matrix.")
        if H.shape[1] != L.shape[1]:
            raise ValueError("H and L must have the same number of columns.")
        if max_iter is True or max_iter is False or not isinstance(
            max_iter, (int, np.integer)
        ) or max_iter < 0:
            raise ValueError("max_iter must be a non-negative integer.")
        if d_max_iter is True or d_max_iter is False or not isinstance(
            d_max_iter, (int, np.integer)
        ) or d_max_iter < 1:
            raise ValueError("d_max_iter must be a positive integer.")
        if omp_thread_count < 1:
            raise ValueError("omp_thread_count must be positive.")
        if ms_scaling_factor < 0:
            raise ValueError("ms_scaling_factor must be non-negative.")
        if random_serial_schedule and serial_schedule_order is not None:
            raise ValueError(
                "random_serial_schedule conflicts with serial_schedule_order."
            )

        self.m = H.shape[0]
        self.n = H.shape[1]
        self.k = L.shape[0]
        self._syndrome.resize(self.m)
        probabilities = _build_error_channel(
            self.n, error_rate, error_channel
        )
        self._error_channel = probabilities
        self._serial_schedule_order = _build_serial_schedule(
            self.n, serial_schedule_order
        )
        parsed_method = _parse_bp_method(bp_method)
        parsed_schedule = _parse_schedule(schedule)
        effective_max_iter = max_iter if max_iter != 0 else self.n

        self.h_pcm = _python_to_bp_sparse(H)
        self.logical_pcm = _python_to_bp_sparse(L)
        self.model = new QuotientModelCpp(
            self.h_pcm[0], self.logical_pcm[0]
        )
        self.h_bpd = new BpDecoderCpp(
            self.h_pcm[0],
            probabilities,
            effective_max_iter,
            parsed_method,
            parsed_schedule,
            ms_scaling_factor,
            omp_thread_count,
            self._serial_schedule_order,
            random_schedule_seed,
            random_serial_schedule,
            SYNDROME,
        )
        self.h_bpd.enable_logical_tracking(self.logical_pcm[0], True)
        self.d_bpd = new DPartitionBpDecoderCpp(
            self.model.d_matrix(),
            probabilities,
            d_max_iter,
            d_damping,
            d_tolerance,
        )

    def __init__(
        self,
        H,
        L,
        error_rate: Optional[float] = None,
        error_channel: Optional[Union[np.ndarray, List[float]]] = None,
        max_iter: int = 60,
        bp_method: Union[str, int] = "minimum_sum",
        ms_scaling_factor: float = 1.0,
        schedule: Union[str, int] = "parallel",
        omp_thread_count: int = 1,
        random_schedule_seed: int = 0,
        serial_schedule_order=None,
        random_serial_schedule=False,
        d_max_iter: int = 60,
        d_damping: float = 0.2,
        d_tolerance: float = 1e-10,
    ):
        pass

    def __dealloc__(self):
        if self.d_bpd != NULL:
            del self.d_bpd
        if self.h_bpd != NULL:
            del self.h_bpd
        if self.model != NULL:
            del self.model
        if self.logical_pcm != NULL:
            del self.logical_pcm
        if self.h_pcm != NULL:
            del self.h_pcm

    def decode(self, syndrome: np.ndarray) -> np.ndarray:
        """Decode a binary syndrome into the selected coset representative."""
        cdef int check
        cdef int candidate_index
        cdef int best_index
        cdef double score
        cdef vector[uint8_t] ignored_h_decoding
        cdef vector[uint8_t] label
        cdef vector[uint8_t] representative

        values = np.asarray(syndrome)
        if values.ndim != 1 or values.size != self.m:
            raise ValueError(
                f"syndrome must be one-dimensional with length {self.m}."
            )
        if np.any((values != 0) & (values != 1)):
            raise ValueError("syndrome must be binary.")
        for check in range(self.m):
            self._syndrome[check] = values[check]

        ignored_h_decoding = self.h_bpd.decode(self._syndrome)
        if self.h_bpd.logical_cosets.empty():
            raise RuntimeError("H-BP produced no logical-coset candidate.")

        self._candidate_representatives.clear()
        self._candidate_log_partitions.clear()
        self._d_iterations.clear()
        self._d_converged.clear()

        for candidate_index in range(self.h_bpd.logical_cosets.size()):
            label = self.h_bpd.logical_cosets[candidate_index]
            representative = self.model.representative(
                self._syndrome, label
            )
            score = self.d_bpd.decode(representative)
            self._candidate_representatives.push_back(representative)
            self._candidate_log_partitions.push_back(score)
            self._d_iterations.push_back(self.d_bpd.iterations)
            self._d_converged.push_back(1 if self.d_bpd.converge else 0)

        best_index = 0
        for candidate_index in range(1, self._candidate_log_partitions.size()):
            if (
                self._candidate_log_partitions[candidate_index]
                > self._candidate_log_partitions[best_index]
            ):
                best_index = candidate_index

        self._decoding = self._candidate_representatives[best_index]
        self._chosen_logical_coset = self.h_bpd.logical_cosets[best_index]
        return _uint8_vector_to_numpy(self._decoding)

    def update_channel_probs(self, error_channel) -> None:
        cdef vector[double] probabilities = _build_error_channel(
            self.n, None, error_channel
        )
        self._error_channel = probabilities
        self.h_bpd.channel_probabilities = probabilities
        self.d_bpd.set_channel_probabilities(probabilities)

    @property
    def error_channel(self) -> np.ndarray:
        return _double_vector_to_numpy(self._error_channel)

    @error_channel.setter
    def error_channel(self, value) -> None:
        self.update_channel_probs(value)

    @property
    def error_rate(self) -> np.ndarray:
        return self.error_channel

    @property
    def decoding(self) -> np.ndarray:
        return _uint8_vector_to_numpy(self._decoding)

    @property
    def chosen_logical_coset(self) -> np.ndarray:
        return _uint8_vector_to_numpy(self._chosen_logical_coset)

    @property
    def trajectory_logical_cosets(self) -> np.ndarray:
        return _uint8_matrix_to_numpy(
            self.h_bpd.trajectory_logical_cosets, self.k
        )

    @property
    def candidate_logical_cosets(self) -> np.ndarray:
        return _uint8_matrix_to_numpy(self.h_bpd.logical_cosets, self.k)

    @property
    def candidate_representatives(self) -> np.ndarray:
        return _uint8_matrix_to_numpy(
            self._candidate_representatives, self.n
        )

    @property
    def candidate_log_partitions(self) -> np.ndarray:
        return _double_vector_to_numpy(self._candidate_log_partitions)

    @property
    def d_iter(self) -> np.ndarray:
        return _int_vector_to_numpy(self._d_iterations)

    @property
    def d_converge(self) -> np.ndarray:
        return _uint8_vector_to_numpy(self._d_converged).astype(np.bool_)

    @property
    def h_log_prob_ratios(self) -> np.ndarray:
        return _double_vector_to_numpy(self.h_bpd.log_prob_ratios)

    @property
    def h_iter(self) -> int:
        return self.h_bpd.iterations

    @property
    def h_converge(self):
        return self.h_bpd.converge

    @property
    def T(self):
        return _sparse_rows_to_scipy(self.model.T, self.m)

    @property
    def K(self):
        return _sparse_rows_to_scipy(
            self.model.K, self.model.kernel_dimension
        )

    @property
    def Q(self):
        return _sparse_rows_to_scipy(
            self.model.Q, self.model.logical_rank
        )

    @property
    def D(self):
        return _sparse_rows_to_scipy(
            self.model.D, self.model.degeneracy_dimension
        )

    @property
    def logical_map(self):
        return _sparse_rows_to_scipy(
            self.model.logical_map, self.model.logical_rank
        )

    @property
    def check_count(self) -> int:
        return self.m

    @property
    def bit_count(self) -> int:
        return self.n

    @property
    def logical_count(self) -> int:
        return self.k

    @property
    def rank_h(self) -> int:
        return self.model.rank_h

    @property
    def logical_rank(self) -> int:
        return self.model.logical_rank

    @property
    def degeneracy_dimension(self) -> int:
        return self.model.degeneracy_dimension

    @property
    def max_iter(self) -> int:
        return self.h_bpd.maximum_iterations

    @max_iter.setter
    def max_iter(self, value) -> None:
        if value is True or value is False or not isinstance(
            value, (int, np.integer)
        ) or value < 1:
            raise ValueError("max_iter must be a positive integer.")
        self.h_bpd.maximum_iterations = value

    @property
    def d_max_iter(self) -> int:
        return self.d_bpd.maximum_iterations

    @d_max_iter.setter
    def d_max_iter(self, value) -> None:
        if value is True or value is False or not isinstance(
            value, (int, np.integer)
        ) or value < 1:
            raise ValueError("d_max_iter must be a positive integer.")
        self.d_bpd.maximum_iterations = value
