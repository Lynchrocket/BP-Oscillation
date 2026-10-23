# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++

from typing import List, Optional, Sequence, Union

import numpy as np
cimport numpy as np
import scipy.sparse

from libcpp.vector cimport vector

import ldpc_modified.helpers.scipy_helpers


cdef BpSparse* _python_to_bp_sparse(matrix):
    """Convert a dense/SciPy binary matrix into the existing C++ GF2Sparse."""
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


cdef vector[int] _build_split_order(int logical_count, object value):
    cdef vector[int] output
    cdef int index
    if value is None:
        output.resize(logical_count)
        for index in range(logical_count):
            output[index] = index
        return output

    order = np.asarray(value)
    if order.ndim != 1 or order.size != logical_count:
        raise ValueError(
            f"split_order must be a permutation of range({logical_count})."
        )
    if set(int(item) for item in order.tolist()) != set(range(logical_count)):
        raise ValueError(
            f"split_order must be a permutation of range({logical_count})."
        )
    output.resize(logical_count)
    for index in range(logical_count):
        output[index] = int(order[index])
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


cdef np.ndarray _leaf_labels_to_numpy(
    const vector[QuotientSplitLeafResultCpp]& leaves,
    int logical_count,
):
    cdef int row_count = leaves.size()
    cdef np.ndarray[np.uint8_t, ndim=2] output = np.zeros(
        (row_count, logical_count), dtype=np.uint8
    )
    cdef int row
    cdef int column
    for row in range(row_count):
        if leaves[row].logical_label.size() != logical_count:
            raise RuntimeError("C++ returned a logical label with an invalid length.")
        for column in range(logical_count):
            output[row, column] = leaves[row].logical_label[column]
    return output


cdef np.ndarray _leaf_corrections_to_numpy(
    const vector[QuotientSplitLeafResultCpp]& leaves,
    int bit_count,
):
    cdef int row_count = leaves.size()
    cdef np.ndarray[np.uint8_t, ndim=2] output = np.zeros(
        (row_count, bit_count), dtype=np.uint8
    )
    cdef int row
    cdef int column
    for row in range(row_count):
        if leaves[row].correction.size() != bit_count:
            raise RuntimeError("C++ returned a correction with an invalid length.")
        for column in range(bit_count):
            output[row, column] = leaves[row].correction[column]
    return output


cdef np.ndarray _leaf_double_field_to_numpy(
    const vector[QuotientSplitLeafResultCpp]& leaves,
    bint log_partition,
):
    cdef int size = leaves.size()
    cdef np.ndarray[np.float64_t, ndim=1] output = np.zeros(
        size, dtype=np.float64
    )
    cdef int index
    for index in range(size):
        output[index] = (
            leaves[index].log_partition
            if log_partition
            else leaves[index].path_score
        )
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


cdef class QuotientSplitDecoder:
    """Quotient successive-cancellation list BP decoder.

    C++ first constructs ``T/K/Q/D`` so that, for the logical-neutral root
    ``x_s``, every syndrome-compatible error has the form

        x = x_s XOR D u XOR Q ell.

    Internal BP runs on ``[D, Q_remaining]``.  Each split fixes one logical bit,
    updates the affine offset, and keeps at most ``list_width`` prefixes by their
    accumulated log marginal.  Surviving leaves are finally ranked by D-only
    Bethe log-partition estimates.
    """

    def __cinit__(
        self,
        H,
        L,
        error_rate: Optional[float] = None,
        error_channel: Optional[Union[np.ndarray, List[float]]] = None,
        list_width: int = 2,
        split_order: Optional[Sequence[int]] = None,
        adaptive_split: bool = False,
        internal_iterations: int = 5,
        leaf_iterations: int = 50,
        damping: float = 0.2,
        tolerance: float = 1e-9,
        probability_floor: float = 1e-300,
    ):
        cdef vector[double] probabilities
        cdef vector[int] parsed_split_order

        self.h_pcm = NULL
        self.logical_pcm = NULL
        self.model = NULL
        self.split_decoder = NULL

        if not hasattr(H, "shape") or len(H.shape) != 2:
            raise ValueError("H must be a two-dimensional matrix.")
        if not hasattr(L, "shape") or len(L.shape) != 2:
            raise ValueError("L must be a two-dimensional matrix.")
        if H.shape[1] != L.shape[1]:
            raise ValueError("H and L must have the same number of columns.")
        if list_width is True or list_width is False or not isinstance(
            list_width, (int, np.integer)
        ) or list_width < 1:
            raise ValueError("list_width must be a positive integer.")
        if internal_iterations is True or internal_iterations is False or not isinstance(
            internal_iterations, (int, np.integer)
        ) or internal_iterations < 0:
            raise ValueError("internal_iterations must be a non-negative integer.")
        if leaf_iterations is True or leaf_iterations is False or not isinstance(
            leaf_iterations, (int, np.integer)
        ) or leaf_iterations < 0:
            raise ValueError("leaf_iterations must be a non-negative integer.")

        self.m = H.shape[0]
        self.n = H.shape[1]
        self.k = L.shape[0]
        self._syndrome.resize(self.m)
        probabilities = _build_error_channel(
            self.n, error_rate, error_channel
        )
        parsed_split_order = _build_split_order(self.k, split_order)

        self.h_pcm = _python_to_bp_sparse(H)
        self.logical_pcm = _python_to_bp_sparse(L)
        self.model = new QuotientModelCpp(
            self.h_pcm[0], self.logical_pcm[0]
        )
        if self.model.logical_rank != self.k:
            raise ValueError(
                "QuotientSplitDecoder requires full logical rank and LQ = I."
            )
        self.split_decoder = new QuotientSplitDecoderCpp(
            self.model[0],
            probabilities,
            list_width,
            parsed_split_order,
            adaptive_split,
            internal_iterations,
            leaf_iterations,
            damping,
            tolerance,
            probability_floor,
        )
        self._error_channel = probabilities

    def __init__(
        self,
        H,
        L,
        error_rate: Optional[float] = None,
        error_channel: Optional[Union[np.ndarray, List[float]]] = None,
        list_width: int = 2,
        split_order: Optional[Sequence[int]] = None,
        adaptive_split: bool = False,
        internal_iterations: int = 5,
        leaf_iterations: int = 50,
        damping: float = 0.2,
        tolerance: float = 1e-9,
        probability_floor: float = 1e-300,
    ):
        pass

    def __dealloc__(self):
        if self.split_decoder != NULL:
            del self.split_decoder
        if self.model != NULL:
            del self.model
        if self.logical_pcm != NULL:
            del self.logical_pcm
        if self.h_pcm != NULL:
            del self.h_pcm

    def decode(self, syndrome: np.ndarray) -> np.ndarray:
        """Decode a binary syndrome into the highest-ranked leaf correction."""
        cdef int check
        values = np.asarray(syndrome)
        if values.ndim != 1 or values.size != self.m:
            raise ValueError(
                f"syndrome must be one-dimensional with length {self.m}."
            )
        if np.any((values != 0) & (values != 1)):
            raise ValueError("syndrome must be binary.")
        for check in range(self.m):
            self._syndrome[check] = values[check]

        self._result = self.split_decoder.decode(self._syndrome)
        return _uint8_vector_to_numpy(self._result.correction)

    def update_channel_probs(self, error_channel) -> None:
        cdef vector[double] probabilities = _build_error_channel(
            self.n, None, error_channel
        )
        self.split_decoder.set_channel_probabilities(probabilities)
        self._error_channel = probabilities

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
        return _uint8_vector_to_numpy(self._result.correction)

    @property
    def best_logical_label(self) -> np.ndarray:
        return _uint8_vector_to_numpy(self._result.best_logical_label)

    @property
    def leaf_logical_labels(self) -> np.ndarray:
        return _leaf_labels_to_numpy(self._result.leaves, self.k)

    @property
    def leaf_corrections(self) -> np.ndarray:
        return _leaf_corrections_to_numpy(self._result.leaves, self.n)

    @property
    def leaf_path_scores(self) -> np.ndarray:
        return _leaf_double_field_to_numpy(self._result.leaves, False)

    @property
    def leaf_log_partitions(self) -> np.ndarray:
        return _leaf_double_field_to_numpy(self._result.leaves, True)

    @property
    def leaf_bp_converged(self) -> np.ndarray:
        cdef int size = self._result.leaves.size()
        cdef np.ndarray[np.uint8_t, ndim=1] output = np.zeros(
            size, dtype=np.uint8
        )
        cdef int index
        for index in range(size):
            output[index] = 1 if self._result.leaves[index].bp_converged else 0
        return output.astype(np.bool_)

    @property
    def leaf_bp_iterations(self) -> np.ndarray:
        cdef int size = self._result.leaves.size()
        cdef np.ndarray[np.int64_t, ndim=1] output = np.zeros(
            size, dtype=np.int64
        )
        cdef int index
        for index in range(size):
            output[index] = self._result.leaves[index].bp_iterations
        return output

    @property
    def surviving_prefix_counts(self) -> np.ndarray:
        return _int_vector_to_numpy(self._result.surviving_prefix_counts)

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
    def list_width(self) -> int:
        return self.split_decoder.list_width

    @property
    def split_order(self) -> np.ndarray:
        return _int_vector_to_numpy(self.split_decoder.split_order)

    @property
    def adaptive_split(self) -> bool:
        return self.split_decoder.adaptive_split

    @property
    def internal_iterations(self) -> int:
        return self.split_decoder.internal_iterations

    @property
    def leaf_iterations(self) -> int:
        return self.split_decoder.leaf_iterations

    @property
    def damping(self) -> float:
        return self.split_decoder.damping

    @property
    def tolerance(self) -> float:
        return self.split_decoder.tolerance

    @property
    def probability_floor(self) -> float:
        return self.split_decoder.probability_floor
