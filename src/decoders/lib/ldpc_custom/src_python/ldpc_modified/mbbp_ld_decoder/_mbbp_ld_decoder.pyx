#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
import numpy as np
import scipy.sparse
from typing import List, Optional, Union

import ldpc_modified.helpers.scipy_helpers


cdef BpSparse* Py2BpSparse(pcm):
    cdef int row
    cdef int column
    cdef int nonzero_count
    cdef unique_ptr[BpSparse] cpcm

    if not isinstance(pcm, (np.ndarray, scipy.sparse.spmatrix)):
        raise TypeError(
            "pcm must be a numpy array or scipy sparse matrix, "
            f"not {type(pcm)}"
        )
    if len(pcm.shape) != 2:
        raise ValueError("pcm must be two-dimensional")
    if isinstance(pcm, np.ndarray):
        if not np.all(np.logical_or(pcm == 0, pcm == 1)):
            raise ValueError("pcm must contain only binary values")
        pcm = np.asarray(pcm, dtype=np.uint8)
    else:
        sparse_pcm = pcm.tocoo(copy=True)
        if not np.all(np.logical_or(sparse_pcm.data == 0, sparse_pcm.data == 1)):
            raise ValueError("pcm must contain only binary values")
        coordinates = (
            sparse_pcm.row.astype(np.int64) * sparse_pcm.shape[1]
            + sparse_pcm.col.astype(np.int64)
        )
        if np.unique(coordinates).size != sparse_pcm.nnz:
            raise ValueError("pcm sparse storage must not contain duplicate entries")
        pcm = sparse_pcm.tocsr().astype(np.uint8)
    pcm = ldpc_modified.helpers.scipy_helpers.convert_to_binary_sparse(pcm)
    nonzero_count = int(np.count_nonzero(pcm)) if isinstance(pcm, np.ndarray) else int(pcm.nnz)
    cpcm = unique_ptr[BpSparse](
        new BpSparse(pcm.shape[0], pcm.shape[1], nonzero_count)
    )

    if isinstance(pcm, np.ndarray):
        for row in range(pcm.shape[0]):
            for column in range(pcm.shape[1]):
                if pcm[row, column]:
                    cpcm.get().insert_entry(row, column)
    else:
        rows, columns = pcm.nonzero()
        for row, column in zip(rows, columns):
            cpcm.get().insert_entry(row, column)
    return cpcm.release()


cdef class MbbpLdDecoder:
    """Multiple-Bases Belief-Propagation List Decoder."""

    def __cinit__(
        self,
        pcm,
        error_rate: Optional[float] = 0.001,
        error_channel: Optional[Union[np.ndarray, List[float]]] = None,
        max_iter: int = 100,
        bp_method: str = "minimum_sum",
        ms_scaling_factor: float = 0.875,
        schedule: str = "serial",
        check_permutations=None,
    ):
        cdef vector[double] channel_probabilities
        cdef vector[vector[int]] native_permutations
        cdef int row
        cdef int column
        cdef int bp_method_value
        cdef int schedule_value

        self.pcm = NULL
        self.decoder = NULL

        if error_rate is None and error_channel is None:
            raise ValueError("Specify either error_rate or error_channel")
        if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
            raise ValueError("max_iter must be a non-negative integer")
        if not isinstance(ms_scaling_factor, (float, int, np.floating, np.integer)):
            raise TypeError("ms_scaling_factor must be numeric")
        if not np.isfinite(ms_scaling_factor) or ms_scaling_factor < 0:
            raise ValueError("ms_scaling_factor must be finite and non-negative")

        self.pcm = Py2BpSparse(pcm)
        self.m = self.pcm.m
        self.n = self.pcm.n
        self._syndrome.resize(self.m)

        channel_probabilities.resize(self.n)
        if error_channel is not None:
            probabilities = np.asarray(error_channel, dtype=float).reshape(-1)
            if probabilities.size != self.n:
                raise ValueError(
                    f"error_channel must have length {self.n}, not {probabilities.size}"
                )
            if not np.all(np.isfinite(probabilities)) or np.any(probabilities <= 0) or np.any(probabilities >= 1):
                raise ValueError("Every error_channel probability must lie strictly between 0 and 1")
            for column in range(self.n):
                channel_probabilities[column] = probabilities[column]
        else:
            if not isinstance(error_rate, (float, int, np.floating, np.integer)):
                raise TypeError("error_rate must be numeric")
            if not np.isfinite(error_rate) or error_rate <= 0 or error_rate >= 1:
                raise ValueError("error_rate must lie strictly between 0 and 1")
            for column in range(self.n):
                channel_probabilities[column] = error_rate

        method = str(bp_method).lower()
        if method in ("product_sum", "prod_sum", "ps", "0"):
            bp_method_value = 0
        elif method in ("minimum_sum", "min_sum", "ms", "1"):
            bp_method_value = 1
        else:
            raise ValueError("bp_method must be 'product_sum' or 'minimum_sum'")

        schedule_name = str(schedule).lower()
        if schedule_name in ("serial", "s", "0"):
            schedule_value = 0
        elif schedule_name in ("parallel", "p", "1"):
            schedule_value = 1
        else:
            raise ValueError("schedule must be 'serial' or 'parallel'")

        if check_permutations is None:
            permutations = np.arange(self.m, dtype=np.int64).reshape(1, self.m)
        else:
            permutations = np.asarray(check_permutations)
            if not np.issubdtype(permutations.dtype, np.integer):
                raise ValueError("check_permutations must contain integer indices")
            if permutations.ndim == 1:
                permutations = permutations.reshape(1, -1)
            if permutations.ndim != 2 or permutations.shape[1] != self.m:
                raise ValueError(
                    f"check_permutations must have shape ({self.m},) or (P, {self.m})"
                )
            if permutations.shape[0] == 0:
                raise ValueError("At least one check permutation is required")
            for row in range(permutations.shape[0]):
                if not np.array_equal(np.sort(permutations[row]), np.arange(self.m)):
                    raise ValueError("Each row of check_permutations must be a permutation of check indices")

        native_permutations.resize(permutations.shape[0])
        for row in range(permutations.shape[0]):
            native_permutations[row].resize(self.m)
            for column in range(self.m):
                native_permutations[row][column] = permutations[row, column]

        self.decoder = new MbbpLdDecoderCpp(
            self.pcm[0],
            channel_probabilities,
            max_iter,
            bp_method_value,
            schedule_value,
            ms_scaling_factor,
            native_permutations,
        )

    def __init__(
        self,
        pcm,
        error_rate: Optional[float] = 0.001,
        error_channel: Optional[Union[np.ndarray, List[float]]] = None,
        max_iter: int = 100,
        bp_method: str = "minimum_sum",
        ms_scaling_factor: float = 0.875,
        schedule: str = "serial",
        check_permutations=None,
    ):
        pass

    def __dealloc__(self):
        if self.decoder != NULL:
            del self.decoder
        if self.pcm != NULL:
            del self.pcm

    def decode(self, syndrome) -> np.ndarray:
        syndrome_array = np.asarray(syndrome, dtype=np.uint8).reshape(-1) % 2
        if syndrome_array.size != self.m:
            raise ValueError(
                f"syndrome must have length {self.m}, not {syndrome_array.size}"
            )
        for index in range(self.m):
            self._syndrome[index] = syndrome_array[index]
        self.decoder.decode(self._syndrome)
        return self.decoding

    @property
    def decoding(self) -> np.ndarray:
        out = np.zeros(self.n, dtype=np.uint8)
        for index in range(self.n):
            out[index] = self.decoder.decoding[index]
        return out

    @property
    def converge(self) -> bool:
        return True if self.decoder.converge else False

    @property
    def basis_count(self) -> int:
        return self.decoder.basis_count()

    @property
    def partition_count(self) -> int:
        return self.decoder.partition_count()

    @property
    def candidate_count(self) -> int:
        return self.decoder.candidate_count

    @property
    def unique_candidate_count(self) -> int:
        return self.decoder.unique_candidate_count

    @property
    def selected_frequency(self) -> int:
        return self.decoder.selected_frequency

    @property
    def selected_weight(self) -> int:
        return self.decoder.selected_weight

    @property
    def partition_rows(self) -> list:
        return [
            [
                [
                    self.decoder.partition_subtree_row(partition, subtree, row)
                    for row in range(self.decoder.partition_subtree_size(partition, subtree))
                ]
                for subtree in range(self.decoder.partition_subtree_count(partition))
            ]
            for partition in range(self.decoder.partition_count())
        ]

    @property
    def basis_rows(self) -> list:
        return [
            [
                self.decoder.basis_source_row(basis, row)
                for row in range(self.decoder.basis_row_count(basis))
            ]
            for basis in range(self.decoder.basis_count())
        ]

    @property
    def basis_candidates(self) -> np.ndarray:
        out = np.zeros((self.decoder.basis_count(), self.n), dtype=np.uint8)
        for basis in range(self.decoder.basis_count()):
            for bit in range(self.n):
                out[basis, bit] = self.decoder.basis_candidate_bit(basis, bit)
        return out

    @property
    def basis_convergence(self) -> np.ndarray:
        out = np.zeros(self.decoder.basis_count(), dtype=np.bool_)
        for basis in range(self.decoder.basis_count()):
            out[basis] = self.decoder.basis_converged_at(basis)
        return out

    @property
    def basis_iterations(self) -> np.ndarray:
        out = np.zeros(self.decoder.basis_count(), dtype=np.int64)
        for basis in range(self.decoder.basis_count()):
            out[basis] = self.decoder.basis_iterations_at(basis)
        return out

    @property
    def basis_edge_counts(self) -> np.ndarray:
        out = np.zeros(self.decoder.basis_count(), dtype=np.int64)
        for basis in range(self.decoder.basis_count()):
            out[basis] = self.decoder.basis_edge_count_at(basis)
        return out

    @property
    def converged_candidates(self) -> np.ndarray:
        out = np.zeros((self.decoder.candidate_count, self.n), dtype=np.uint8)
        for candidate in range(self.decoder.candidate_count):
            for bit in range(self.n):
                out[candidate, bit] = self.decoder.converged_candidate_bit(candidate, bit)
        return out

    @property
    def unique_candidates(self) -> np.ndarray:
        out = np.zeros((self.decoder.unique_candidate_count, self.n), dtype=np.uint8)
        for candidate in range(self.decoder.unique_candidate_count):
            for bit in range(self.n):
                out[candidate, bit] = self.decoder.unique_candidate_bit(candidate, bit)
        return out

    @property
    def unique_candidate_frequencies(self) -> np.ndarray:
        out = np.zeros(self.decoder.unique_candidate_count, dtype=np.int64)
        for candidate in range(self.decoder.unique_candidate_count):
            out[candidate] = self.decoder.unique_candidate_frequency_at(candidate)
        return out

    @property
    def unique_candidate_weights(self) -> np.ndarray:
        out = np.zeros(self.decoder.unique_candidate_count, dtype=np.int64)
        for candidate in range(self.decoder.unique_candidate_count):
            out[candidate] = self.decoder.unique_candidate_weight_at(candidate)
        return out

    @property
    def channel_probabilities(self) -> np.ndarray:
        out = np.zeros(self.n, dtype=float)
        for bit in range(self.n):
            out[bit] = self.decoder.channel_probability(bit)
        return out
