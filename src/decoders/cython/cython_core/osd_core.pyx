# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

from itertools import combinations

import numpy as np
cimport numpy as cnp
from libc.math cimport log


OSD_VERSION = "cpp_cs_v1"


def osd_candidate_count(information_bits, order):
    assert information_bits >= 0 and order >= 0
    width = min(int(order), int(information_bits))
    return 1 if order == 0 else 1 + information_bits + width * (width - 1) // 2


ctypedef cnp.uint8_t uint8_t
ctypedef cnp.int32_t int32_t


cdef cnp.ndarray[cnp.int32_t, ndim=1] _select_info_set(cnp.ndarray[cnp.uint8_t, ndim=2] H_arr):
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] work_arr = np.ascontiguousarray(H_arr.copy() % 2, dtype=np.uint8)
    cdef uint8_t[:, ::1] work = work_arr
    cdef int m = work_arr.shape[0]
    cdef int n = work_arr.shape[1]
    cdef list pivots = []
    cdef int row = 0
    cdef int col, r, c, pivot
    cdef uint8_t tmp

    for col in range(n):
        pivot = -1
        for r in range(row, m):
            if work[r, col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue

        if pivot != row:
            for c in range(n):
                tmp = work[row, c]
                work[row, c] = work[pivot, c]
                work[pivot, c] = tmp

        for r in range(m):
            if r != row and work[r, col] != 0:
                for c in range(col, n):
                    work[r, c] ^= work[row, c]

        pivots.append(col)
        row += 1
        if row == m:
            break

    return np.asarray(pivots, dtype=np.int32)


cdef cnp.ndarray[cnp.uint8_t, ndim=2] _gf2_inv(cnp.ndarray[cnp.uint8_t, ndim=2] A_arr):
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] A_work = np.ascontiguousarray(A_arr.copy() % 2, dtype=np.uint8)
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] I_arr = np.ascontiguousarray(np.eye(A_work.shape[0], dtype=np.uint8))
    cdef uint8_t[:, ::1] A = A_work
    cdef uint8_t[:, ::1] I = I_arr
    cdef int n = A_work.shape[0]
    cdef int i, j, k, pivot
    cdef uint8_t tmp

    for i in range(n):
        pivot = i
        if A[i, i] == 0:
            pivot = -1
            for j in range(i + 1, n):
                if A[j, i] != 0:
                    pivot = j
                    break
            if pivot < 0:
                raise np.linalg.LinAlgError("matrix is singular over GF(2)")

        if pivot != i:
            for k in range(n):
                tmp = A[i, k]
                A[i, k] = A[pivot, k]
                A[pivot, k] = tmp
                tmp = I[i, k]
                I[i, k] = I[pivot, k]
                I[pivot, k] = tmp

        for j in range(n):
            if j != i and A[j, i] != 0:
                for k in range(n):
                    A[j, k] ^= A[i, k]
                    I[j, k] ^= I[i, k]

    return I_arr


cpdef cnp.ndarray[cnp.uint8_t, ndim=1] ordered_statistics_postprocess(
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] llr_arr,
    cnp.ndarray[cnp.uint8_t, ndim=1, mode="c"] syndrome_arr,
    cnp.ndarray[cnp.uint8_t, ndim=2, mode="c"] H_arr,
    int osd_order,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] priors_arr,
    osd_method="osd_cs",
):
    cdef int m = H_arr.shape[0]
    cdef int n = H_arr.shape[1]

    assert osd_method == "osd_cs", "Only the benchmark osd_cs search is supported"
    assert osd_order >= 0, "osd_order must be non-negative"
    assert llr_arr.shape[0] == n, "llr has length inconsistent with H"
    assert syndrome_arr.shape[0] == m, "syndrome has length inconsistent with H"
    assert priors_arr.shape[0] == n, "priors have length inconsistent with H"
    assert np.isfinite(llr_arr).all() and np.isfinite(priors_arr).all()
    assert np.all((priors_arr > 0) & (priors_arr < 1))
    assert np.all(H_arr <= 1) and np.all(syndrome_arr <= 1)
    # The benchmark sorts signed BP values. Experimental rank scores can be
    # supplied explicitly, but channel priors always score the candidates.
    cdef cnp.ndarray[cnp.int64_t, ndim=1] reliability_order = np.ascontiguousarray(np.argsort(llr_arr, kind="stable"), dtype=np.int64)
    cdef cnp.ndarray[cnp.int64_t, ndim=1] inverse_order = np.ascontiguousarray(np.argsort(reliability_order), dtype=np.int64)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] costs = np.empty(n, dtype=np.float64)
    cdef int original_index
    for original_index in range(n):
        costs[original_index] = log(1.0 / priors_arr[original_index])
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] H_ordered = np.ascontiguousarray(H_arr[:, reliability_order], dtype=np.uint8)
    
    cdef cnp.ndarray[cnp.int32_t, ndim=1] J = _select_info_set(H_ordered)
    cdef int rank = J.shape[0]

    # Redundant checks do not change the coset. Solve on independent rows,
    # then verify the returned solution against every original constraint.
    cdef cnp.ndarray[cnp.int32_t, ndim=1] independent_rows = _select_info_set(np.ascontiguousarray(H_ordered.T))

    cdef cnp.ndarray[cnp.uint8_t, ndim=1] is_pivot_arr = np.zeros(n, dtype=np.uint8)
    cdef uint8_t[::1] is_pivot = is_pivot_arr
    cdef int idx, pos, w
    for idx in range(rank):
        is_pivot[J[idx]] = 1

    cdef cnp.ndarray[cnp.int32_t, ndim=1] J_bar = np.empty(n - rank, dtype=np.int32)
    pos = 0
    for idx in range(n):
        if is_pivot[idx] == 0:
            J_bar[pos] = idx
            pos += 1
    
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] Lambda_J = np.ascontiguousarray(H_ordered[np.ix_(independent_rows, J)], dtype=np.uint8)
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] Lambda_J_bar = np.ascontiguousarray(H_ordered[np.ix_(independent_rows, J_bar)], dtype=np.uint8)
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] Lambda_J_inv = _gf2_inv(Lambda_J)

    cdef cnp.ndarray[cnp.uint8_t, ndim=1] e0 = (Lambda_J_inv @ syndrome_arr[independent_rows]) % 2
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] best_ordered = np.zeros(n, dtype=np.uint8)
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] trial_ordered
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] e_J_bar
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] e_J
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] rhs
    cdef double best_cost = 0.0
    cdef double candidate_cost
    cdef tuple flip_indices
    cdef list patterns = []

    for idx in range(rank):
        best_ordered[J[idx]] = e0[idx]
    # Match the C++ score, including summation in original column order.
    for original_index in range(n):
        if best_ordered[inverse_order[original_index]]:
            best_cost += costs[original_index]

    if osd_order > 0:
        patterns = [(idx,) for idx in range(n - rank)]
        patterns.extend(combinations(range(min(osd_order, n - rank)), 2))
    for flip_indices in patterns:
        e_J_bar = np.zeros(n - rank, dtype=np.uint8)
        for idx in flip_indices:
            e_J_bar[idx] = 1
        rhs = (Lambda_J_bar @ e_J_bar) % 2
        e_J = (e0 + ((Lambda_J_inv @ rhs) % 2)) % 2
        trial_ordered = np.zeros(n, dtype=np.uint8)
        for idx in range(rank):
            trial_ordered[J[idx]] = e_J[idx]
        for idx in range(n - rank):
            trial_ordered[J_bar[idx]] = e_J_bar[idx]
        candidate_cost = 0.0
        for original_index in range(n):
            if trial_ordered[inverse_order[original_index]]:
                candidate_cost += costs[original_index]
        if candidate_cost < best_cost:
            best_ordered = trial_ordered
            best_cost = candidate_cost

    cdef cnp.ndarray[cnp.uint8_t, ndim=1] correction = np.empty(n, dtype=np.uint8)
    for idx in range(n):
        correction[reliability_order[idx]] = best_ordered[idx]
    if not np.array_equal(H_arr @ correction % 2, syndrome_arr):
        raise ValueError("inconsistent syndrome for the parity-check matrix")
    return correction
