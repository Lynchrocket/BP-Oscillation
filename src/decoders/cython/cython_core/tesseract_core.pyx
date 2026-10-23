# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import heapq
import math

import numpy as np
cimport numpy as cnp


ctypedef cnp.uint8_t uint8_t


cdef tuple _build_error_tables(cnp.ndarray[cnp.uint8_t, ndim=2, mode="c"] H_arr, double error_rate):
    cdef int m = H_arr.shape[0]
    cdef int n = H_arr.shape[1]
    cdef uint8_t[:, ::1] H = H_arr
    cdef int detector, error
    cdef int degree
    cdef double cost = -math.log(error_rate / (1.0 - error_rate))
    cdef dict d2e = {}
    cdef list edets = []
    cdef list touched
    cdef object error_cost = np.empty(n, dtype=np.float64)
    cdef object error_det_cost = np.empty(n, dtype=np.float64)

    for error in range(n):
        touched = []
        for detector in range(m):
            if H[detector, error] != 0:
                touched.append(detector)
                d2e.setdefault(detector, []).append(error)

        degree = len(touched)
        edets.append(tuple(touched))
        error_cost[error] = cost
        error_det_cost[error] = cost / max(1, degree)

    for detector, touched in d2e.items():
        touched.sort(key=lambda error_id: (error_det_cost[error_id], error_id))

    # The canonical expansion rule always chooses the first residual detector in
    # this order.  The wrapper currently uses natural detector order.
    cdef list det_order = sorted(d2e.keys())
    return d2e, edets, error_cost, error_det_cost, det_order


cdef set _initial_detectors(cnp.ndarray[cnp.uint8_t, ndim=1, mode="c"] syndrome_arr):
    cdef uint8_t[::1] syndrome = syndrome_arr
    cdef int detector
    cdef int m = syndrome_arr.shape[0]
    cdef set detectors = set()

    for detector in range(m):
        if syndrome[detector] != 0:
            detectors.add(detector)
    return detectors


cdef void _flip_error_detectors(set detectors, tuple error_detectors):
    cdef int detector

    for detector in error_detectors:
        if detector in detectors:
            detectors.remove(detector)
        else:
            detectors.add(detector)


cdef tuple _chain_errors(list parent, list node_error, int node_idx):
    cdef list errors = []
    cdef int idx = node_idx

    while idx > 0:
        errors.append(node_error[idx])
        idx = parent[idx]
    errors.reverse()
    return tuple(errors)


cdef tuple _replay_chain(
    set start_detectors,
    list parent,
    list node_error,
    list node_min_detector,
    list edets,
    dict d2e,
    int node_idx,
):
    cdef set detectors = set(start_detectors)
    cdef set selected = set()
    cdef set forbidden = set()
    cdef tuple errors = _chain_errors(parent, node_error, node_idx)
    cdef int error
    cdef int min_detector
    cdef int candidate

    for error in errors:
        _flip_error_detectors(detectors, edets[error])
        selected.add(error)

    # Precedence pruning: when a chain step selected `error` while expanding
    # `min_detector`, every earlier candidate for that same detector is
    # forbidden for descendants.  This canonicalizes equivalent paths through
    # the implicit error-set graph and is the key Tesseract tree construction.
    cdef list chain_indices = []
    cdef int idx = node_idx
    while idx > 0:
        chain_indices.append(idx)
        idx = parent[idx]
    chain_indices.reverse()

    for idx in chain_indices:
        min_detector = node_min_detector[idx]
        error = node_error[idx]
        for candidate in d2e.get(min_detector, []):
            if candidate == error:
                break
            forbidden.add(candidate)
        forbidden.add(error)

    return detectors, selected, forbidden, errors


cdef int _first_residual_detector(set detectors, list det_order):
    cdef int detector

    for detector in det_order:
        if detector in detectors:
            return detector
    return -1


cdef double _detector_cost(
    dict d2e,
    list edets,
    object error_cost,
    int detector,
    set residual_detectors,
    set forbidden_errors,
):
    cdef double best = math.inf
    cdef int error
    cdef int hit_count
    cdef int touched_detector

    for error in d2e.get(detector, []):
        if error in forbidden_errors:
            continue

        hit_count = 0
        for touched_detector in edets[error]:
            if touched_detector in residual_detectors:
                hit_count += 1

        if hit_count > 0:
            best = min(best, error_cost[error] / hit_count)

    return best


cdef double _heuristic(
    dict d2e,
    list edets,
    object error_cost,
    set residual_detectors,
    set forbidden_errors,
):
    cdef double total = 0.0
    cdef double cost
    cdef int detector

    for detector in residual_detectors:
        cost = _detector_cost(d2e, edets, error_cost, detector, residual_detectors, forbidden_errors)
        if math.isinf(cost):
            return math.inf
        total += cost
    return total


cdef cnp.ndarray[cnp.uint8_t, ndim=1] _correction_from_errors(int n, tuple errors):
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] correction = np.zeros(n, dtype=np.uint8)
    cdef int error

    for error in errors:
        correction[error] ^= 1
    return correction


cpdef tuple search_sector_core(
    cnp.ndarray[cnp.uint8_t, ndim=2, mode="c"] H_arr,
    cnp.ndarray[cnp.uint8_t, ndim=1, mode="c"] syndrome_arr,
    double error_rate,
    int beam,
    int pq_limit,
    double detector_penalty,
    bint no_revisit_dets=False,
):
    """Run one Tesseract sector search.

    Returns
    -------
    correction, cost, search_path, low_confidence
        `low_confidence` is true when beam / pq limits or graph exhaustion stop
        the search before reaching a syndrome-consistent error set.
    """
    cdef int n = H_arr.shape[1]
    cdef dict d2e
    cdef list edets
    cdef object error_cost
    cdef object error_det_cost
    cdef list det_order
    cdef set start_detectors
    cdef set detectors
    cdef set selected
    cdef set forbidden
    cdef set new_detectors
    cdef tuple errors
    cdef list pq = []
    cdef list search_path = []
    cdef list parent = [-1]
    cdef list node_error = [-1]
    cdef list node_min_detector = [-1]
    cdef list node_g = [0.0]
    cdef list visited_detectors = []
    cdef object visited = set()
    cdef int node_idx
    cdef int new_idx
    cdef int target_detector
    cdef int error
    cdef int pushed = 1
    cdef int best_num_detectors
    cdef double f
    cdef double g
    cdef double h
    cdef double g2
    cdef double h2
    cdef tuple det_key
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] zero_correction = np.zeros(n, dtype=np.uint8)

    assert beam >= 0
    assert pq_limit >= 1

    d2e, edets, error_cost, error_det_cost, det_order = _build_error_tables(H_arr, error_rate)
    start_detectors = _initial_detectors(syndrome_arr)
    h = _heuristic(d2e, edets, error_cost, start_detectors, set())
    if not math.isinf(h):
        h += detector_penalty * len(start_detectors)

    heapq.heappush(pq, (h, 0.0, 0))
    best_num_detectors = len(start_detectors)
    search_path.append(
        {
            "event": "start",
            "f": h,
            "g": 0.0,
            "h": h,
            "num_detectors": len(start_detectors),
            "detectors": tuple(sorted(start_detectors)),
            "errors": tuple(),
        }
    )

    while pq:
        f, g, node_idx = heapq.heappop(pq)
        detectors, selected, forbidden, errors = _replay_chain(
            start_detectors, parent, node_error, node_min_detector, edets, d2e, node_idx
        )

        search_path.append(
            {
                "event": "pop",
                "node": node_idx,
                "f": f,
                "g": g,
                "num_detectors": len(detectors),
                "detectors": tuple(sorted(detectors)),
                "errors": errors,
            }
        )

        if not detectors:
            search_path.append({"event": "goal", "node": node_idx, "cost": g, "errors": errors})
            return _correction_from_errors(n, errors), g, search_path, False

        if len(detectors) > best_num_detectors + beam:
            search_path.append(
                {
                    "event": "prune_beam",
                    "node": node_idx,
                    "num_detectors": len(detectors),
                    "best_num_detectors": best_num_detectors,
                }
            )
            continue

        best_num_detectors = min(best_num_detectors, len(detectors))

        det_key = tuple(sorted(detectors))
        if no_revisit_dets:
            if det_key in visited:
                search_path.append({"event": "prune_revisit", "node": node_idx, "detectors": det_key})
                continue
            visited.add(det_key)

        target_detector = _first_residual_detector(detectors, det_order)
        if target_detector < 0:
            search_path.append({"event": "dead_end", "node": node_idx, "detectors": det_key})
            continue

        search_path.append({"event": "expand", "node": node_idx, "target_detector": target_detector})

        for error in d2e.get(target_detector, []):
            if error in forbidden:
                continue

            new_detectors = set(detectors)
            _flip_error_detectors(new_detectors, edets[error])
            g2 = g + error_cost[error]

            # The candidate itself becomes part of the chain.  Forbidden errors
            # for the child include all lower-precedence candidates for the
            # target detector, so compute the child's heuristic after adding
            # that canonicalization constraint.
            new_idx = len(parent)
            parent.append(node_idx)
            node_error.append(error)
            node_min_detector.append(target_detector)
            node_g.append(g2)

            detectors, selected, forbidden, errors = _replay_chain(
                start_detectors, parent, node_error, node_min_detector, edets, d2e, new_idx
            )
            h2 = _heuristic(d2e, edets, error_cost, detectors, forbidden)

            if math.isinf(h2) and detectors:
                search_path.append(
                    {
                        "event": "prune_no_completion",
                        "node": new_idx,
                        "parent": node_idx,
                        "error": error,
                        "num_detectors": len(detectors),
                        "detectors": tuple(sorted(detectors)),
                        "errors": errors,
                    }
                )
                continue

            if not math.isinf(h2):
                h2 += detector_penalty * len(detectors)
            heapq.heappush(pq, (g2 + h2, g2, new_idx))
            pushed += 1
            search_path.append(
                {
                    "event": "push",
                    "node": new_idx,
                    "parent": node_idx,
                    "error": error,
                    "f": g2 + h2,
                    "g": g2,
                    "h": h2,
                    "num_detectors": len(detectors),
                    "detectors": tuple(sorted(detectors)),
                    "errors": errors,
                }
            )

            if pushed > pq_limit:
                search_path.append(
                    {
                        "event": "low_confidence",
                        "reason": "pq_limit",
                        "pushed": pushed,
                        "pq_limit": pq_limit,
                    }
                )
                return zero_correction, math.inf, search_path, True

    search_path.append({"event": "low_confidence", "reason": "exhausted"})
    return zero_correction, math.inf, search_path, True
