# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import numpy as np
cimport numpy as cnp
from libc.math cimport fabs, tanh, atanh
from time import perf_counter

ctypedef cnp.uint8_t uint8_t
ctypedef cnp.int32_t int32_t


cdef inline double _sign(double x) noexcept:
    if x >= 0.0:
        return 1.0
    else:
        return -1.0


cdef bint _matches_syndrome(
    int32_t[::1] head,
    int32_t[::1] head_nxt,
    int32_t[::1] to,
    uint8_t[::1] syndrome,
    uint8_t[::1] correction,
    int m,
    int n,
) noexcept:
    cdef int i, j, ce, var
    cdef unsigned char acc

    for i in range(m):
        acc = 0
        ce = head[i]
        while ce != -1:
            var = to[ce]
            acc ^= correction[var - m]
            ce = head_nxt[ce]
        if acc != syndrome[i]:
            return False
    return True


cdef inline void _llr_from_messages(
    double[::1] llr,
    double[::1] msg,
    int32_t[::1] head,
    int32_t[::1] head_nxt,
    double[::1] channel_llr,
    int m,
    int n,
) noexcept:
    cdef int j, ve
    cdef double total

    for j in range(n):
        total = channel_llr[j] if channel_llr.shape[0] > 1 else channel_llr[0]
        ve = head[j + m]
        while ve != -1:
            total += msg[ve-1]
            ve = head_nxt[ve]
        llr[j] = total


cdef inline void _hard_decision(uint8_t[::1] correction, double[::1] llr, int n) noexcept:
    cdef int j
    for j in range(n):
        correction[j] = 1 if llr[j] < 0.0 else 0


cpdef tuple min_sum_trace(
    cnp.ndarray[cnp.uint8_t, ndim=1, mode="c"] syndrome_arr,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] msg_arr,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] head_arr,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] head_nxt_arr,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] to_arr,
    int m,
    int n,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] channel_llr_arr,
    int max_iter,
    double scaling_factor=1.0,
    bint need_tracking=True,
):
    cdef uint8_t[::1] syndrome = syndrome_arr
    cdef double[::1] msg = msg_arr
    cdef int32_t[::1] head = head_arr
    cdef int32_t[::1] head_nxt = head_nxt_arr
    cdef int32_t[::1] to = to_arr
    cdef double[::1] channel_llr = channel_llr_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=2] history_arr = np.empty((max_iter, n), dtype=np.float64)
    cdef double[:, ::1] history = history_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] llr_arr = np.empty(n, dtype=np.float64)
    cdef double[::1] llr = llr_arr
    cdef uint8_t[::1] correction = np.empty(n, dtype=np.uint8)
    cdef int i, j, ce, ve, iter_idx
    cdef double cllr, sign_prod, min_abs, min_abs2, abs_value, total, this_sign_prod, this_min_abs
    cdef int degree
    cdef double[::1] incoming = np.empty(msg.shape[0] // 2, dtype=np.float64)
    cdef double[::1] outgoing = np.empty(msg.shape[0] // 2, dtype=np.float64)

    if channel_llr.shape[0] == 1:
        cllr = channel_llr[0]
        for i in range(0, msg.shape[0], 2):
            msg[i] = 0.0
            msg[i+1] = cllr
    else:
        for i in range(0, msg.shape[0], 2):
            msg[i] = 0.0
        for j in range(n):
            ve = head[j + m]
            cllr = channel_llr[j]
            while ve != -1:
                msg[ve] = cllr
                ve = head_nxt[ve]

    for iter_idx in range(max_iter):
        for i in range(m):
            degree = 0
            ce = head[i]
            while ce != -1:
                incoming[degree] = msg[ce + 1]
                degree += 1
                ce = head_nxt[ce]
            _min_sum_node(incoming, outgoing, degree, syndrome[i], False)
            degree = 0
            ce = head[i]
            while ce != -1:
                msg[ce] = scaling_factor * outgoing[degree]
                degree += 1
                ce = head_nxt[ce]

        for j in range(n):
            ve = head[j + m] # 1st edge from var j
            total = channel_llr[j] if channel_llr.shape[0] > 1 else channel_llr[0]
            # sweep through all edges from var j to compute the total message to check neighbors
            while ve != -1:
                total += msg[ve-1]
                ve = head_nxt[ve] # next edge from var j
            ve = head[j + m] # reset to 1st edge from var j
            while ve != -1:
                msg[ve] = total - msg[ve-1] # message from var j to check neighbor = total message - message from check neighbor to var j
                ve = head_nxt[ve]

        _llr_from_messages(llr, msg, head, head_nxt, channel_llr, m, n)
        if need_tracking:
            for j in range(n):
                history[iter_idx, j] = llr[j]
        _hard_decision(correction, llr, n)
        if _matches_syndrome(head, head_nxt, to, syndrome, correction, m, n):
            return llr_arr, history_arr[:iter_idx + 1]

    return llr_arr, history_arr


cpdef tuple product_sum_trace(
    cnp.ndarray[cnp.uint8_t, ndim=1, mode="c"] syndrome_arr,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] msg_arr,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] head_arr,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] head_nxt_arr,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] to_arr,
    int m,
    int n,    
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] channel_llr_arr,
    int max_iter,
    bint need_tracking=True,
):
    cdef uint8_t[::1] syndrome = syndrome_arr
    cdef double[::1] msg = msg_arr
    cdef int32_t[::1] head = head_arr
    cdef int32_t[::1] head_nxt = head_nxt_arr
    cdef int32_t[::1] to = to_arr
    cdef double[::1] channel_llr = channel_llr_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=2] history_arr = np.empty((max_iter, n), dtype=np.float64)
    cdef double[:, ::1] history = history_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] llr_arr = np.empty(n, dtype=np.float64)
    cdef double[::1] llr = llr_arr
    cdef uint8_t[::1] correction = np.empty(n, dtype=np.uint8)
    cdef int i, j, ce, ve, iter_idx, zero_cnt
    cdef double cllr, syndrome_sign, tanh_prod, total, this_tanh_prod

    if channel_llr.shape[0] == 1:
        cllr = channel_llr[0]
        for i in range(0, msg.shape[0], 2):
            msg[i] = 0.0
            msg[i+1] = cllr
    else:
        for i in range(0, msg.shape[0], 2):
            msg[i] = 0.0
        for j in range(n):
            ve = head[j + m]
            cllr = channel_llr[j]
            while ve != -1:
                msg[ve] = cllr
                ve = head_nxt[ve]

    for iter_idx in range(max_iter):
        for i in range(m):
            syndrome_sign = -1.0 if syndrome[i] != 0 else 1.0
            tanh_prod = 1.0
            # sweep through all edges from check i to compute the product of signs and the minimum absolute value
            ce = head[i] # 1st edge from check i
            zero_cnt = 0
            while ce != -1:
                if msg[ce+1] == 0.0:
                    zero_cnt += 1
                else:
                    tanh_prod *= tanh(msg[ce+1] / 2.0)
                ce = head_nxt[ce] # next edge from check i
            # compute messages from check i to its neighboring variables
            ce = head[i] # reset to 1st edge from check i
            while ce != -1:
                if zero_cnt >= 2:
                    this_tanh_prod = 0.0
                elif zero_cnt == 1:
                    this_tanh_prod = tanh_prod if msg[ce+1] == 0.0 else 0.0
                else:
                    this_tanh_prod = tanh_prod / tanh(msg[ce+1] / 2.0)
                # Clip the extrinsic product actually passed to atanh.
                if this_tanh_prod > 1.0 - 1e-12:
                    this_tanh_prod = 1.0 - 1e-12
                elif this_tanh_prod < -1.0 + 1e-12:
                    this_tanh_prod = -1.0 + 1e-12
                msg[ce] = 2.0 * syndrome_sign * atanh(this_tanh_prod)
                ce = head_nxt[ce]

        for j in range(n):
            ve = head[j + m] # 1st edge from var j
            total = channel_llr[j] if channel_llr.shape[0] > 1 else channel_llr[0]
            # sweep through all edges from var j to compute the total message to check neighbors
            while ve != -1:
                total += msg[ve-1]
                ve = head_nxt[ve] # next edge from var j
            ve = head[j + m] # reset to 1st edge from var j
            while ve != -1:
                msg[ve] = total - msg[ve-1] # message from var j to check neighbor = total message - message from check neighbor to var j
                ve = head_nxt[ve]

        _llr_from_messages(llr, msg, head, head_nxt, channel_llr, m, n)
        if need_tracking:
            for j in range(n):
                history[iter_idx, j] = llr[j]
        _hard_decision(correction, llr, n)
        if _matches_syndrome(head, head_nxt, to, syndrome, correction, m, n):
            return llr_arr, history_arr[:iter_idx + 1]

    return llr_arr, history_arr


cdef void _min_sum_node(
    double[::1] incoming,
    double[::1] outgoing,
    int degree,
    unsigned char syndrome,
    bint finite_degree_one,
) noexcept:
    cdef int k
    cdef double signs = -1.0 if syndrome else 1.0
    cdef double minimum = 1.7976931348623157e308
    cdef double second = 1.7976931348623157e308
    cdef double value, magnitude
    for k in range(degree):
        signs *= _sign(incoming[k])
        value = fabs(incoming[k])
        if value < minimum:
            second = minimum
            minimum = value
        elif value < second:
            second = value
    if finite_degree_one and degree == 1:
        second = minimum
    for k in range(degree):
        magnitude = minimum if fabs(incoming[k]) > minimum else second
        outgoing[k] = signs * _sign(incoming[k]) * magnitude


cdef void _research_product_node(
    double[::1] incoming,
    double[::1] outgoing,
    double[::1] values,
    int degree,
    unsigned char syndrome,
    double epsilon,
) noexcept:
    cdef int k, j
    cdef double product, value
    for k in range(degree):
        value = tanh(incoming[k] / 2.0)
        values[k] = min(1.0 - epsilon, max(-1.0 + epsilon, value))
    # Match the research clipping and excluded-edge product, including zero messages.
    for k in range(degree):
        product = -1.0 if syndrome else 1.0
        for j in range(degree):
            if j != k:
                product *= values[j]
        product = min(1.0 - epsilon, max(-1.0 + epsilon, product))
        outgoing[k] = 2.0 * atanh(product)


cdef double _variable_total(
    int variable,
    double[::1] prior,
    double[::1] c2v,
    int32_t[::1] variable_ptr,
    int32_t[::1] variable_edges,
    int mode,
    double[::1] scratch,
):
    cdef int k, begin = variable_ptr[variable], end = variable_ptr[variable + 1]
    cdef double total = prior[variable] if mode == 0 else 0.0
    if mode == 2:
        for k in range(begin, end):
            scratch[k - begin] = c2v[variable_edges[k]]
        return prior[variable] + np.sum(np.asarray(scratch)[:end - begin])
    for k in range(begin, end):
        total += c2v[variable_edges[k]]
    return total if mode == 0 else prior[variable] + total


cpdef tuple research_trace(
    uint8_t[::1] syndrome,
    double[::1] prior,
    int32_t[::1] check_ptr,
    int32_t[::1] edge_variables,
    int32_t[::1] variable_ptr,
    int32_t[::1] variable_edges,
    int iterations,
    double[::1] scaling,
    double[::1] damping,
    double offset,
    bint product_sum,
    int32_t[:, ::1] check_orders,
    bint serial,
    double[::1] initial_v2c,
    int accumulation,
    int serial_accumulation,
    bint include_initial_valid,
    bint stop_early,
    bint record_messages,
    bint combine_scaling_damping,
):
    """Run a research trajectory in row-major edge order; all shot state is local."""
    cdef double started = perf_counter()
    cdef int m = syndrome.shape[0], n = prior.shape[0], edges = edge_variables.shape[0]
    cdef double[::1] c2v = np.zeros(edges, dtype=np.float64)
    cdef double[::1] v2c = np.array(initial_v2c, dtype=np.float64, copy=True)
    cdef double[::1] incoming = np.empty(edges, dtype=np.float64)
    cdef double[::1] outgoing = np.empty(edges, dtype=np.float64)
    cdef double[::1] values = np.empty(edges, dtype=np.float64)
    cdef double[::1] scratch = np.empty(edges, dtype=np.float64)
    cdef cnp.ndarray c_history_arr = np.empty((iterations + 1 if record_messages else 0, edges), dtype=np.float64)
    cdef cnp.ndarray v_history_arr = np.empty((iterations + 1 if record_messages else 0, edges), dtype=np.float64)
    cdef double[:, ::1] c_history = c_history_arr, v_history = v_history_arr
    cdef cnp.ndarray llr_arr = np.empty((iterations + 1, n), dtype=np.float64)
    cdef double[:, ::1] llr = llr_arr
    cdef cnp.ndarray hard_arr = np.empty((iterations + 1, n), dtype=np.uint8)
    cdef uint8_t[:, ::1] hard = hard_arr
    cdef cnp.ndarray valid_arr = np.zeros(iterations + 1, dtype=np.bool_)
    cdef cnp.ndarray runtime_arr = np.zeros(iterations + 1, dtype=np.float64)
    cdef double[::1] runtime = runtime_arr
    cdef int t, k, check, e, j, variable, degree, begin, end, first = -1
    cdef int completed = iterations, order_row
    cdef double proposed, total
    cdef unsigned char parity
    cdef bint valid

    for j in range(n):
        llr[0, j] = prior[j]
        hard[0, j] = 1 if prior[j] < 0.0 else 0
    valid = True
    for check in range(m):
        parity = 0
        for e in range(check_ptr[check], check_ptr[check + 1]):
            parity ^= hard[0, edge_variables[e]]
        if parity != syndrome[check]:
            valid = False
    valid_arr[0] = valid
    if include_initial_valid and valid:
        first = 0
    if record_messages:
        for e in range(edges):
            c_history[0, e] = 0.0
            v_history[0, e] = v2c[e]

    for t in range(1, iterations + 1):
        order_row = 0 if check_orders.shape[0] == 1 else t - 1
        for k in range(m):
            check = check_orders[order_row, k]
            begin = check_ptr[check]
            end = check_ptr[check + 1]
            degree = end - begin
            for e in range(begin, end):
                incoming[e - begin] = v2c[e]
            if product_sum:
                _research_product_node(incoming, outgoing, values, degree, syndrome[check], 1e-15)
            else:
                _min_sum_node(incoming, outgoing, degree, syndrome[check], True)
            for e in range(begin, end):
                proposed = outgoing[e - begin]
                if combine_scaling_damping and not product_sum:
                    c2v[e] = (1.0 - damping[e]) * scaling[e] * proposed + damping[e] * c2v[e]
                    continue
                if not product_sum:
                    proposed = _sign(proposed) * max(scaling[e] * fabs(proposed) - offset, 0.0)
                c2v[e] = (1.0 - damping[e]) * proposed + damping[e] * c2v[e]
            if serial:
                for e in range(begin, end):
                    variable = edge_variables[e]
                    total = _variable_total(variable, prior, c2v, variable_ptr, variable_edges, serial_accumulation, scratch)
                    for j in range(variable_ptr[variable], variable_ptr[variable + 1]):
                        v2c[variable_edges[j]] = total - c2v[variable_edges[j]]

        for variable in range(n):
            total = _variable_total(variable, prior, c2v, variable_ptr, variable_edges, accumulation, scratch)
            llr[t, variable] = total
            hard[t, variable] = 1 if total < 0.0 else 0
        if not serial:
            for e in range(edges):
                v2c[e] = llr[t, edge_variables[e]] - c2v[e]
        valid = True
        for check in range(m):
            parity = 0
            for e in range(check_ptr[check], check_ptr[check + 1]):
                parity ^= hard[t, edge_variables[e]]
            if parity != syndrome[check]:
                valid = False
        valid_arr[t] = valid
        if first == -1 and valid:
            first = t
        if record_messages:
            for e in range(edges):
                c_history[t, e] = c2v[e]
                v_history[t, e] = v2c[e]
        runtime[t] = perf_counter() - started
        if stop_early and valid:
            completed = t
            break
    end = completed + 1
    return (c_history_arr[:end], v_history_arr[:end], llr_arr[:end],
            hard_arr[:end], valid_arr[:end], runtime_arr[:end],
            None if first == -1 else first)
