import numpy as np



# def _build_T(
#     rules: list[int],
#     periodic: bool = False,
# ) -> np.ndarray:
#     """Build the transition matrix T for one step of CA evolution."""
#     rule_offsets = {
#         90: (-1, 1),
#         150: (-1, 0, 1),
#         102: (0, 1),
#     }
#     n = len(rules)

#     T = np.zeros((n, n), dtype=np.uint8)
#     for row in range(n):
#         offset = rule_offsets.get(rules[row], (0,))
#         for o in offset:
#             col = row + o
#             if periodic:
#                 col %= n
#                 T[row, col] ^= 1
#             else:
#                 if 0 <= col < n:
#                     T[row, col] ^= 1
#     # print(f'T:\n{T}')
    
#     return T

# def _ca_parity_matrix(
#     rules: list[int], 
#     k: int,
#     steps: int,
#     periodic: bool = False,
# ) -> np.ndarray:
#     """Build the classical parity-check matrix induced by a local CA rule."""
#     n = len(rules)
#     r = n - k

#     T = _build_T(rules, periodic)
#     # compute T^steps
#     Tt = np.eye(n, dtype=np.uint8)
#     for _ in range(steps):
#         Tt = (Tt @ T) % 2

#     print(f'Tt:\n{Tt}')

#     # G = [I_k | P]  where P = Tt[:k, k:].T  (shape k×r after transposing)
#     # encode e_i -> Tt @ e_i, take last r=n-k components as check bits.
#     G_rows = []
#     for i in range(k):
#         e = np.zeros(n, dtype=np.uint8)
#         e[i] = 1
#         codeword = Tt @ e % 2
#         G_rows.append(codeword)
#     G = np.array(G_rows)          # shape (k, n)
#     print(f'G:\n{G}')

#     # G = [I_k | P] ->  H = [P^T | I_r]
#     P = G[:, k:]                   # shape (k, r)
#     H = np.hstack([P.T, np.eye(r, dtype=np.uint8)]).astype(np.uint8)
#     print(f'H:\n{H}')

#     return H                       # shape (r, n)

# def _logical_ops(H: np.ndarray) -> np.ndarray:
#     ...

# class CACode:
#     """
#     cellular automaton code.

#     Parameters
#     ----------
#     n        : total number of physical qubits (= number of CA cells).
#     k        : number of information bits for the CA classical code.
#     d        : distance of the code.
#     rules    : list of n ints, each 90 or 150 or 102.
#     periodic : if circulant in constructing transition matrix.
#     """

#     def __init__(
#         self,
#         n: int,
#         k: int,
#         d: int,
#         rules: list[int],
#         periodic: bool = False,
#     ):
#         self.n = n
#         self.k = k
#         self.d = d
#         self.rules = rules
#         self.periodic = periodic

#         steps = d - 2
#         self.H = _ca_parity_matrix(rules, k, steps, periodic)
#         self.L = _logical_ops(self.H)

#     def get_matrices(self) -> tuple[np.ndarray, np.ndarray]:
#         return self.H, self.L

#     def __repr__(self) -> str:
#         return f"CACode: [[{self.n}, {self.k}, {self.d}]] | rules={self.rules}"


####################################################################################
# Use only the code below for the triangular code family.
####################################################################################

class CA_5_2_3():
    """[[5, 2, 3]] CA code"""
    def __init__(self):
        self.H = np.asarray([
            [1,0,1,0,0],
            [1,1,0,1,0],
            [0,1,0,0,1]
        ], dtype=np.uint8)
        self.L_X = np.asarray([
            [1, 0, 1, 1, 0],
            [0, 1, 0, 1, 1]
        ], dtype=np.uint8)
        self.L_Z = np.asarray([
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0]
        ], dtype=np.uint8)
        self.k, self.n = self.L_X.shape
        self.d = 3
    
    def stabilizer_coords(self):
        return [(2, 0), (4, 0), (6, 0)]
    
    def qubit_coords(self):
        return [(3, 0), (5, 0), (1, 0), (4, 1), (7, 0)]

class CA_9_2_5():
    """[[9, 2, 5]] CA code"""
    def __init__(self):
        self.H = np.asarray([
            [1,0,1,0,0,0,0,0,0],
            [1,1,0,1,0,0,0,0,0],
            [0,1,0,0,1,0,0,0,0],
            [0,0,1,0,0,1,0,0,0],
            [0,0,1,1,0,0,1,0,0],
            [0,0,0,1,1,0,0,1,0],
            [0,0,0,0,1,0,0,0,1],
        ], dtype=np.uint8)
        self.L_X = np.asarray([
            [1, 0, 1, 1, 0, 1, 0, 1, 0],
            [0, 1, 0, 1, 1, 0, 1, 0, 1]
        ], dtype=np.uint8)
        self.L_Z = np.asarray([
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0]
        ], dtype=np.uint8)
        self.k, self.n = self.L_X.shape
        self.d = 5

    def stabilizer_coords(self):
        return [(2, 0), (4, 0), (6, 0), (1, 2), (3, 2), (5, 2), (7, 2)]
    
    def qubit_coords(self):
        return [(3, 0), (5, 0), (2, 1), (4, 1), (6, 1), (0, 3), (3, 3), (5, 3), (8, 3)]

class TriangularCode():
    """
    Recursive triangular code family.

    level = 1 -> [[5,  2,  3]]
    level = 2 -> [[9,  2,  5]]
    level = 3 -> [[14, 2,  9]]
    level = 4 -> [[20, 2, 11]]
    level = 5 -> [[27, 2, 15]]
    ...
    """

    def __init__(self, level: int):
        self.level = level

        T1 = CA_5_2_3()
        self.H = T1.H
        self.L_X = T1.L_X
        self.L_Z = T1.L_Z
        self.k, self.n = self.L_X.shape

        # In T_1:
        #
        # bottom row:
        #     q0 q1
        #
        # frontier:
        #     q2 q3 q4
        #
        self.rows = [
            [0, 1],
            [2, 3, 4],
        ]

        self.frontier = self.rows[-1]

        # Grow recursively
        for _ in range(level - 1):
            self._grow_once()

        self.k, self.n = self.L_X.shape
        self.d = min(np.count_nonzero(x) for x in [self.L_X[0], self.L_X[1], self.L_X[0] ^ self.L_X[1]])

    def _grow_once(self):
        """
        Add one new triangular code layer.

        If old frontier has m qubits

            v0 v1 ... v_{m-1}

        the new frontier has m+1 qubits

            u0 u1 ... u_m

        and checks are

            v0 + u0 = 0

            v_{j-1} + v_j + u_j = 0

            v_{m-1} + u_m = 0.
        """

        old_frontier = self.frontier

        m = len(old_frontier)
        old_n = self.H.shape[1]

        # new layer with m+1 qubits
        new_frontier = list(range(old_n, old_n + m + 1))
        new_n = old_n + m + 1
        H_old = np.pad(self.H, pad_width=((0, 0), (0, m + 1)))

        # ------------------------------------------------------------
        # Construct new stabilizer checks
        # ------------------------------------------------------------

        new_checks = []

        # Left boundary:
        #
        # v0 --- u0
        #
        row = np.zeros(new_n, dtype=np.uint8)

        row[old_frontier[0]] = 1
        row[new_frontier[0]] = 1

        new_checks.append(row)

        # Interior triangular checks:
        #
        # v_{j-1} --- v_j
        #       \     /
        #         u_j
        #
        for j in range(1, m):

            row = np.zeros(new_n, dtype=np.uint8)

            row[old_frontier[j - 1]] = 1
            row[old_frontier[j]] = 1
            row[new_frontier[j]] = 1

            new_checks.append(row)

        # Right boundary
        row = np.zeros(new_n, dtype=np.uint8)

        row[old_frontier[-1]] = 1
        row[new_frontier[-1]] = 1

        new_checks.append(row)

        self.H = np.vstack([H_old, np.asarray(new_checks, dtype=np.uint8)])

        # ------------------------------------------------------------
        # Extend logical X operators
        #
        # Same CA rule:
        #
        # new[0] = old[0]
        # new[j] = old[j-1] XOR old[j]
        # new[m] = old[m-1]
        # ------------------------------------------------------------

        old_values = self.L_X[:, old_frontier]

        new_values = np.zeros((self.k, m + 1), dtype=np.uint8)

        new_values[:, 0] = old_values[:, 0]
        new_values[:, -1] = old_values[:, -1]

        if m > 1:
            new_values[:, 1:-1] = old_values[:, :-1] ^ old_values[:, 1:]

        self.L_X = np.hstack([self.L_X, new_values])
        self.L_Z = np.pad(self.L_Z, pad_width=((0, 0), (0, m + 1)))

        self.rows.append(new_frontier)
        self.frontier = new_frontier

    def check(self, verbose=True):
        """
        Basic consistency checks.
        """
        rank = np.linalg.matrix_rank(self.H)

        hx = (self.H @ self.L_X.T) % 2
        xz = (self.L_X @ self.L_Z.T) % 2

        ok = (
            np.all(hx == 0)
            and np.array_equal(xz, np.eye(self.k, dtype=np.uint8))
            and rank == self.n - self.k
        )

        if verbose:
            print(f"level       = {self.level}")
            print(f"n           = {self.n}")
            print(f"k           = {self.k}")
            print(f"d           = {self.d}")
            print(f"H shape     = {self.H.shape}")
            print(f"rank(H)     = {rank}")
            print(f"n - rank(H) = {self.n - rank}")
            print(f"H L_X^T = 0 : {np.all(hx == 0)}")
            print(f"L_X L_Z^T=I : {np.array_equal(xz, np.eye(self.k, dtype=np.uint8))}")
            print(f"valid       = {ok}")

        return ok

    def __repr__(self):
        return (f"Triangular(level={self.level}, [[{self.n}, {self.k}, {self.d}]])")


if __name__ == "__main__":
    for level in range(1, 4):
        code = TriangularCode(level)
        code.check()
        print(code.H)
        print(code.L_X)