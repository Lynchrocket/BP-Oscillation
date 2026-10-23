"""
Union-Find Decoder

Almost-linear time decoding algorithm for topological codes (Nicolas et al.)
"""

import numpy as np
import scipy.sparse as sp
from typing import Literal
from collections import deque
import heapq

from ..base import Decoder


class UF(Decoder):
    def __init__(
        self, 
        code=None, 
        error_rate=0.001,
        uf_method = "union_find",
        **kwargs
    ):
        super().__init__(code=code, error_rate=error_rate, **kwargs)
        
        self.uf_method = uf_method
        self.uf = None

        if self.code is not None:
            self._init_backend()

    def _init_backend(self):
        self._parse_pcm(self.code)

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        """
        Union-Find Decoder:
        1. Initialize clusters
        2. Grow odd clusters until they collide
        3. Merge clusters
        4. Construct local correction within each cluster
        """
        if not hasattr(self, "H"):
            self._init_backend()

        m, n = self.H.shape

        synd_indices = np.where(syndrome)[0]
        if len(synd_indices) == 0:
            return np.zeros(n, dtype=int)

        uf = UnionFind(synd_indices)

        qubit_adj = {i: set() for i in range(n)}
        for check in range(m):
            qubits = np.where(self.H[check])[0]
            for q in qubits:
                qubit_adj[q].update(set(qubits) - {q})

        # Cluster growth using BFS
        queue = deque(synd_indices)
        visited = set(synd_indices)
        while queue:
            node = queue.popleft()
            cluster = uf.get_cluster(node)
            # Skip even clusters
            if cluster.parity == 0:
                continue
            # Grow to neighbors
            for neighbor in qubit_adj[node]:
                if neighbor not in visited:
                    queue.append(neighbor)
                    visited.add(neighbor)
                # Merge if neighbor belongs to another cluster
                if neighbor in synd_indices:
                    uf.union(node, neighbor)

        # Peeling using MST within cluster
        correction = np.zeros(n, dtype=int)
        for cluster in uf.clusters.values():
            members = list(cluster.members)
            if len(members) <= 1:
                # Single syndrome: just flip one qubit
                for q in members:
                    correction[q] ^= 1
                continue

            # Compute all pair shortest paths using BFS
            dist_map = {}
            path_map = {}
            for u in members:
                d, p = self._bfs_all_paths(u, qubit_adj)
                for v in members:
                    if v != u:
                        dist_map[(u, v)] = d[v]
                        path_map[(u, v)] = p[v]

            # Build MST using Prim's algorithm
            in_tree = {members[0]}
            edge_candidates = []
            for v in members[1:]:
                heapq.heappush(edge_candidates, (dist_map[(members[0], v)], members[0], v))

            while len(in_tree) < len(members):
                while True:
                    dist, u, v = heapq.heappop(edge_candidates)
                    if v not in in_tree:
                        break
                # Add correction along shortest path
                path = path_map[(u, v)]
                for q in path:
                    correction[q] ^= 1
                in_tree.add(v)
                # Add new edges from v to remaining members
                for w in members:
                    if w not in in_tree:
                        heapq.heappush(edge_candidates, (dist_map[(v, w)], v, w))

        return correction
    
    def _bfs_all_paths(self, start, adj):
        """
        Compute shortest distance and path from start to all reachable nodes
        """
        dist = {start: 0}
        path = {start: [start]}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for neighbor in adj[node]:
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    path[neighbor] = path[node] + [neighbor]
                    queue.append(neighbor)
        return dist, path


class Cluster:
    """
    manage syndrome points, size, and parity
    """
    def __init__(self, root):
        self.root = root        # root node in the union-find structure
        self.members = {root}   # set of syndrome points in the cluster
        self.parity = 1         # parity of the cluster (1 if odd number of syndrome points, 0 otherwise)

    def merge(self, other):
        self.members |= other.members
        self.parity = (self.parity + other.parity) % 2


class UnionFind:
    def __init__(self, nodes):
        self.parent = {node: node for node in nodes}
        self.size = {node: 1 for node in nodes}
        self.clusters = {node: Cluster(node) for node in nodes}

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return rx
        # union by size
        if self.size[rx] < self.size[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        self.size[rx] += self.size[ry]
        self.clusters[rx].merge(self.clusters[ry])
        del self.clusters[ry]
        return rx

    def get_cluster(self, x):
        return self.clusters[self.find(x)]
