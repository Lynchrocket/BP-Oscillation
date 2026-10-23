"""
Search-based Decoder

A search based decoder based on the paper "Tesseract: A Search-Based Decoder for Quantum Error Correction".
"""

import numpy as np
import math
import heapq
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict, List, Set
import time

from ..base import Decoder


class Tesseract(Decoder):
    @dataclass
    class _SearchInfo:
        d2e: Dict[int, List[int]]                       # detector -> list of error ids (errors that trigger this detector)
        edets: Dict[int, Set[int]]                      # error -> detectors (detectors triggered by this error)
        error_cost: Dict[int, Tuple[float, float]]      # error -> tuple(heuristic cost, actual cost)
        det_order: List[int]                            # precedence order for detectors
    
    @dataclass(order=True)
    class _Node:
        f: float                                        # f = g + h
        g: float                                        # actual cost already incurred
        detectors: Set[int] = field(compare=False)      # unsatisfied detectors in this node
        errors: Set[int] = field(compare=False)         # Set of error ids selected in this node
        blocked: Set[int] = field(compare=False)        # Set of error ids that are blocked (already selected in this path or pruned)


    def __init__(
        self, 
        code=None, 
        error_rate: float = 0.001,
        beam: int = 10,
        pq_limit: int = 1000000,
        det_penalty: float = 0.0,
        **kwargs
    ):
        """
        Initialize Search-based decoder.
        
        Parameters
        -----------
        code : CSS
            The CSS code to decode
        error_rate : float
            Physical error probability
        beam : int
            Beam width for search pruning
        pq_limit : int
            Priority queue size limit
        det_penalty : float
            Additional penalty added to the heuristic cost of each unsatisfied detector to encourage faster convergence
        """
        super().__init__(code=code, error_rate=error_rate, **kwargs)

        self.beam = beam
        self.pq_limit = pq_limit
        self.detector_penalty = det_penalty

        self.search_info_x = self._init_structures(code.H_Z, code.n)
        self.search_info_z = self._init_structures(code.H_X, code.n)

    def _init_backend(self):
        self._parse_pcm(self.code)


    def _init_structures(self, H: np.ndarray, num: int, det_order: Optional[List[int]] = None) -> _SearchInfo:
        """
        Parameters
        ----------
        H : np.ndarray
            Stabilizer check matrix (H_X or H_Z)
        num : int
            Number of physical qubits (columns of H)
        det_order : Optional[List[int]]
            Optional predefined order of detectors to prioritize

        Returns
        -------
        _SearchInfo
        """
        d2e: Dict[int, List[int]] = {}
        edets: Dict[int, Set[int]] = {}
        error_cost = {}

        cost = -np.log(self.error_rate / (1 - self.error_rate))

        for id in range(num):
            dets = set(np.where(H[:, id] == 1)[0])

            edets[id] = dets

            for d in dets:
                d2e.setdefault(d, []).append(id)

            error_cost[id] = (cost / max(1, len(dets)), cost)

        # sort detector->errors
        for d, elist in d2e.items():
            elist.sort(key=lambda e_id: error_cost[e_id][0])

        if det_order is None:
            det_order = sorted(d2e.keys(), key=lambda d: min(error_cost[e_id][0] for e_id in d2e[d]))

        return self._SearchInfo(d2e=d2e, edets=edets, error_cost=error_cost, det_order=det_order)

    def _get_detcost(self, si: _SearchInfo, d: int, detectors_left: Set[int], blocked_errors: Set[int]) -> float:
        """
        Compute the heuristic cost of satisfying a single detector.

        Parameters
        ----------
        si : _SearchInfo
            Precomputed search structures
        d : int
            Detector index to compute cost for
        detectors_left : set[int]
            Remaining detectors that need to be satisfied
        blocked_errors : set[int]
            Errors that are currently blocked in this branch

        Returns
        -------
        float
            Heuristic cost for this detector (including detector penalty)
        """
        best = math.inf
        for ei in si.d2e.get(d, []):
            if ei in blocked_errors:
                continue
            dets = si.edets[ei] & detectors_left
            if not dets:
                continue
            _, cost = si.error_cost[ei]
            best = min(best, cost / len(dets))
        return best + self.detector_penalty

    def _heuristic(self, si: _SearchInfo, detectors_left: Set[int], blocked_errors: Set[int]) -> float:
        return sum(
            self._get_detcost(si, d, detectors_left, blocked_errors) for d in detectors_left
        )
    
    def _search(self, si: _SearchInfo, syndrome: np.ndarray) -> Tuple[Set[int], float]:
        """
        A* Search-based decoder for a given syndrome and search info.

        Parameters
        ----------
        si : _SearchInfo
            Precomputed structure used for the search
        syndrome : np.ndarray
            Syndrome vector for a given type of error (X or Z)

        Returns
        -------
        errors : Set[int]
            Set of error indices (columns of H) that form the correction
        total_cost : float
            Total actual cost (sum of error probabilities in log-likelihood)
        """
        start_detectors = set(np.where(syndrome == 1)[0])
       
        h0 = self._heuristic(si, start_detectors, set())
        start = self._Node(
            f=h0,
            g=0.0,
            detectors=start_detectors,
            errors=set(),
            blocked=set(),
        )

        pq = [start]
        pushed = 1
        best_num_dets = len(start_detectors)

        while pq:
            node = heapq.heappop(pq)

            # ------------------------------
            # Goal: 
            # If no detectors left, syndrome fully explained
            # Return selected errors and total cost
            # ------------------------------
            if not node.detectors:
                return set(node.errors), node.g

            # ------------------------------
            # Beam pruning
            # ------------------------------
            if len(node.detectors) > best_num_dets + self.beam:
                continue

            best_num_dets = min(best_num_dets, len(node.detectors))

            # ------------------------------
            # Select target detector to expand
            # precedence: pick smallest detector
            # ------------------------------
            for d in si.det_order:
                if d in node.detectors:
                    target_det = d
                    break
            else:
                continue # skip the detectors that are not in det_order

            # check all candidate errors for this detector
            for ei in si.d2e.get(target_det, []):
                if ei in node.blocked or ei in node.errors:
                    continue

                # apply error
                new_detectors = set(node.detectors)
                for d in si.edets[ei]:
                    if d in new_detectors:
                        new_detectors.remove(d)
                    else:
                        new_detectors.add(d)

                new_errors = node.errors.union({ei})
                new_blocked = node.blocked.union({ei})

                g2 = node.g + si.error_cost[ei][1]
                h2 = self._heuristic(si, new_detectors, new_blocked)

                new_node = self._Node(
                    f=g2 + h2,
                    g=g2,
                    detectors=new_detectors,
                    errors=new_errors,
                    blocked=new_blocked,
                )

                heapq.heappush(pq, new_node)
                pushed += 1
                if pushed > self.pq_limit:
                    return set(node.errors), node.g

        return set(), math.inf

    def decode_to_correction(self, syndrome: np.ndarray) -> np.ndarray:
        """
        Decode using SB.
        
        Parameters
        -----------
        syndrome : np.ndarray, shape (c_z + c_x,)
            
        Returns
        --------
        correction : np.ndarray, shape (2*n,)
        """
        start_time = time.time()

        syndrome_x = syndrome[:self.c_z]
        syndrome_z = syndrome[self.c_z:]

        errors_X, cost_X = self._search(
            self.search_info_x,
            syndrome_x
        )

        errors_Z, cost_Z = self._search(
            self.search_info_z,
            syndrome_z
        )

        correction_x = np.zeros(self.n, dtype=int)
        for e in errors_X:
            correction_x[e] = 1

        correction_z = np.zeros(self.n, dtype=int)
        for e in errors_Z:
            correction_z[e] = 1

        correction = np.concatenate([correction_x, correction_z])

        return correction