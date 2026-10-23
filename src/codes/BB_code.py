import numpy as np

from lightstim.ir.qec_system import QECSystem
from lightstim.ir.tracker import SyndromeTracker
from lightstim.ir.builder import CircuitBuilder

from lightstim.qec_code.BB_code import (
    BBCode,
    BBCodeExtractionBlock,
)

from lightstim.noise.config import NoiseConfig


def build_bb_memory_circuit(
    l=12,
    m=6,
    d=12,
    A=[[3,0],[0,1],[0,2]],
    B=[[0,3],[1,0],[2,0]],
    rounds=12,
    p=1e-3,
):
    """
    Examples:
    ---------
    # [[144,12,12]] gross code
    >>> code = BBCode(l=12, m=6, A=[[3,0],[0,1],[0,2]], B=[[0,3],[1,0],[2,0]])

    # [[72,12,6]] code
    >>> code = BBCode(l=6, m=6, A=[[3,0],[0,1],[0,2]], B=[[0,3],[1,0],[2,0]], d=6)
    """
    
    code = BBCode(l=l, m=m, A=A, B=B, d=d)

    system = QECSystem()
    system.add_patch(code, name="bb")

    tracker = SyndromeTracker(system.num_qubits, expected_num_logicals=system.num_logicals)

    builder = CircuitBuilder(tracker, system)
    extraction = BBCodeExtractionBlock(system)

    builder.write_coordinates()

    builder.initialize(
        {q: "Z" for q in system.data_indices},
        n=system.num_qubits,
    )
    builder.apply_syndrome_extraction(extraction.circuit, rounds=rounds)
    builder.apply_data_readout({q: "Z" for q in system.data_indices})

    noise = NoiseConfig(
        p_1q=p,
        p_2q=p,
        p_meas=p,
        p_reset=p,
    )

    noisy_circuit = builder.build_noisy_circuit(noise, noise_model="circuit_level")

    return code, noisy_circuit, 


if __name__ == "__main__":

    code, circuit = build_bb_memory_circuit(rounds=3)
    dem = circuit.detector_error_model(decompose_errors=False)

    print(code.get_info())
    print("Circuit qubits:", circuit.num_qubits)
    print("DEM detectors:", dem.num_detectors)
    print("DEM observables:", dem.num_observables)

    # print(dem)