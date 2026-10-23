# Implementation and model checks

These notebooks answer different implementation questions and run independently. Use conda `qec`, set the kernel working directory to the `Decoders/` repository root, and execute the chosen notebook from its initialization cell. The [configuration reference](config/README.md) explains every path-manifest field and lists its current values.

## Adaptive BP cases

[adaptive.ipynb](notebooks/adaptive.ipynb) has two parts:

1. Build a one-round circuit-level detector error model for the fixed 9-bit, distance-5 triangular code at `p=0.0002848035868435`. Decode the two 14-bit syndromes `11011000000110` and `11010000000010`. Their Adaptive BP settings are `max_iter=60`, `bp_method='minimum_sum'`, and `adaptive_k=8` or `6`, respectively. Inspect correction indices, primary/auxiliary LLRs and both convergence flags.
2. Read the shared weight-2 table and rerun Adaptive BP on each saved syndrome using that row's p. This part uses `p_meas=0`, one circuit round, `max_iter=60`, `adaptive_k=8`, and compares the predicted logical observables with the saved truth. It displays the failed/total count and the failed rows; the earlier filter also identifies BPOSD-fail/Tesseract-success examples.

The input is `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv`. Required fields are `p`, `err`, `det`, `corr_b`, `corr_t`, `pred_b`, `pred_t` and `true`; bit-string columns are read as strings to preserve leading zeros. The scan evaluates the saved fixed-weight population, not a newly sampled physical-noise LER. This notebook currently displays its results in memory and does not write an aggregate CSV or figure file. Its manifest reserves separate shared stats/trace/config and local plot directories for explicit future outputs.

## Triangular-code circuits and matching graphs

[matching_graph.ipynb](notebooks/matching_graph.ipynb) checks how the fixed 5-bit/distance-3 and 9-bit/distance-5 triangular codes appear as circuits and graphs. The notebook uses the existing fixed-code classes `CA_5_2_3` and `CA_9_2_5` directly. These classes describe triangular codes; their exported names and the module paths `codes.ca_code` and `noises_ca` retain their existing basenames for API compatibility. The recursive `TriangularCode` class is not substituted because it does not provide the coordinate methods used here.

The notebook contains:

- A local `build_circuit(code, error_rate=0.001, rounds=3)` example: data X errors, noisy CNOTs, repeated ancilla measurements, final data measurements, detectors and logical observables. The illustrated 5-bit circuit uses one and three rounds; the 9-bit circuit uses three rounds.
- `CircuitBuilder` examples for two-round phenomenological noise and one-round circuit-level noise at `p=0.01`. Circuit-level examples set `p_meas=0`; the 5-bit example also sets `p_idle=0`. A final 9-bit example uses CNOT `schedule=1`.
- Stim timeline and detector matching-graph diagrams, a Tanner graph built from the DEM, and matrix displays comparing the DEM check matrix with `code.H`. A final Tanner graph uses the code's own PCM and coordinates.

There are no input data files and no decoder sampling loop. Selected timeline/matching-graph SVGs are written under `plots/matching_graph/mg_test/`; other graph and matrix views are inline. The historical shared-data producer path `experiment/data/graph_analysis/ca_matching_graphs/root/main/` remains in the manifest for compatibility; this visualization notebook does not generate stats/trace tables there. A detector matching graph and a Tanner graph are separate representations and should be interpreted separately.

## Noise-model reference comparison

[check_noise_model.ipynb](notebooks/check_noise_model.ipynb) reads two saved reference tables for BPOSD and Tesseract at distance 5 and three rounds. Each table supplies `p_noise`, `p_best`, `p_low` and `p_high`. The p grids are taken separately from the corresponding reference table.

For each p, the notebook builds a three-round circuit-level model with `ErrorRate(p=p)`; unspecified measurement, CNOT, preparation and idle error probabilities inherit p. It constructs BPOSD with `max_iter=60`, minimum-sum BP, `osd_order=10` and `osd_method='osd_cs'`; Tesseract uses `beam=8`, `pq_limit=200000` and the existing `det_penalty=0.0` constructor argument. The sampling budget is at most `1,000,000,000` shots or `2,000` logical errors per task, using `multiprocessing.cpu_count()` workers and the CSV resume path in the manifest.

The resulting CSV is `experiment/data/noise_validation/noise_model_sanity_check/root/main/stats/bposd_tesseract_d=5_r=3.csv`. The parser groups BPOSD/Tesseract results and calls `sinter.fit_binomial(..., max_likelihood_factor=1000.0)` to obtain point estimates and likelihood intervals. The two inline log–log plots compare these intervals with the saved reference intervals; they are not declared 95% confidence intervals. This is a consistency comparison, not a formal equivalence test between noise models. The historical `MPP == Z meas ?` note is not itself a completed equivalence experiment.

## Reading results and configuring execution

The three notebooks use separate path manifests and output roots. Configuration loading checks the schema and repository-relative paths; missing named inputs raise errors. The shared manifest schema and exact input/output mappings are documented in [config/README.md](config/README.md). Scientific settings remain explicit in code cells rather than being selected through hidden defaults.

The check notebooks do not currently read a `print_level` configuration: Adaptive BP prints its final failed/total count, matching-graph cells display diagrams, and the noise comparison enables Sinter progress directly. There is no logging setting to change in their path manifests. Running the noise collection cell can launch a large simulation; read the saved CSV and execute the parser/plot cells to analyze existing results without collection.

## Available evidence and its limits

The saved noise-model check CSV has uneven sampling across settings and includes settings without observed errors. Its fitted intervals can be compared with the reference curves, but the saved observations do not establish formal noise-model equivalence or a general decoder ranking.

Matching-graph examples illustrate circuit and graph construction; they are not decoding-performance measurements. Adaptive-BP examples and fixed-support scans expose individual behavior, but no saved aggregate benchmark here establishes a whole-population advantage. Aggregate Adaptive-BP performance evaluation: 该部分实验未进行，可参考data_packedup自行分析
