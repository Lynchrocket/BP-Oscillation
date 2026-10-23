# Benchmark configuration and path reference

Scientific parameters are assigned in [decoder_benchmark.ipynb](../notebooks/decoder_benchmark.ipynb), not parsed from a scientific-settings JSON. The notebook's main, initial-analysis and final-figure initialization cells all read the same two local files.

## Runtime settings

[runtime.json](runtime.json) contains one key:

| Key | Default | Allowed values and effect |
| --- | --- | --- |
| `print_level` | `2` | `1`: suppress configured status/progress output; `2`: key status output; `3`: detailed output, including Sinter progress in collection calls that test this value. |

`print_level` does not change circuits, decoder parameters, sampling budgets, CSV contents or file selection. Scientific displays remain notebook outputs. Collection/plot helpers must be initialized before their cells are run; changing logging does not create a new scientific run.

## Explicit paths

[paths.json](paths.json) supplies dataset and figure paths. The complete local schema and current values follow.

All path manifests contain the following required fields; no field has an implicit fallback:

| Field | Current/default value | Meaning and constraints |
| --- | --- | --- |
| `schema_version` | `1` | Integer schema version; only 1 is accepted. |
| `inputs` | Mapping listed below | Named existing inputs; `{}` when no dataset is read. |
| `outputs` | Mapping listed below | Explicit output files/directories; includes `run`, `stats`, `trace`, `config` and `plots`. |
| `templates` | Mapping listed below | Parameterized output filenames; `{}` when unused. |

Every mapping key and value must be a non-empty string. Values are paths relative to `Decoders/`; absolute paths, parent traversal and paths escaping the project are rejected. Named inputs must exist when read. Outputs are resolved explicitly; these path fields do not create scientific settings or enable sampling. `run` is the default main dataset root, `stats` is for aggregate results, `trace` for per-case/iteration data, `config` for effective settings, and `plots` for local figures. Supplying another run name uses the path helper's run-specific directory; the notebooks shown here use their stated main paths. Do not run concurrent writers into the same result file.
### `inputs`

Empty for this manifest.

### `outputs`

| Key | Repository-relative value |
| --- | --- |
| `run` | `experiment/data/benchmark/decoder_benchmark/root/main` |
| `stats` | `experiment/data/benchmark/decoder_benchmark/root/main/stats` |
| `trace` | `experiment/data/benchmark/decoder_benchmark/root/main/trace` |
| `config` | `experiment/data/benchmark/decoder_benchmark/root/main/config` |
| `plots` | `experiment/non_exploratory/decoder_benchmarks/plots/main` |
| `benchmark` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark` |
| `2026_06_22_circuit_level_3_shots_1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-3,shots=1000000000.csv` |
| `benchmark_code_capacity_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/benchmark_code_capacity.csv` |
| `2026-07-16,circuit_level-3,force OSD,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-16,circuit_level-3,force OSD,shots=1000000000.csv` |
| `2026-08-30,circuit_level-3,adaptive BP,shots=50000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-30,circuit_level-3,adaptive BP,shots=50000000.csv` |
| `2026-08-07,circuit_level-3,quotient split,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-07,circuit_level-3,quotient split,shots=1000000000.csv` |
| `2026-07-23,circuit_level-3,alpha=0_8,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-23,circuit_level-3,alpha=0.8,shots=1000000000.csv` |
| `2026-06-22,circuit_level-1,perfect meas idle,random,serial,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-1,perfect meas idle,random,serial,shots=1000000000.csv` |
| `2026-08-29,circuit_level-3,adaptive BP,MBBP-LD,shots=50000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-29,circuit_level-3,adaptive BP,MBBP-LD,shots=50000000.csv` |
| `2026-05-28,circuit_level-1,perfect meas,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-05-28,circuit_level-1,perfect meas,shots=1000000000.csv` |
| `2026-07-23,circuit_level-1,perfect meas,alpha=1_2,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-23,circuit_level-1,perfect meas,alpha=1.2,shots=1000000000.csv` |
| `2026-06-22,circuit_level-1,perfect meas reset idle,serial,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-1,perfect meas reset idle,serial,shots=1000000000.csv` |
| `2026-07-16,circuit_level-1,force OSD,perfect meas,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-16,circuit_level-1,force OSD,perfect meas,shots=1000000000.csv` |
| `2026-06-17,code_capacity-1,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-17,code_capacity-1,shots=1000000000.csv` |
| `2026-07-23,circuit_level-1,perfect meas,alpha=0_8,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-23,circuit_level-1,perfect meas,alpha=0.8,shots=1000000000.csv` |
| `2026-06-16,circuit_level-1,perfect meas idle,random,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-16,circuit_level-1,perfect meas idle,random,shots=1000000000.csv` |
| `2026-06-17,phenomenological-2,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-17,phenomenological-2,shots=1000000000.csv` |
| `2026-06-22,circuit_level-1,perfect meas reset idle,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-1,perfect meas reset idle,shots=1000000000.csv` |
| `2026-09-09,circuit_level-3,force OSD,d=9,shots=50000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-09-09,circuit_level-3,force OSD,d=9,shots=50000000.csv` |
| `2026-06-22,circuit_level-1,perfect meas idle,serial,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-1,perfect meas idle,serial,shots=1000000000.csv` |
| `2026-06-22,circuit_level-1,perfect meas,serial,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-1,perfect meas,serial,shots=1000000000.csv` |
| `2026-07-19,circuit_level-3,adaptive OSD,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-19,circuit_level-3,adaptive OSD,shots=1000000000.csv` |
| `2026-09-01,circuit_level-3,expanded adaptive BP,shots=50000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-09-01,circuit_level-3,expanded adaptive BP,shots=50000000.csv` |
| `2026-06-17,phenomenological-1,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-17,phenomenological-1,shots=1000000000.csv` |
| `2026-06-22,circuit_level-1,perfect meas reset idle,random,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,circuit_level-1,perfect meas reset idle,random,shots=1000000000.csv` |
| `2026-06-15,circuit_level-1,perfect meas idle,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-15,circuit_level-1,perfect meas idle,shots=1000000000.csv` |
| `2026-08-24,circuit_level-3,adaptive BP,shots=1000000000_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-24,circuit_level-3,adaptive BP,shots=1000000000.csv` |

### `templates`

| Key | Repository-relative value |
| --- | --- |
| `2026_05_28_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-05-28,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_15_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_idle_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-15,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas idle,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_reset_idle_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas reset idle,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_17_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-17,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_16_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_idle_random_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-16,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas idle,random,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_reset_idle_random_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas reset idle,random,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_serial_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas,serial,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_idle_serial_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas idle,serial,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_idle_random_serial_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas idle,random,serial,shots={decoder_benchmark_max_shots}.csv` |
| `2026_06_22_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_reset_idle_serial_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-06-22,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas reset idle,serial,shots={decoder_benchmark_max_shots}.csv` |
| `2026_07_23_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_perfect_meas_alpha_alpha_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-23,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},perfect meas,alpha={alpha},shots={decoder_benchmark_max_shots}.csv` |
| `2026_07_23_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_alpha_alpha_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-23,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},alpha={alpha},shots={decoder_benchmark_max_shots}.csv` |
| `2026_07_16_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_force_OSD_perfect_meas_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-16,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},force OSD,perfect meas,shots={decoder_benchmark_max_shots}.csv` |
| `2026_07_16_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_force_OSD_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-16,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},force OSD,shots={decoder_benchmark_max_shots}.csv` |
| `2026_09_09_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_force_OSD_d_9_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-09-09,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},force OSD,d=9,shots={decoder_benchmark_max_shots}.csv` |
| `2026_07_19_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_adaptive_OSD_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-07-19,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},adaptive OSD,shots={decoder_benchmark_max_shots}.csv` |
| `2026_08_24_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_adaptive_BP_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-24,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},adaptive BP,shots={decoder_benchmark_max_shots}.csv` |
| `2026_08_29_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_adaptive_BP_MBBP_LD_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-29,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},adaptive BP,MBBP-LD,shots={decoder_benchmark_max_shots}.csv` |
| `2026_08_30_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_adaptive_BP_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-30,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},adaptive BP,shots={decoder_benchmark_max_shots}.csv` |
| `2026_09_01_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_expanded_adaptive_BP_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-09-01,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},expanded adaptive BP,shots={decoder_benchmark_max_shots}.csv` |
| `2026_08_07_decoder_benchmark_noise_model_decoder_benchmark_d_rounds_quotient_split_shots_decoder_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-08-07,{decoder_benchmark_noise_model}-{decoder_benchmark_d_rounds},quotient split,shots={decoder_benchmark_max_shots}.csv` |
| `2026_05_28_initial_benchmark_noise_model_initial_benchmark_d_rounds_perfect_meas_shots_initial_benchmark_max_shots_csv` | `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/2026-05-28,{initial_benchmark_noise_model}-{initial_benchmark_d_rounds},perfect meas,shots={initial_benchmark_max_shots}.csv` |

## Template arguments and dependencies

A template is selected by its exact key through `data_paths.output(key, **arguments)`. The argument names must match the placeholders in that template exactly; unknown/missing arguments assert. These identifiers retain historical section prefixes for filename compatibility and are not extra configuration fields.

| Placeholder | Value supplied by the notebook | Meaning |
| --- | --- | --- |
| `decoder_benchmark_noise_model` | Section's `noise_model` | `circuit_level`, `phenomenological` or `code_capacity`. |
| `decoder_benchmark_d_rounds` | Section's `d_rounds` | Positive circuit round count; shown cases use 1, 2 or 3. |
| `decoder_benchmark_max_shots` | Section's `max_shots` | Sampling cap encoded in the filename; typically 1,000,000,000 or 50,000,000 as explicitly assigned. |
| `alpha` | Section's `alpha` | Minimum-sum scaling in the selected alpha experiment; 0.8 or 1.2. Present only in alpha templates. |
| `initial_benchmark_noise_model` | `initial_benchmark_noise_model='circuit_level'` | Baseline analysis noise label. |
| `initial_benchmark_d_rounds` | `initial_benchmark_d_rounds=1` | Baseline analysis circuit round count. |
| `initial_benchmark_max_shots` | `initial_benchmark_max_shots=1000000000` | Baseline filename budget label. |

The initial-baseline template resolves to the same CSV as the main one-round perfect-measurement baseline at matching values. The initial section reads that file and does not collect another population. Static CSV keys select the specifically dated experiments listed above. Plot keys select local destinations. JSON path resolution checks path/schema/template consistency; it does not validate that a historical filename describes every currently edited decoder setting. Keep settings, task metadata and the chosen resume file consistent.
