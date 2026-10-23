# Configuration for 02 candidates

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [paths.json](paths.json) | Named directly by the current notebook code. |
| [runtime.json](runtime.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.failure_cases_with_errors` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors` |
| `inputs.ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bfts_cases.csv` |
| `inputs.ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bstf_cases.csv` |
| `inputs.ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bs_cases.csv` |
| `inputs.ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bf_cases.csv` |
| `outputs.run` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main` |
| `outputs.stats` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/stats` |
| `outputs.trace` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace` |
| `outputs.config` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/config` |
| `outputs.plots` | `experiment/research/trapping_sets/stages/02_candidates/plots/main` |
| `outputs.circuit_level` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace/circuit_level` |
| `outputs.candidate_error_weight_distribution_csv` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/stats/candidate_error_weight_distribution.csv` |
| `outputs.cell_0016_1_png` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0016_1.png` |
| `outputs.cell_0019_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0019_1.html` |
| `outputs.cell_0020_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0020_1.html` |
| `outputs.cell_0021_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0021_1.html` |
| `outputs.cell_0022_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0022_1.html` |
| `outputs.cell_0023_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0023_1.html` |
| `outputs.cell_0024_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0024_1.html` |
| `outputs.cell_0025_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0025_1.html` |
| `outputs.cell_0026_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0026_1.html` |
| `outputs.cell_0028_1_html` | `experiment/research/trapping_sets/stages/02_candidates/plots/main/cell_0028_1.html` |
| `outputs.TS_failure_cases_feather` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace/circuit_level/trapping_set/CA_9_2_5,p=0.1,(a_max=6,b_max=10)/TS_failure_cases.feather` |
| `outputs.TS_overview_csv` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/stats/circuit_level/trapping_set/CA_9_2_5,p=0.1,(a_max=6,b_max=10)/TS_overview.csv` |
| `templates.parameter_1_p_parameter_2_a_max_a_max_b_max_b_max` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace/circuit_level/trapping_set/{parameter_1},p={parameter_2},(a_max={a_max},b_max={b_max})` |
| `templates.TS_overview_csv` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace/circuit_level/trapping_set/{parameter_1},p={parameter_2},(a_max={a_max},b_max={b_max})/TS_overview.csv` |
| `templates.TS_failure_cases_feather` | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace/circuit_level/trapping_set/{parameter_1},p={parameter_2},(a_max={a_max},b_max={b_max})/TS_failure_cases.feather` |

Placeholder arguments: `a_max`, `b_max`, `parameter_1`, `parameter_2`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## runtime.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
