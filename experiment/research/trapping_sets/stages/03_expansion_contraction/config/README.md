# Configuration for 03 expansion contraction

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
| `inputs.failure_cases` | `experiment/data/input_snapshots/subgraph_structure/root/failure_cases` |
| `inputs.noise_model_d_rounds_perfect_meas_shots_expansion_contraction_max_shots` | `experiment/data/input_snapshots/subgraph_structure/root/failure_cases/{noise_model}-{d_rounds},perfect meas,shots={expansion_contraction_max_shots}` |
| `inputs.bfts_cases_csv` | `experiment/data/input_snapshots/subgraph_structure/root/failure_cases/{noise_model}-{d_rounds},perfect meas,shots={expansion_contraction_max_shots}/bfts_cases.csv` |
| `outputs.run` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main` |
| `outputs.stats` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/stats` |
| `outputs.trace` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace` |
| `outputs.config` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/config` |
| `outputs.plots` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main` |
| `outputs.cell_0015_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0015_1.html` |
| `outputs.cell_0025_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0025_1.html` |
| `outputs.trapping_set` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace/trapping_set` |
| `outputs.cell_0050_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0050_1.html` |
| `outputs.cell_0051_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0051_1.html` |
| `outputs.cell_0052_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0052_1.html` |
| `outputs.cell_0053_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0053_1.html` |
| `outputs.cell_0062_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0062_1.html` |
| `outputs.cell_0063_1_html` | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main/cell_0063_1.html` |
| `outputs.ts_depth=2_w=4_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace/trapping_set/ts_depth=2_w=4.csv` |
| `outputs.ts_depth=2_w=3_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace/trapping_set/ts_depth=2_w=3.csv` |
| `outputs.ts_depth=1_w=4_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace/trapping_set/ts_depth=1_w=4.csv` |
| `outputs.ts_depth=2_w=2_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace/trapping_set/ts_depth=2_w=2.csv` |
| `outputs.ts_depth=2_w=2_summary_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/stats/ts_depth=2_w=2_summary.csv` |
| `outputs.ts_depth=2_w=3_summary_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/stats/ts_depth=2_w=3_summary.csv` |
| `outputs.ts_depth=2_w=4_summary_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/stats/ts_depth=2_w=4_summary.csv` |
| `templates.ts_depth_depth__w_max_error_weight_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace/trapping_set/ts_depth={depth}_w={max_error_weight}.csv` |
| `templates.parameter_1__summary_csv` | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/stats/{parameter_1}_summary.csv` |

Placeholder arguments: `d_rounds`, `depth`, `expansion_contraction_max_shots`, `max_error_weight`, `noise_model`, `parameter_1`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## runtime.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
