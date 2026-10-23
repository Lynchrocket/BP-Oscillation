# Configuration for 04 low weight structure

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
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `inputs.weight_4_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_4_data.csv` |
| `outputs.run` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main` |
| `outputs.stats` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats` |
| `outputs.trace` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace` |
| `outputs.config` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/config` |
| `outputs.plots` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main` |
| `outputs.subgraph` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph` |
| `outputs.weight_2_err_subgraph_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph/weight-2,err_subgraph.csv` |
| `outputs.weight_2_det_subgraph_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph/weight-2,det_subgraph.csv` |
| `outputs.cell_0093_1_png` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0093_1.png` |
| `outputs.cell_0098_1_png` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0098_1.png` |
| `outputs.cell_0103_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0103_1.html` |
| `outputs.cell_0104_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0104_1.html` |
| `outputs.cell_0105_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0105_1.html` |
| `outputs.cell_0108_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0108_1.html` |
| `outputs.cell_0109_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0109_1.html` |
| `outputs.cell_0110_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0110_1.html` |
| `outputs.weight_3_err_subgraph_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph/weight-3,err_subgraph.csv` |
| `outputs.weight_3_det_subgraph_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph/weight-3,det_subgraph.csv` |
| `outputs.cell_0153_1_png` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0153_1.png` |
| `outputs.cell_0158_1_png` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0158_1.png` |
| `outputs.cell_0165_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0165_1.html` |
| `outputs.cell_0166_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0166_1.html` |
| `outputs.cell_0167_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0167_1.html` |
| `outputs.weight_4_err_subgraph_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph/weight-4,err_subgraph.csv` |
| `outputs.weight_4_det_subgraph_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace/subgraph/weight-4,det_subgraph.csv` |
| `outputs.cell_0193_1_png` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0193_1.png` |
| `outputs.cell_0198_1_png` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0198_1.png` |
| `outputs.cell_0205_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0205_1.html` |
| `outputs.cell_0206_1_html` | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main/cell_0206_1.html` |
| `outputs.weight-2,error_induced_summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-2,error_induced_summary.csv` |
| `outputs.weight-2,detector_induced_summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-2,detector_induced_summary.csv` |
| `outputs.weight-3,error_induced_summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-3,error_induced_summary.csv` |
| `outputs.weight-3,detector_induced_summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-3,detector_induced_summary.csv` |
| `outputs.weight-4,error_induced_summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-4,error_induced_summary.csv` |
| `outputs.weight-4,detector_induced_summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-4,detector_induced_summary.csv` |
| `templates.weight_2_population__summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-2,{population}_summary.csv` |
| `templates.weight_3_population__summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-3,{population}_summary.csv` |
| `templates.weight_4_population__summary_csv` | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats/weight-4,{population}_summary.csv` |

Placeholder arguments: `population`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## runtime.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
