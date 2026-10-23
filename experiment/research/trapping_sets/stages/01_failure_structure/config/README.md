# Configuration for 01 failure structure

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
| `inputs.failure_cases` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases` |
| `inputs.ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bfts_cases.csv` |
| `inputs.ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bstf_cases.csv` |
| `inputs.ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bs_cases.csv` |
| `inputs.ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bf_cases.csv` |
| `outputs.run` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main` |
| `outputs.stats` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats` |
| `outputs.trace` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace` |
| `outputs.config` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/config` |
| `outputs.plots` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main` |
| `outputs.circuit_level` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level` |
| `outputs.cell_0011_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0011_1.png` |
| `outputs.cell_0013_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0013_1.png` |
| `outputs.graph_analysis` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis` |
| `outputs.event_graph_features_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/event_graph_features.csv` |
| `outputs.cell_0019_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0019_1.png` |
| `outputs.cell_0021_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0021_1.png` |
| `outputs.cell_0022_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0022_1.png` |
| `outputs.cell_0026_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0026_1.png` |
| `outputs.cell_0030_1_html` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0030_1.html` |
| `outputs.cell_0031_1_html` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0031_1.html` |
| `outputs.trapping_set` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set` |
| `outputs.bfts_trapping_param_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_trapping_param_summary.csv` |
| `outputs.bfts_error_node_parity_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_error_node_parity_summary.csv` |
| `outputs.bfts_wl_hash_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_wl_hash_summary.csv` |
| `outputs.bfts_wl_trapping_cross_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_wl_trapping_cross_summary.csv` |
| `outputs.cell_0042_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0042_1.png` |
| `outputs.bfts_classical_trapping_candidates_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_classical_trapping_candidates.csv` |
| `outputs.bfts_classical_trapping_by_wl_hash_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_classical_trapping_by_wl_hash.csv` |
| `outputs.cell_0046_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0046_1.png` |
| `outputs.bfts_symmetric_stabilizer_by_wl_hash_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_symmetric_stabilizer_by_wl_hash.csv` |
| `outputs.bfts_symmetric_stabilizer_candidates_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_symmetric_stabilizer_candidates.csv` |
| `outputs.cell_0049_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0049_1.png` |
| `outputs.bfts_wl_feature_distance_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_wl_feature_distance.csv` |
| `outputs.cell_0054_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0054_1.png` |
| `outputs.bfts_feature_distance_clusters_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_feature_distance_clusters.csv` |
| `outputs.bfts_representative_ged_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_representative_ged.csv` |
| `outputs.cell_0060_1_png` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0060_1.png` |
| `outputs.bfts_ged_similarity_clusters_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/bfts_ged_similarity_clusters.csv` |
| `outputs.adversarial_trapping_events_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/adversarial_trapping_events.csv` |
| `outputs.adversarial_trapping_decoder_results_long_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/adversarial_trapping_decoder_results_long.csv` |
| `outputs.adversarial_trapping_decoder_results_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/adversarial_trapping_decoder_results.csv` |
| `outputs.adversarial_trapping_decoder_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/adversarial_trapping_decoder_summary.csv` |
| `outputs.adversarial_trapping_csv_dtypes_json` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/trapping_set/adversarial_trapping_csv_dtypes.json` |
| `outputs.cell_0072_1_html` | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main/cell_0072_1.html` |
| `outputs.graph_analysis_event_graph_features_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/event_graph_features.csv` |
| `outputs.trapping_set_bfts_wl_feature_distance_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/bfts_wl_feature_distance.csv` |
| `outputs.trapping_set_adversarial_trapping_csv_dtypes_json` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/adversarial_trapping_csv_dtypes.json` |
| `outputs.trapping_set_bfts_representative_ged_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/bfts_representative_ged.csv` |
| `outputs.trapping_set_bfts_ged_similarity_clusters_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/bfts_ged_similarity_clusters.csv` |
| `outputs.trapping_set_bfts_classical_trapping_candidates_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/bfts_classical_trapping_candidates.csv` |
| `outputs.trapping_set_adversarial_trapping_decoder_results_long_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/adversarial_trapping_decoder_results_long.csv` |
| `outputs.trapping_set_bfts_feature_distance_clusters_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/bfts_feature_distance_clusters.csv` |
| `outputs.trapping_set_adversarial_trapping_decoder_results_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/adversarial_trapping_decoder_results.csv` |
| `outputs.trapping_set_bfts_symmetric_stabilizer_candidates_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/bfts_symmetric_stabilizer_candidates.csv` |
| `outputs.trapping_set_adversarial_trapping_events_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/trapping_set/adversarial_trapping_events.csv` |
| `outputs.trapping_set_bfts_wl_hash_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/circuit_level/graph_analysis/trapping_set/bfts_wl_hash_summary.csv` |
| `outputs.trapping_set_bfts_classical_trapping_by_wl_hash_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/circuit_level/graph_analysis/trapping_set/bfts_classical_trapping_by_wl_hash.csv` |
| `outputs.trapping_set_adversarial_trapping_decoder_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/phenomenological/graph_analysis/trapping_set/adversarial_trapping_decoder_summary.csv` |
| `outputs.trapping_set_bfts_symmetric_stabilizer_by_wl_hash_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/circuit_level/graph_analysis/trapping_set/bfts_symmetric_stabilizer_by_wl_hash.csv` |
| `outputs.trapping_set_bfts_trapping_param_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/circuit_level/graph_analysis/trapping_set/bfts_trapping_param_summary.csv` |
| `outputs.trapping_set_bfts_error_node_parity_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/circuit_level/graph_analysis/trapping_set/bfts_error_node_parity_summary.csv` |
| `outputs.trapping_set_bfts_wl_trapping_cross_summary_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats/circuit_level/graph_analysis/trapping_set/bfts_wl_trapping_cross_summary.csv` |
| `templates` | `{}` (no named entries) |

## runtime.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
