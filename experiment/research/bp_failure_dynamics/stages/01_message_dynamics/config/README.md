# Configuration for 01 message dynamics

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [bp_stability_syndrome_logical.json](bp_stability_syndrome_logical.json) | Named directly by the current notebook code. |
| [message.json](message.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## bp_stability_syndrome_logical.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `benchmark_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `diagnostic_iterations` | `50` | Iteration budget/window; retain the ordering relationships listed below. |
| `robustness_iterations` | `100` | Iteration budget/window; retain the ordering relationships listed below. |
| `ms_scaling_factor` | `1.0` | Stage-specific ms scaling factor; used in the notebook's corresponding analysis. |
| `late_fraction` | `0.2` | Stage-specific late fraction; used in the notebook's corresponding analysis. |
| `message_zero_tolerance` | `1e-10` | Stage-specific message zero tolerance; used in the notebook's corresponding analysis. |
| `stability_thresholds` | `[1e-08,1e-06,0.0001,0.001,0.01]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `threshold_quantile_count` | `40` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `validation_cases` | `12` | Stage-specific validation cases; used in the notebook's corresponding analysis. |
| `native_validation_tolerance` | `1e-10` | Stage-specific native validation tolerance; used in the notebook's corresponding analysis. |
| `diagnostic_candidate_limit` | `null` | Stage-specific diagnostic candidate limit; used in the notebook's corresponding analysis. |
| `matched_cases_per_class` | `500` | Stage-specific matched cases per class; used in the notebook's corresponding analysis. |
| `random_seed` | `20260906` | Seed for the configured random support/decoder experiment; retain it with the result. |

## message.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `diagnostic_iterations` | `50` | Iteration budget/window; retain the ordering relationships listed below. |
| `benchmark_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `ms_scaling_factor` | `1.0` | Stage-specific ms scaling factor; used in the notebook's corresponding analysis. |
| `late_fraction` | `0.2` | Stage-specific late fraction; used in the notebook's corresponding analysis. |
| `message_zero_tolerance` | `1e-10` | Stage-specific message zero tolerance; used in the notebook's corresponding analysis. |
| `hidden_belief_threshold` | `0.01` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `hidden_message_threshold` | `0.1` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `unstable_normalized_residual_floor` | `0.001` | Stage-specific unstable normalized residual floor; used in the notebook's corresponding analysis. |
| `top_edge_count` | `20` | Stage-specific top edge count; used in the notebook's corresponding analysis. |
| `unstable_edge_quantile` | `0.9` | Stage-specific unstable edge quantile; used in the notebook's corresponding analysis. |
| `cases_per_class` | `500` | Stage-specific cases per class; used in the notebook's corresponding analysis. |
| `candidate_scan_limit` | `10000` | Stage-specific candidate scan limit; used in the notebook's corresponding analysis. |
| `random_seed` | `20260904` | Seed for the configured random support/decoder experiment; retain it with the result. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `outputs.run` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main` |
| `outputs.stats` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats` |
| `outputs.trace` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace` |
| `outputs.config` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/config` |
| `outputs.plots` | `experiment/research/bp_failure_dynamics/stages/01_message_dynamics/plots/main` |
| `outputs.weight2_data_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/weight2_data_audit.csv` |
| `outputs.representative_case_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/representative_case_summary.csv` |
| `outputs.matched_sample_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/matched_sample_metrics.csv` |
| `outputs.failure_success_group_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/failure_success_group_summary.csv` |
| `outputs.unique_trajectory_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/unique_trajectory_metrics.csv` |
| `outputs.matched_failure_enrichment_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/matched_failure_enrichment.csv` |
| `outputs.representative_hidden_variables_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/representative_hidden_variables.csv` |
| `outputs.hidden_instability_group_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/hidden_instability_group_summary.csv` |
| `outputs.representative_unstable_subgraph_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/representative_unstable_subgraph_summary.csv` |
| `outputs.stability_input_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/stability_input_audit.csv` |
| `outputs.instrumented_native_validation_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/instrumented_native_validation_summary.csv` |
| `outputs.diagnostic_population_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/diagnostic_population_metrics.csv` |
| `outputs.exact_error_match_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/exact_error_match_audit.csv` |
| `outputs.syndrome_consistent_vs_residual_threshold_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/syndrome_consistent_vs_residual_threshold.csv` |
| `outputs.four_outcome_group_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/four_outcome_group_summary.csv` |
| `outputs.syndrome_vs_logical_separation_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/syndrome_vs_logical_separation.csv` |
| `outputs.matched_logical_wrong_enrichment_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/matched_logical_wrong_enrichment.csv` |
| `outputs.stable_valid_logical_wrong_fraction_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/stable_valid_logical_wrong_fraction.csv` |
| `outputs.weight2_failure_trajectory_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/weight2_failure_trajectory_summary.csv` |
| `outputs.weight2_failure_iteration_trajectories_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/weight2_failure_iteration_trajectories.csv` |
| `outputs.residual_prediction_auc_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/residual_prediction_auc.csv` |
| `outputs.exact_error_match_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/exact_error_match_summary.csv` |
| `outputs.success_edge_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/success_edge_metrics.csv` |
| `outputs.failure_message_history_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/failure_message_history.csv` |
| `outputs.failure_edge_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/failure_edge_metrics.csv` |
| `outputs.success_variable_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/success_variable_metrics.csv` |
| `outputs.success_message_history_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/success_message_history.csv` |
| `outputs.failure_variable_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/failure_variable_metrics.csv` |
| `templates.name__message_history_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/{name}_message_history.csv` |
| `templates.name__edge_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/{name}_edge_metrics.csv` |
| `templates.name__variable_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace/{name}_variable_metrics.csv` |

Placeholder arguments: `name`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- diagnostic_iterations: at least 3.
- benchmark_iterations: positive.
- ms_scaling_factor: positive.
- late_fraction: in (0, 1].
- message_zero_tolerance: non-negative.
- top_edge_count: positive.
- hidden_belief_threshold: non-negative.
- hidden_message_threshold: positive.
- unstable_normalized_residual_floor: positive.
- unstable_edge_quantile: in (0, 1).
- Expected relationship: candidate_scan_limit >= 2 * cases_per_class.
- Expected relationship: benchmark_iterations >= diagnostic_iterations >= 3.
- Expected relationship: robustness_iterations >= benchmark_iterations.
- stability_thresholds: each item is positive.
- threshold_quantile_count: at least 2.
- validation_cases: positive.
- native_validation_tolerance: positive.
- diagnostic_candidate_limit: None or positive.
- matched_cases_per_class: positive.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.
