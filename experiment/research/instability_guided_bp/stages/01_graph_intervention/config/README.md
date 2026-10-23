# Configuration for 01 graph intervention

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [instability_guided_graph_intervention_experiments.json](instability_guided_graph_intervention_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## instability_guided_graph_intervention_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `benchmark_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `diagnostic_iterations` | `50` | Iteration budget/window; retain the ordering relationships listed below. |
| `long_iterations` | `100` | Iteration budget/window; retain the ordering relationships listed below. |
| `probe_iterations` | `8` | Iteration budget/window; retain the ordering relationships listed below. |
| `candidate_test_iterations` | `[3,5]` | Iteration budget/window; retain the ordering relationships listed below. |
| `primary_candidate_test_iterations` | `5` | Iteration budget/window; retain the ordering relationships listed below. |
| `late_fraction` | `0.2` | Stage-specific late fraction; used in the notebook's corresponding analysis. |
| `message_zero_tolerance` | `1e-10` | Stage-specific message zero tolerance; used in the notebook's corresponding analysis. |
| `native_validation_tolerance` | `1e-10` | Stage-specific native validation tolerance; used in the notebook's corresponding analysis. |
| `ms_scaling_factor` | `1.0` | Stage-specific ms scaling factor; used in the notebook's corresponding analysis. |
| `stability_threshold` | `1e-08` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `top_k_values` | `[4,8,12,16]` | Stage-specific top k values; used in the notebook's corresponding analysis. |
| `primary_top_k` | `8` | Stage-specific primary top k; used in the notebook's corresponding analysis. |
| `composite_weights` | `[0.5,0.25,0.25]` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `candidate_max_row_weights` | `[4,5,6,8]` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `primary_max_row_weight` | `8` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `candidate_screen_limit` | `30` | Stage-specific candidate screen limit; used in the notebook's corresponding analysis. |
| `samples_per_outcome` | `100` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `random_support_repetitions` | `100` | Stage-specific random support repetitions; used in the notebook's corresponding analysis. |
| `targeting_random_seeds` | `5` | Randomness/repetition setting; use the initializer's stated seed semantics. |
| `gamma_min_values` | `[0.0,0.1]` | Stage-specific gamma min values; used in the notebook's corresponding analysis. |
| `gamma_max_values` | `[0.3,0.5,0.7,0.9]` | Stage-specific gamma max values; used in the notebook's corresponding analysis. |
| `primary_gamma_min` | `0.1` | Stage-specific primary gamma min; used in the notebook's corresponding analysis. |
| `primary_gamma_max` | `0.7` | Stage-specific primary gamma max; used in the notebook's corresponding analysis. |
| `alpha_min_values` | `[0.5,0.7,0.8]` | Stage-specific alpha min values; used in the notebook's corresponding analysis. |
| `alpha_max` | `1.0` | Stage-specific alpha max; used in the notebook's corresponding analysis. |
| `primary_alpha_min` | `0.7` | Stage-specific primary alpha min; used in the notebook's corresponding analysis. |
| `adaptation_rounds` | `[2,3,4]` | Stage-specific adaptation rounds; used in the notebook's corresponding analysis. |
| `cross_distance_levels` | `[1,2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_error_weights` | `[2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_cases_per_stratum` | `40` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_p_limit` | `6` | Code-family levels/distances evaluated or selected by this stage. |
| `bootstrap_repetitions` | `2000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `cv_folds` | `5` | Cross-validation split setting; preserve grouping and development/test separation. |
| `random_seed` | `20260908` | Seed for the configured random support/decoder experiment; retain it with the result. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `outputs.run` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main` |
| `outputs.stats` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats` |
| `outputs.trace` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace` |
| `outputs.config` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/config` |
| `outputs.plots` | `experiment/research/instability_guided_bp/stages/01_graph_intervention/plots/main` |
| `outputs.transient_instability_summary_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/transient_instability_summary.csv` |
| `outputs.graph_intervention_selection_audit_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/graph_intervention_selection_audit.csv` |
| `outputs.instability_residual_overlap_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/instability_residual_overlap.csv` |
| `outputs.candidate_pool_summary_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/candidate_pool_summary.csv` |
| `outputs.candidate_row_gain_summary_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/candidate_row_gain_summary.csv` |
| `outputs.repair_vs_instability_gain_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/repair_vs_instability_gain.csv` |
| `outputs.targeting_ablation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/targeting_ablation.csv` |
| `outputs.damping_ablation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/damping_ablation.csv` |
| `outputs.scaling_ablation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/scaling_ablation.csv` |
| `outputs.unified_intervention_comparison_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/unified_intervention_comparison.csv` |
| `outputs.iterative_adaptation_summary_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/iterative_adaptation_summary.csv` |
| `outputs.repair_by_failure_type_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/repair_by_failure_type.csv` |
| `outputs.cross_distance_input_audit_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/cross_distance_input_audit.csv` |
| `outputs.cross_distance_instability_targeting_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/cross_distance_instability_targeting.csv` |
| `outputs.graph_intervention_native_validation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/graph_intervention_native_validation.csv` |
| `outputs.targeting_paired_comparison_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/targeting_paired_comparison.csv` |
| `outputs.repair_gain_model_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/repair_gain_model.csv` |
| `outputs.instability_residual_overlap_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/instability_residual_overlap_cases.csv` |
| `outputs.transient_variable_instability_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/transient_variable_instability.csv` |
| `outputs.unified_intervention_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/unified_intervention_cases.csv` |
| `outputs.transient_instability_iterations_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/transient_instability_iterations.csv` |
| `outputs.cross_distance_instability_targeting_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/cross_distance_instability_targeting_cases.csv` |
| `outputs.graph_intervention_population_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/graph_intervention_population_cases.csv` |
| `outputs.candidate_row_trials_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/candidate_row_trials.csv` |
| `outputs.damping_intervention_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/damping_intervention_cases.csv` |
| `outputs.iterative_graph_trajectories_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/iterative_graph_trajectories.csv` |
| `outputs.targeting_ablation_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/targeting_ablation_cases.csv` |
| `outputs.scaling_intervention_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/scaling_intervention_cases.csv` |
| `outputs.iterative_adaptation_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace/iterative_adaptation_cases.csv` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- Expected relationship: long_iterations >= benchmark_iterations >= diagnostic_iterations >=
- probe_iterations >= 2.
- primary_candidate_test_iterations: chosen from candidate_test_iterations.
- candidate_test_iterations: each item is within the bounds 2 <= value <= probe_iterations.
- late_fraction: in (0, 1].
- message_zero_tolerance: non-negative.
- native_validation_tolerance: positive.
- ms_scaling_factor: finite and positive.
- stability_threshold: positive.
- primary_top_k: chosen from top_k_values.
- top_k_values: each item is positive.
- Expected relationship: number of entries in composite_weights == 3.
- sum of composite_weights: approximately 1.
- composite_weights: each item is non-negative.
- primary_max_row_weight: chosen from candidate_max_row_weights.
- candidate_screen_limit: positive.
- samples_per_outcome: positive.
- random_support_repetitions: positive.
- targeting_random_seeds: positive.
- cross_distance_levels: values exactly {1, 2, 3}.
- cross_distance_error_weights: values exactly {2, 3}.
- bootstrap_repetitions: positive.
- cv_folds: at least 2.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.
