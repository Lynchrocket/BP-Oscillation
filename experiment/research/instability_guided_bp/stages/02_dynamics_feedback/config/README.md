# Configuration for 02 dynamics feedback

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [next_config.json](next_config.json) | Named directly by the current notebook code. |
| [next_instability_guided_bp.json](next_instability_guided_bp.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## next_config.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `probe_iterations` | `8` | Iteration budget/window; retain the ordering relationships listed below. |
| `candidate_test_iterations` | `[3,5]` | Iteration budget/window; retain the ordering relationships listed below. |
| `primary_candidate_test_iterations` | `5` | Iteration budget/window; retain the ordering relationships listed below. |
| `primary_top_k` | `8` | Stage-specific primary top k; used in the notebook's corresponding analysis. |
| `candidate_max_row_weight` | `8` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `candidate_screen_limit` | `30` | Stage-specific candidate screen limit; used in the notebook's corresponding analysis. |
| `unbiased_candidate_limit` | `50` | Stage-specific unbiased candidate limit; used in the notebook's corresponding analysis. |
| `samples_per_outcome` | `100` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `random_control_seeds` | `5` | Randomness/repetition setting; use the initializer's stated seed semantics. |
| `hybrid_alphas` | `[0.0,0.25,0.5,0.75,1.0]` | Stage-specific hybrid alphas; used in the notebook's corresponding analysis. |
| `development_fraction` | `0.5` | Stage-specific development fraction; used in the notebook's corresponding analysis. |
| `early_probe_iterations` | `[2,4,6,8,12,16]` | Iteration budget/window; retain the ordering relationships listed below. |
| `iterative_max_rounds` | `4` | Stage-specific iterative max rounds; used in the notebook's corresponding analysis. |
| `plateau_delta_values` | `[0.0,0.0001,0.001,0.01]` | Stage-specific plateau delta values; used in the notebook's corresponding analysis. |
| `plateau_patience` | `2` | Stage-specific plateau patience; used in the notebook's corresponding analysis. |
| `gating_quantiles` | `[0.5,0.75,0.9,0.95]` | Stage-specific gating quantiles; used in the notebook's corresponding analysis. |
| `candidate_budget_values` | `[3,5,10,15,20,30]` | Stage-specific candidate budget values; used in the notebook's corresponding analysis. |
| `progressive_gain_quantiles` | `[0.5,0.75,0.9]` | Stage-specific progressive gain quantiles; used in the notebook's corresponding analysis. |
| `progressive_patience_values` | `[3,5]` | Stage-specific progressive patience values; used in the notebook's corresponding analysis. |
| `damping_gamma_min` | `0.1` | Stage-specific damping gamma min; used in the notebook's corresponding analysis. |
| `damping_gamma_max` | `0.7` | Stage-specific damping gamma max; used in the notebook's corresponding analysis. |
| `cross_distance_levels` | `[1,2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_error_weights` | `[2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_p_limit` | `6` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_cases_per_stratum` | `40` | Code-family levels/distances evaluated or selected by this stage. |
| `stochastic_shots_per_point` | `2000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `stochastic_batch_size` | `200` | Stage-specific stochastic batch size; used in the notebook's corresponding analysis. |
| `longer_bp_iterations` | `240` | Iteration budget/window; retain the ordering relationships listed below. |
| `osd_order` | `10` | OSD search order; compare with the configured OSD method and decoder budget. |
| `adaptive_k` | `8` | Stage-specific adaptive k; used in the notebook's corresponding analysis. |
| `adaptive_max_views` | `4` | Stage-specific adaptive max views; used in the notebook's corresponding analysis. |
| `native_validation_cases` | `12` | Stage-specific native validation cases; used in the notebook's corresponding analysis. |
| `native_validation_tolerance` | `1e-10` | Stage-specific native validation tolerance; used in the notebook's corresponding analysis. |
| `next_random_seed` | `20260909` | Randomness/repetition setting; use the initializer's stated seed semantics. |

## next_instability_guided_bp.json

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
| `inputs.instability_graph` | `experiment/data/input_snapshots/instability_guided_bp/02_dynamics_feedback/instability_graph` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.trace` | `experiment/data/input_snapshots/instability_guided_bp/02_dynamics_feedback/instability_graph/trace` |
| `inputs.graph_intervention_population_cases_csv` | `experiment/data/input_snapshots/instability_guided_bp/02_dynamics_feedback/instability_graph/trace/graph_intervention_population_cases.csv` |
| `outputs.run` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main` |
| `outputs.stats` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats` |
| `outputs.trace` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace` |
| `outputs.config` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/config` |
| `outputs.plots` | `experiment/research/instability_guided_bp/stages/02_dynamics_feedback/plots/main` |
| `outputs.next_instability_input_audit_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/next_instability_input_audit.csv` |
| `outputs.raw_candidate_trials_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/raw_candidate_trials.csv` |
| `outputs.raw_gain_vs_repair_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/raw_gain_vs_repair.csv` |
| `outputs.instability_metric_ablation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/instability_metric_ablation.csv` |
| `outputs.overlap_gain_joint_analysis_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/overlap_gain_joint_analysis.csv` |
| `outputs.hybrid_targeting_ablation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/hybrid_targeting_ablation.csv` |
| `outputs.unbiased_candidate_validation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/unbiased_candidate_validation.csv` |
| `outputs.trajectory_redirection_analysis_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/trajectory_redirection_analysis.csv` |
| `outputs.early_intervention_trials_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/early_intervention_trials.csv` |
| `outputs.early_intervention_stable_wrong_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/early_intervention_stable_wrong.csv` |
| `outputs.iterative_feedback_trajectories_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/iterative_feedback_trajectories.csv` |
| `outputs.fixed_region_feedback_summary_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/fixed_region_feedback_summary.csv` |
| `outputs.feedback_stopping_cases_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/feedback_stopping_cases.csv` |
| `outputs.feedback_stopping_rules_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/feedback_stopping_rules.csv` |
| `outputs.intervention_gating_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/intervention_gating.csv` |
| `outputs.candidate_budget_trials_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/candidate_budget_trials.csv` |
| `outputs.candidate_budget_ablation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/candidate_budget_ablation.csv` |
| `outputs.final_decoder_development_tuning_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/final_decoder_development_tuning.csv` |
| `outputs.candidate_structure_analysis_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/candidate_structure_analysis.csv` |
| `outputs.failure_type_mechanism_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/failure_type_mechanism.csv` |
| `outputs.cross_distance_trials_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/cross_distance_trials.csv` |
| `outputs.raw_gain_cross_distance_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/raw_gain_cross_distance.csv` |
| `outputs.final_stochastic_trials_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace/final_stochastic_trials.csv` |
| `outputs.final_stochastic_benchmark_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/final_stochastic_benchmark.csv` |
| `outputs.next_native_validation_csv` | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/next_native_validation.csv` |
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
- Expected relationship: probe_iterations >= maximum of candidate_test_iterations >= 3.
- primary_top_k: positive.
- candidate_max_row_weight: at least 2.
- Expected relationship: candidate_screen_limit >= maximum of candidate_budget_values > 0.
- unbiased_candidate_limit: positive.
- random_control_seeds: positive.
- development_fraction: in (0, 1).
- early_probe_iterations: each item is at least 2.
- Expected relationship: iterative_max_rounds >= plateau_patience > 0.
- plateau_delta_values: each item is non-negative.
- Expected relationship: 0 <= damping_gamma_min <= damping_gamma_max < 1.
- stochastic_shots_per_point: positive.
- stochastic_batch_size: positive.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.
