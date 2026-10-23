# Configuration for 01 reliability dynamics

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [bp_message_dynamics_for_osd_reliability_experiments.json](bp_message_dynamics_for_osd_reliability_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## bp_message_dynamics_for_osd_reliability_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `schema_version` | `3` | Supported schema/version identifier; preserve the version accepted by the reader. |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `smoke_mode` | `false` | Select the stage's reduced-budget verification path where implemented; not scientific benchmark evidence. |
| `reuse_existing_data` | `true` | Reuse available cached outputs when true; backend/algorithm version checks still apply where implemented. Configuration hashes are not checked. |
| `levels` | `[2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `primary_level` | `2` | Code-family levels/distances evaluated or selected by this stage. |
| `instance_weights_by_level` | `{"2":[2,3],"3":[4,5]}` | Code-family levels/distances evaluated or selected by this stage. |
| `event_p_values` | `[0.01,0.02,0.03,0.05]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `event_gap_targets_by_level` | `{"2":150,"3":75}` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `event_batch_size` | `2000` | Stage-specific event batch size; used in the notebook's corresponding analysis. |
| `event_max_shots_per_point` | `100000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `unbiased_p_values` | `[0.002,0.005,0.01,0.02]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `unbiased_shots_per_point` | `5000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `bp_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `snapshot_iterations` | `[1,2,3,4,5,6,8,10,12,16,20,30,40,60]` | Iteration budget/window; retain the ordering relationships listed below. |
| `osd_orders_by_level` | `{"2":[0,1,3],"3":[0,1]}` | Code-family levels/distances evaluated or selected by this stage. |
| `top_k_values` | `[5,10,20]` | Stage-specific top k values; used in the notebook's corresponding analysis. |
| `temporal_windows` | `[3,5,8,12,20]` | Stage-specific temporal windows; used in the notebook's corresponding analysis. |
| `instability_beta_values` | `[0.1,0.25,0.5,1.0,2.0,4.0]` | Stage-specific instability beta values; used in the notebook's corresponding analysis. |
| `sign_gamma_values` | `[0.5,1.0,2.0]` | Stage-specific sign gamma values; used in the notebook's corresponding analysis. |
| `conflict_beta_values` | `[0.25,0.5,1.0,2.0]` | Stage-specific conflict beta values; used in the notebook's corresponding analysis. |
| `rank_threshold_values` | `[0.9,0.95,0.98]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `rank_persistence_values` | `[2,3,4]` | Stage-specific rank persistence values; used in the notebook's corresponding analysis. |
| `snapshot_stability_window` | `3` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `message_residual_quantiles` | `[0.1,0.25,0.5]` | Stage-specific message residual quantiles; used in the notebook's corresponding analysis. |
| `combined_weight_values` | `[0.0,0.5,1.0]` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `development_fraction` | `0.5` | Stage-specific development fraction; used in the notebook's corresponding analysis. |
| `bootstrap_repetitions` | `1000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `random_seed` | `20260909` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `scaling_values` | `[0.7,0.8,1.0,1.2]` | Stage-specific scaling values; used in the notebook's corresponding analysis. |
| `damping_values` | `[0.0,0.2,0.5]` | Stage-specific damping values; used in the notebook's corresponding analysis. |
| `schedule_values` | `["parallel","serial"]` | Stage-specific schedule values; used in the notebook's corresponding analysis. |
| `iteration_values` | `[10,20,40,60]` | Iteration budget/window; retain the ordering relationships listed below. |
| `hyperparameter_case_limit` | `200` | Stage-specific hyperparameter case limit; used in the notebook's corresponding analysis. |
| `message_zero_tolerance` | `1e-10` | Stage-specific message zero tolerance; used in the notebook's corresponding analysis. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/osd/stages/01_reliability_dynamics/plots/principle_fixes_v1` |
| `outputs.config_json` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/config/config.json` |
| `outputs.event_screen_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/event_screen.csv` |
| `outputs.event_cohort_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/event_cohort.csv` |
| `outputs.experiment_status_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/experiment_status.csv` |
| `outputs.event_case_metrics_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/event_case_metrics.csv` |
| `outputs.bp_llr_trajectories_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/bp_llr_trajectories.csv` |
| `outputs.bp_rank_trajectories_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/bp_rank_trajectories.csv` |
| `outputs.bp_message_trajectories_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/bp_message_trajectories.csv` |
| `outputs.variable_dynamics_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/variable_dynamics.csv` |
| `outputs.osd_snapshot_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/osd_snapshot_trials.csv` |
| `outputs.osd_success_vs_iteration_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/osd_success_vs_iteration.csv` |
| `outputs.reliability_rank_dynamics_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/reliability_rank_dynamics.csv` |
| `outputs.instability_rank_correlation_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/instability_rank_correlation.csv` |
| `outputs.relevant_variable_ranks_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/relevant_variable_ranks.csv` |
| `outputs.reliability_ablation_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/reliability_ablation_trials.csv` |
| `outputs.combined_reliability_development_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/combined_reliability_development_trials.csv` |
| `outputs.selected_reliability_parameters_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/selected_reliability_parameters.csv` |
| `outputs.event_reliability_method_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/event_reliability_method_trials.csv` |
| `outputs.handoff_rule_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/handoff_rule_trials.csv` |
| `outputs.handoff_rule_benchmark_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/handoff_rule_benchmark.csv` |
| `outputs.snapshot_selection_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/snapshot_selection_trials.csv` |
| `outputs.snapshot_selection_benchmark_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/snapshot_selection_benchmark.csv` |
| `outputs.bp_parameter_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/bp_parameter_trials.csv` |
| `outputs.bp_parameter_osd_objective_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/bp_parameter_osd_objective.csv` |
| `outputs.unbiased_cases_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/unbiased_cases.csv` |
| `outputs.unbiased_reliability_trials_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace/unbiased_reliability_trials.csv` |
| `outputs.reliability_method_osd_benchmark_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/reliability_method_osd_benchmark.csv` |
| `outputs.bposd_tesseract_gap_closure_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/bposd_tesseract_gap_closure.csv` |
| `outputs.temporal_mean_reliability_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/temporal_mean_reliability.csv` |
| `outputs.instability_penalty_ablation_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/instability_penalty_ablation.csv` |
| `outputs.snr_reliability_ablation_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/snr_reliability_ablation.csv` |
| `outputs.sign_consistency_ablation_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/sign_consistency_ablation.csv` |
| `outputs.temporal_median_reliability_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/temporal_median_reliability.csv` |
| `outputs.conflict_reliability_ablation_csv` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/conflict_reliability_ablation.csv` |
| `templates.output_name` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/{output_name}` |
| `templates.name` | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats/{name}` |

Placeholder arguments: `name`, `output_name`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- primary_level: chosen from levels.
- event_p_values, unbiased_p_values: each item is in (0, 1).
- event_batch_size: positive.
- event_max_shots_per_point: positive.
- unbiased_shots_per_point: positive.
- Expected relationship: bp_iterations >= maximum of snapshot_iterations.
- temporal_windows, top_k_values: each item is positive.
- snapshot_stability_window: positive.
- scaling_values: each item is positive.
- development_fraction: in (0, 1).
- bootstrap_repetitions: positive.
- message_zero_tolerance: non-negative.

Environment variables read literally in these cells: `BP_OSD_RELIABILITY_REBUILD`, `BP_OSD_RELIABILITY_SMOKE`. Initialization code defines their precedence and parsing.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

`algorithm_versions.osd` is computed as `cpp_cs_v1`. Candidate costs use original channel priors even when an experiment deliberately supplies a different ordering score. For `k = n - rank(H)` and order `q`, the tested candidate budget is `1` for `q=0`, otherwise `1 + k + choose(min(q,k),2)`. Negative orders, non-finite scores and priors outside `(0,1)` are rejected. The postprocessor supports the benchmark `osd_cs` method; unsupported methods are rejected rather than silently mapped to another search.
