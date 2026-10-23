# Configuration for 01 online handoff

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [bp_dynamics_online_osd_handoff_experiments.json](bp_dynamics_online_osd_handoff_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## bp_dynamics_online_osd_handoff_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `schema_version` | `3` | Supported schema/version identifier; preserve the version accepted by the reader. |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `smoke_mode` | `false` | Select the stage's reduced-budget verification path where implemented; not scientific benchmark evidence. |
| `reuse_existing_data` | `true` | Reuse available cached outputs when true; backend/algorithm version checks still apply where implemented. Configuration hashes are not checked. |
| `primary_level` | `2` | Code-family levels/distances evaluated or selected by this stage. |
| `mechanism_levels` | `[2]` | Code-family levels/distances evaluated or selected by this stage. |
| `optional_generalization_levels` | `[3]` | Code-family levels/distances evaluated or selected by this stage. |
| `mechanism_p_values` | `[0.01,0.02,0.03,0.05]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `mechanism_gap_target` | `150` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `mechanism_control_target` | `150` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `mechanism_batch_size` | `1000` | Stage-specific mechanism batch size; used in the notebook's corresponding analysis. |
| `mechanism_max_shots_per_point` | `100000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `minimum_d9_strict_gap_events` | `50` | Stage-specific minimum d9 strict gap events; used in the notebook's corresponding analysis. |
| `final_levels` | `[2]` | Code-family levels/distances evaluated or selected by this stage. |
| `final_p_values` | `[0.002,0.005,0.01,0.02]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `final_min_shots_per_point` | `5000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `final_max_shots_per_point` | `100000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `final_batch_size` | `1000` | Stage-specific final batch size; used in the notebook's corresponding analysis. |
| `final_min_gap_events` | `100` | Stage-specific final min gap events; used in the notebook's corresponding analysis. |
| `bp_max_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `full_osd_method` | `"osd_cs"` | Stage-specific full osd method; used in the notebook's corresponding analysis. |
| `full_osd_order_by_distance` | `{"5":3,"9":1}` | Code-family levels/distances evaluated or selected by this stage. |
| `full_osd_combination_sweep_setting` | `"all_information_set_subsets_weight_leq_order"` | Stage-specific full osd combination sweep setting; used in the notebook's corresponding analysis. |
| `full_osd_candidate_budget_by_distance` | `{"5":1160,"9":33}` | Code-family levels/distances evaluated or selected by this stage. |
| `feature_windows` | `[3,5,8]` | Stage-specific feature windows; used in the notebook's corresponding analysis. |
| `primary_window` | `5` | Stage-specific primary window; used in the notebook's corresponding analysis. |
| `top_k_values` | `[5,10,15]` | Stage-specific top k values; used in the notebook's corresponding analysis. |
| `primary_top_k` | `10` | Stage-specific primary top k; used in the notebook's corresponding analysis. |
| `handoff_thresholds` | `[0.5,0.6,0.7,0.8,0.9,0.95]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `handoff_persistence` | `[1,2,3]` | Stage-specific handoff persistence; used in the notebook's corresponding analysis. |
| `fixed_handoff_iterations` | `[2,4,6,8,10,12,16]` | Iteration budget/window; retain the ordering relationships listed below. |
| `maximum_wait_iterations` | `[8,12,16,60]` | Iteration budget/window; retain the ordering relationships listed below. |
| `prefix_buffer_sizes` | `[3,5,8,16,-1]` | Stage-specific prefix buffer sizes; used in the notebook's corresponding analysis. |
| `adaptive_order_thresholds` | `[0.5,0.6,0.7,0.8,0.9,0.95]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `scaling_values` | `[0.7,0.8,1.0,1.2]` | Stage-specific scaling values; used in the notebook's corresponding analysis. |
| `damping_values` | `[0.0,0.2,0.5]` | Stage-specific damping values; used in the notebook's corresponding analysis. |
| `schedule` | `"parallel"` | Stage-specific schedule; used in the notebook's corresponding analysis. |
| `parameter_scan_case_limit` | `50` | Stage-specific parameter scan case limit; used in the notebook's corresponding analysis. |
| `parameter_scan_unbiased_limit` | `500` | Stage-specific parameter scan unbiased limit; used in the notebook's corresponding analysis. |
| `development_fraction` | `0.5` | Stage-specific development fraction; used in the notebook's corresponding analysis. |
| `cv_folds` | `5` | Cross-validation split setting; preserve grouping and development/test separation. |
| `bootstrap_repetitions` | `1000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `temporal_shuffle_repetitions` | `10` | Stage-specific temporal shuffle repetitions; used in the notebook's corresponding analysis. |
| `random_seed` | `20260910` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/osd/stages/01_online_handoff/plots/principle_fixes_v1` |
| `outputs.config_json` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/config/config.json` |
| `outputs.mechanism_population_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/mechanism_population.csv` |
| `outputs.mechanism_cohort_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/mechanism_cohort.csv` |
| `outputs.final_unbiased_cases_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/final_unbiased_cases.csv` |
| `outputs.qmin_snapshots_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/qmin_snapshots.csv` |
| `outputs.online_snapshot_features_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/online_snapshot_features.csv` |
| `outputs.qmin_temporal_summary_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/qmin_temporal_summary.csv` |
| `outputs.prefix_oracle_ceiling_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/prefix_oracle_ceiling.csv` |
| `outputs.snapshot_feature_ablation_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/snapshot_feature_ablation.csv` |
| `outputs.temporal_shuffle_control_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/temporal_shuffle_control.csv` |
| `outputs.online_handoff_rules_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/online_handoff_rules.csv` |
| `outputs.prefix_buffer_selection_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/prefix_buffer_selection.csv` |
| `outputs.adaptive_osd_order_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/adaptive_osd_order.csv` |
| `outputs.instability_value_ablation_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/instability_value_ablation.csv` |
| `outputs.bp_soft_generator_parameter_scan_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/bp_soft_generator_parameter_scan.csv` |
| `outputs.endpoint_vs_final_osd_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/endpoint_vs_final_osd.csv` |
| `outputs.final_gap_benchmark_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/final_gap_benchmark.csv` |
| `outputs.accuracy_cost_frontier_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/accuracy_cost_frontier.csv` |
| `outputs.prefix_buffer_trials_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/prefix_buffer_trials.csv` |
| `outputs.online_handoff_trials_csv` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace/online_handoff_trials.csv` |
| `templates.name` | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats/{name}` |

Placeholder arguments: `name`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- mechanism_p_values, final_p_values: each item is in (0, 1).
- primary_window: chosen from feature_windows.
- primary_top_k: chosen from top_k_values.
- scaling_values: each item is positive.
- development_fraction: in (0, 1).

Environment variables read literally in these cells: `BP_OSD_HANDOFF_REBUILD`, `BP_OSD_HANDOFF_SMOKE`. Initialization code defines their precedence and parsing.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

`algorithm_versions.osd` is computed as `cpp_cs_v1`. Candidate costs use original channel priors even when an experiment deliberately supplies a different ordering score. For `k = n - rank(H)` and order `q`, the tested candidate budget is `1` for `q=0`, otherwise `1 + k + choose(min(q,k),2)`. Negative orders, non-finite scores and priors outside `(0,1)` are rejected. The postprocessor supports the benchmark `osd_cs` method; unsupported methods are rejected rather than silently mapped to another search.
