# Configuration for 02 failure mechanisms

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [bp_failure_mechanisms.json](bp_failure_mechanisms.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## bp_failure_mechanisms.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `benchmark_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `diagnostic_iterations` | `50` | Iteration budget/window; retain the ordering relationships listed below. |
| `long_iterations` | `100` | Iteration budget/window; retain the ordering relationships listed below. |
| `ms_scaling_factor` | `1.0` | Stage-specific ms scaling factor; used in the notebook's corresponding analysis. |
| `late_fraction` | `0.2` | Stage-specific late fraction; used in the notebook's corresponding analysis. |
| `message_zero_tolerance` | `1e-10` | Stage-specific message zero tolerance; used in the notebook's corresponding analysis. |
| `stability_thresholds` | `[1e-08,1e-06,0.0001,0.001,0.01]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `primary_stability_threshold` | `1e-08` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `fixed_patience_values` | `[1,2,4,8,16]` | Stage-specific fixed patience values; used in the notebook's corresponding analysis. |
| `logical_stability_values` | `[2,3,4,5]` | Stage-specific logical stability values; used in the notebook's corresponding analysis. |
| `residual_thresholds` | `[1e-08,1e-06,0.0001,0.001,0.01]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `adaptive_k` | `8` | Stage-specific adaptive k; used in the notebook's corresponding analysis. |
| `adaptive_max_views` | `4` | Stage-specific adaptive max views; used in the notebook's corresponding analysis. |
| `d9_weight4_bp_failure_samples_per_p` | `1000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `d9_weight4_reuse_cached_bp_failures` | `true` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `equal_budget_runtime_repeats` | `1` | Stage-specific equal budget runtime repeats; used in the notebook's corresponding analysis. |
| `random_view_seed` | `20260907` | Randomness/repetition setting; use the initializer's stated seed semantics. |
| `cross_distance_levels` | `[1,2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_error_weights` | `[2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `cross_distance_weight3_sample_limit` | `2000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `cross_distance_reuse_distance5` | `true` | Code-family levels/distances evaluated or selected by this stage. |
| `validation_cases` | `12` | Stage-specific validation cases; used in the notebook's corresponding analysis. |
| `native_validation_tolerance` | `1e-10` | Stage-specific native validation tolerance; used in the notebook's corresponding analysis. |
| `random_seed` | `20260906` | Seed for the configured random support/decoder experiment; retain it with the result. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `outputs.run` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/bp_failure_dynamics/stages/02_failure_mechanisms/plots/principle_fixes_v1` |
| `outputs.next_stage_input_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/next_stage_input_audit.csv` |
| `outputs.next_stage_native_validation_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/next_stage_native_validation_summary.csv` |
| `outputs.failure_type_population_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/failure_type_population.csv` |
| `outputs.failure_input_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/failure_input_summary.csv` |
| `outputs.logical_sector_trajectories_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/logical_sector_trajectories.csv` |
| `outputs.residual_operator_metrics_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/residual_operator_metrics.csv` |
| `outputs.failure_type_population_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/failure_type_population_summary.csv` |
| `outputs.logical_delay_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/logical_delay_summary.csv` |
| `outputs.logical_sector_transition_matrix_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/logical_sector_transition_matrix.csv` |
| `outputs.stopping_rule_benchmark_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/stopping_rule_benchmark.csv` |
| `outputs.stopping_rule_repair_by_failure_type_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/stopping_rule_repair_by_failure_type.csv` |
| `outputs.stable_fixed_point_comparison_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/stable_fixed_point_comparison.csv` |
| `outputs.stable_structure_matched_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/stable_structure_matched_cases.csv` |
| `outputs.stable_wrong_structure_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/stable_wrong_structure_summary.csv` |
| `outputs.adaptive_repair_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/adaptive_repair_cases.csv` |
| `outputs.adaptive_repair_by_failure_type_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/adaptive_repair_by_failure_type.csv` |
| `outputs.d9_weight4_bp_failure_sample_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/d9_weight4_bp_failure_sample.csv` |
| `outputs.d9_weight4_bp_failure_sampling_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/d9_weight4_bp_failure_sampling_audit.csv` |
| `outputs.d9_weight4_bp_failure_type_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/d9_weight4_bp_failure_type_summary.csv` |
| `outputs.d9_weight4_bp_adaptive_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/d9_weight4_bp_adaptive_cases.csv` |
| `outputs.d9_weight4_bp_adaptive_repair_by_failure_type_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/d9_weight4_bp_adaptive_repair_by_failure_type.csv` |
| `outputs.adaptive_equal_budget_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/adaptive_equal_budget_cases.csv` |
| `outputs.adaptive_equal_budget_comparison_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/adaptive_equal_budget_comparison.csv` |
| `outputs.adaptive_equal_budget_by_failure_type_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/adaptive_equal_budget_by_failure_type.csv` |
| `outputs.cross_distance_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace/cross_distance_cases.csv` |
| `outputs.cross_distance_input_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/cross_distance_input_audit.csv` |
| `outputs.failure_type_vs_distance_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/failure_type_vs_distance.csv` |
| `outputs.cross_distance_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/cross_distance_summary.csv` |
| `outputs.cross_distance_adaptive_repair_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats/cross_distance_adaptive_repair.csv` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- Expected relationship: long_iterations >= benchmark_iterations >= diagnostic_iterations >= 3.
- ms_scaling_factor: positive.
- late_fraction: in (0, 1].
- primary_stability_threshold: chosen from stability_thresholds.
- stability_thresholds, residual_thresholds: each item is positive.
- fixed_patience_values, logical_stability_values: each item is positive.
- adaptive_k: non-negative.
- adaptive_max_views: non-negative.
- d9_weight4_bp_failure_samples_per_p: positive.
- equal_budget_runtime_repeats: at least 1.
- cross_distance_error_weights: values exactly {2, 3}.
- cross_distance_weight3_sample_limit: positive.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

| Key | Allowed values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `restart_message_sigma` | Finite non-negative number | `0.05` | Standard deviation of additive Gaussian noise on initial V2C messages. Zero recovers the old deterministic initialization; priors and sweep budgets do not change. |
| `algorithm_versions.restart` | `gaussian_messages_v1` | Computed | Identifies randomized restart semantics; incompatible with unmarked historical results. |

Seeds are derived from `random_seed`, sample ID and restart ID through NumPy `SeedSequence`. Per-case output adds `restart_message_sigma` and `restart_seeds_json` for the randomized control; other variants store no restart seeds. Existing variant names and budget fields remain available.
