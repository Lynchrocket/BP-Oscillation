# Configuration for 00 gap exploration

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [bposd_tesseract_gap_next_experiments.json](bposd_tesseract_gap_next_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## bposd_tesseract_gap_next_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `smoke_mode` | `false` | Select the stage's reduced-budget verification path where implemented; not scientific benchmark evidence. |
| `reuse_existing_data` | `true` | Reuse available cached outputs when true; backend/algorithm version checks still apply where implemented. Configuration hashes are not checked. |
| `levels` | `[2,3]` | Code-family levels/distances evaluated or selected by this stage. |
| `p_values` | `[0.0003,0.0006,0.001,0.002,0.003,0.005]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `shots_per_point` | `2000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `random_seed` | `20260909` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `bp_max_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `long_bp_iterations` | `100` | Iteration budget/window; retain the ordering relationships listed below. |
| `snapshot_iterations` | `[1,2,4,6,8,10,12,16,20,30,40,60]` | Iteration budget/window; retain the ordering relationships listed below. |
| `osd_orders_by_level` | `{"1":[0,1,9],"2":[0,1,3],"3":[0,1]}` | Code-family levels/distances evaluated or selected by this stage. |
| `snapshot_cases_per_stratum` | `200` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `trajectory_cases_per_stratum` | `500` | Stage-specific trajectory cases per stratum; used in the notebook's corresponding analysis. |
| `patience_values` | `[1,2,4,8,16]` | Stage-specific patience values; used in the notebook's corresponding analysis. |
| `top_k_values` | `[5,10,20]` | Stage-specific top k values; used in the notebook's corresponding analysis. |
| `development_fraction` | `0.5` | Stage-specific development fraction; used in the notebook's corresponding analysis. |
| `cv_folds` | `5` | Cross-validation split setting; preserve grouping and development/test separation. |
| `trigger_quantiles` | `[0.5,0.6,0.7,0.8,0.9,0.95]` | Stage-specific trigger quantiles; used in the notebook's corresponding analysis. |
| `scaling_values` | `[0.7,0.8,1.0,1.2]` | Stage-specific scaling values; used in the notebook's corresponding analysis. |
| `damping_values` | `[0.0,0.2,0.5]` | Stage-specific damping values; used in the notebook's corresponding analysis. |
| `schedule_values` | `["parallel","serial"]` | Stage-specific schedule values; used in the notebook's corresponding analysis. |
| `iteration_values` | `[10,20,40,60]` | Iteration budget/window; retain the ordering relationships listed below. |
| `handoff_grid_cases` | `200` | Stage-specific handoff grid cases; used in the notebook's corresponding analysis. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/gap_mechanisms/stages/00_gap_exploration/plots/principle_fixes_v1` |
| `outputs.config_json` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/config/config.json` |
| `outputs.gap_decomposition_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/gap_decomposition.csv` |
| `outputs.reliability_snapshots_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/reliability_snapshots.csv` |
| `outputs.bp_valid_candidates_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/bp_valid_candidates.csv` |
| `outputs.forced_osd_order_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/forced_osd_order_trials.csv` |
| `outputs.gap_class_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/gap_class_summary.csv` |
| `outputs.gap_closure_fraction_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/gap_closure_fraction.csv` |
| `outputs.gap_decoder_ler_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/gap_decoder_ler.csv` |
| `outputs.forced_osd_order_ablation_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/forced_osd_order_ablation.csv` |
| `outputs.reliability_snapshot_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/reliability_snapshot_summary.csv` |
| `outputs.relative_snapshot_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/relative_snapshot_summary.csv` |
| `outputs.reliability_ranking_quality_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/reliability_ranking_quality.csv` |
| `outputs.bp_trajectory_oracle_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/bp_trajectory_oracle_summary.csv` |
| `outputs.bp_stopping_rules_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/bp_stopping_rules.csv` |
| `outputs.false_stop_models_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/false_stop_models.csv` |
| `outputs.selective_osd_frontier_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/selective_osd_frontier.csv` |
| `outputs.osd0_audit_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/osd0_audit.csv` |
| `outputs.osd0_escalation_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/osd0_escalation.csv` |
| `outputs.bp_soft_information_value_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/bp_soft_information_value.csv` |
| `outputs.experiment_status_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/experiment_status.csv` |
| `outputs.bp_handoff_optimization_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/bp_handoff_optimization.csv` |
| `outputs.instability_aware_handoff_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/instability_aware_handoff.csv` |
| `outputs.final_benchmark_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/final_benchmark.csv` |
| `outputs.selected_selective_policy_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats/selected_selective_policy.csv` |
| `outputs.bp_trajectory_oracle_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/bp_trajectory_oracle.csv` |
| `outputs.final_decoder_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/final_decoder_trials.csv` |
| `outputs.bp_stopping_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/bp_stopping_trials.csv` |
| `outputs.bp_handoff_grid_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/bp_handoff_grid_trials.csv` |
| `outputs.false_stop_scores_csv` | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace/false_stop_scores.csv` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- p_values: each item is in (0, 1).
- shots_per_point: positive.
- Expected relationship: bp_max_iterations >= maximum of snapshot_iterations.
- Expected relationship: long_bp_iterations >= bp_max_iterations.
- development_fraction: in (0, 1).
- cv_folds: at least 2.
- patience_values, top_k_values, iteration_values: each item is positive.
- scaling_values: each item is positive.
- snapshot_cases_per_stratum: positive.
- trajectory_cases_per_stratum: positive.
- handoff_grid_cases: positive.
- tesseract_beam: positive.
- tesseract_pq_limit: positive.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

`algorithm_versions.osd` is computed as `cpp_cs_v1`. Candidate costs use original channel priors even when an experiment deliberately supplies a different ordering score. For `k = n - rank(H)` and order `q`, the tested candidate budget is `1` for `q=0`, otherwise `1 + k + choose(min(q,k),2)`. Negative orders, non-finite scores and priors outside `(0,1)` are rejected. The postprocessor supports the benchmark `osd_cs` method; unsupported methods are rejected rather than silently mapped to another search.
