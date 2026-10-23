# Configuration for 02 intervention gating

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [paths.json](paths.json) | Named directly by the current notebook code. |
| [two_stage_bposd_intervention_snapshot_experiments.json](two_stage_bposd_intervention_snapshot_experiments.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/osd/stages/02_intervention_gating/plots/principle_fixes_v1` |
| `outputs.config_json` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/config/config.json` |
| `outputs.unbiased_intervention_cases_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace/unbiased_intervention_cases.csv` |
| `outputs.forced_osd_reconciliation_trials_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace/forced_osd_reconciliation_trials.csv` |
| `outputs.forced_osd_mismatches_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace/forced_osd_mismatches.csv` |
| `outputs.forced_osd_reconciliation_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/forced_osd_reconciliation.csv` |
| `outputs.intervention_label_summary_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/intervention_label_summary.csv` |
| `outputs.stage1_gate_ablation_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/stage1_gate_ablation.csv` |
| `outputs.stage1_calibration_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/stage1_calibration.csv` |
| `outputs.risk_stratified_intervention_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/risk_stratified_intervention.csv` |
| `outputs.intervention_threshold_frontier_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/intervention_threshold_frontier.csv` |
| `outputs.provisional_full_osd_gate_selection_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/provisional_full_osd_gate_selection.csv` |
| `outputs.stage2_snapshot_features_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace/stage2_snapshot_features.csv` |
| `outputs.stage2_snapshot_prediction_diagnostic_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/stage2_snapshot_prediction_diagnostic.csv` |
| `outputs.stage2_selector_ablation_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/stage2_selector_ablation.csv` |
| `outputs.stage2_prefix_oracle_fixed_t_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/stage2_prefix_oracle_fixed_t.csv` |
| `outputs.stage2_temporal_shuffle_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/stage2_temporal_shuffle.csv` |
| `outputs.prefix_buffer_recovery_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/prefix_buffer_recovery.csv` |
| `outputs.cost_sensitive_gate_selection_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/cost_sensitive_gate_selection.csv` |
| `outputs.deployed_two_stage_gate_frontier_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/deployed_two_stage_gate_frontier.csv` |
| `outputs.gated_adaptive_osd_order_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/gated_adaptive_osd_order.csv` |
| `outputs.recovery_policy_comparison_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/recovery_policy_comparison.csv` |
| `outputs.endpoint_vs_final_transitions_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/endpoint_vs_final_transitions.csv` |
| `outputs.bp_parameter_two_stage_scan_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/bp_parameter_two_stage_scan.csv` |
| `outputs.final_two_stage_benchmark_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/final_two_stage_benchmark.csv` |
| `outputs.two_stage_instability_ablation_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/two_stage_instability_ablation.csv` |
| `outputs.accuracy_cost_frontier_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/accuracy_cost_frontier.csv` |
| `outputs.experiment_status_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/experiment_status.csv` |
| `outputs.stage2_selector_trials_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace/stage2_selector_trials.csv` |
| `outputs.two_stage_decoder_trials_csv` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace/two_stage_decoder_trials.csv` |
| `templates.name` | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats/{name}` |

Placeholder arguments: `name`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## two_stage_bposd_intervention_snapshot_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `schema_version` | `8` | Supported schema/version identifier; preserve the version accepted by the reader. |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `smoke_mode` | `false` | Select the stage's reduced-budget verification path where implemented; not scientific benchmark evidence. |
| `reuse_existing_data` | `true` | Reuse available cached outputs when true; backend/algorithm version checks still apply where implemented. Configuration hashes are not checked. |
| `levels` | `[2]` | Code-family levels/distances evaluated or selected by this stage. |
| `optional_d9_generalization` | `false` | Boolean stage option; true enables and false disables the named behavior in the code. |
| `p_values` | `[0.005,0.01,0.02,0.03,0.05]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `minimum_shots_per_point` | `20000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `maximum_shots_per_point` | `100000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `sampling_batch_size` | `2000` | Stage-specific sampling batch size; used in the notebook's corresponding analysis. |
| `minimum_unsafe_events` | `100` | Stage-specific minimum unsafe events; used in the notebook's corresponding analysis. |
| `reconciliation_shots_per_point` | `5000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `reference_audit_shots_per_point` | `20` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `bp_max_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `full_osd_order_by_distance` | `{"5":3,"9":1}` | Code-family levels/distances evaluated or selected by this stage. |
| `full_osd_candidate_budget_by_distance` | `{"5":23,"9":33}` | Code-family levels/distances evaluated or selected by this stage. |
| `full_osd_combination_sweep` | `"cpp_combination_sweep"` | Stage-specific full osd combination sweep; used in the notebook's corresponding analysis. |
| `fixed_snapshot_iterations` | `[1,2,3,4,5,6,8,10]` | Iteration budget/window; retain the ordering relationships listed below. |
| `feature_windows` | `[3,5,8]` | Stage-specific feature windows; used in the notebook's corresponding analysis. |
| `primary_window` | `5` | Stage-specific primary window; used in the notebook's corresponding analysis. |
| `top_k_values` | `[5,10,15]` | Stage-specific top k values; used in the notebook's corresponding analysis. |
| `primary_top_k` | `10` | Stage-specific primary top k; used in the notebook's corresponding analysis. |
| `risk_percentile_edges` | `[0.0,0.5,0.8,0.9,0.95,0.99,1.0]` | Stage-specific risk percentile edges; used in the notebook's corresponding analysis. |
| `gate_thresholds` | `[0.5,0.7,0.8,0.9,0.95,0.98,0.99]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `utility_penalties` | `[1,2,5,10,20]` | Stage-specific utility penalties; used in the notebook's corresponding analysis. |
| `primary_utility_penalty` | `10` | Stage-specific primary utility penalty; used in the notebook's corresponding analysis. |
| `prefix_buffer_sizes` | `[1,3,5,8,-1]` | Stage-specific prefix buffer sizes; used in the notebook's corresponding analysis. |
| `adaptive_order_thresholds` | `[0.5,0.7,0.8,0.9,0.95]` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `development_fit_fraction` | `0.5` | Stage-specific development fit fraction; used in the notebook's corresponding analysis. |
| `development_calibration_fraction` | `0.25` | Stage-specific development calibration fraction; used in the notebook's corresponding analysis. |
| `stage2_selector_fit_fraction` | `0.5` | Stage-specific stage2 selector fit fraction; used in the notebook's corresponding analysis. |
| `calibration_bins` | `10` | Stage-specific calibration bins; used in the notebook's corresponding analysis. |
| `bootstrap_repetitions` | `1000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `temporal_shuffle_repetitions` | `10` | Stage-specific temporal shuffle repetitions; used in the notebook's corresponding analysis. |
| `scaling_values` | `[0.7,0.8,1.0,1.2]` | Stage-specific scaling values; used in the notebook's corresponding analysis. |
| `damping_values` | `[0.0,0.2,0.5]` | Stage-specific damping values; used in the notebook's corresponding analysis. |
| `parameter_scan_case_limit` | `100` | Stage-specific parameter scan case limit; used in the notebook's corresponding analysis. |
| `random_seed` | `20260911` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- p_values: each item is in (0, 1).
- Expected relationship: minimum_shots_per_point <= maximum_shots_per_point.
- primary_window: chosen from feature_windows.
- primary_top_k: chosen from top_k_values.
- primary_utility_penalty: chosen from utility_penalties.
- Expected relationship: development_fit_fraction + development_calibration_fraction < 1.
- stage2_selector_fit_fraction: in (0, 1).

Environment variables read literally in these cells: `TWO_STAGE_BPOSD_REBUILD`, `TWO_STAGE_BPOSD_SMOKE`. Initialization code defines their precedence and parsing.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

`algorithm_versions.osd` is computed as `cpp_cs_v1`. Candidate costs use original channel priors even when an experiment deliberately supplies a different ordering score. For `k = n - rank(H)` and order `q`, the tested candidate budget is `1` for `q=0`, otherwise `1 + k + choose(min(q,k),2)`. Negative orders, non-finite scores and priors outside `(0,1)` are rejected. The postprocessor supports the benchmark `osd_cs` method; unsupported methods are rejected rather than silently mapped to another search.

`full_osd_combination_sweep` is `cpp_combination_sweep`; `full_osd_candidate_budget_by_distance` is `{"5":23,"9":33}` for the current contexts/orders. The context initializer asserts the actual budget. The main reference uses signed BP values; an explicitly requested absolute-order audit remains a separate convention.
