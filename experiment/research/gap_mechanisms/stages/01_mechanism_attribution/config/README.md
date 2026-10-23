# Configuration for 01 mechanism attribution

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [bposd_tesseract_gap_mechanism_attribution_experiments.json](bposd_tesseract_gap_mechanism_attribution_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## bposd_tesseract_gap_mechanism_attribution_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `schema_version` | `9` | Supported schema/version identifier; preserve the version accepted by the reader. |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `smoke_mode` | `false` | Select the stage's reduced-budget verification path where implemented; not scientific benchmark evidence. |
| `reuse_existing_data` | `true` | Reuse available cached outputs when true; backend/algorithm version checks still apply where implemented. Configuration hashes are not checked. |
| `primary_levels` | `[2]` | Code-family levels/distances evaluated or selected by this stage. |
| `validation_levels` | `[3]` | Code-family levels/distances evaluated or selected by this stage. |
| `p_values` | `[0.005,0.01,0.02,0.03,0.05]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `validation_p_values` | `[0.02,0.03,0.05]` | Physical-noise probabilities evaluated by this stage; comparisons are conditioned on p. |
| `minimum_shots_per_point` | `20000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `maximum_shots_per_point` | `100000` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `sampling_batch_size` | `2000` | Stage-specific sampling batch size; used in the notebook's corresponding analysis. |
| `target_gap_cases_per_point` | `100` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `analysis_gap_limit_per_point` | `100` | Stage-specific analysis gap limit per point; used in the notebook's corresponding analysis. |
| `bp_standard_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `continued_iterations` | `100` | Iteration budget/window; retain the ordering relationships listed below. |
| `configured_full_osd_order_by_distance` | `{"5":3,"9":3}` | Code-family levels/distances evaluated or selected by this stage. |
| `osd_orders` | `[0,1,2]` | Stage-specific osd orders; used in the notebook's corresponding analysis. |
| `pcm_views` | `20` | Stage-specific pcm views; used in the notebook's corresponding analysis. |
| `selected_pcm_views` | `50` | Stage-specific selected pcm views; used in the notebook's corresponding analysis. |
| `selected_pcm_case_limit` | `20` | Stage-specific selected pcm case limit; used in the notebook's corresponding analysis. |
| `pcm_rows_per_view` | `[1,2,3]` | Stage-specific pcm rows per view; used in the notebook's corresponding analysis. |
| `candidate_max_row_weight` | `10` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `restart_count` | `50` | Stage-specific restart count; used in the notebook's corresponding analysis. |
| `restart_message_sigma` | `0.05` | Stage-specific restart message sigma; used in the notebook's corresponding analysis. |
| `restart_prior_sigma` | `0.02` | Stage-specific restart prior sigma; used in the notebook's corresponding analysis. |
| `nms_alpha_values` | `[0.6,0.7,0.8,0.9,1.0]` | Stage-specific nms alpha values; used in the notebook's corresponding analysis. |
| `oms_beta_values` | `[0.1,0.25,0.5,1.0]` | Stage-specific oms beta values; used in the notebook's corresponding analysis. |
| `scaling_values` | `[0.6,0.7,0.8,0.9,1.0]` | Stage-specific scaling values; used in the notebook's corresponding analysis. |
| `damping_values` | `[0.0,0.2,0.5,0.7]` | Stage-specific damping values; used in the notebook's corresponding analysis. |
| `prior_decode_multipliers` | `[0.75,0.9,1.0,1.1,1.25]` | Stage-specific prior decode multipliers; used in the notebook's corresponding analysis. |
| `prior_variable_sigma` | `0.03` | Stage-specific prior variable sigma; used in the notebook's corresponding analysis. |
| `prior_perturbations` | `20` | Stage-specific prior perturbations; used in the notebook's corresponding analysis. |
| `ambiguity_top_k` | `8` | Stage-specific ambiguity top k; used in the notebook's corresponding analysis. |
| `pcm_ensemble_sizes` | `[2,4,8,16]` | Stage-specific pcm ensemble sizes; used in the notebook's corresponding analysis. |
| `adaptive_max_views` | `4` | Stage-specific adaptive max views; used in the notebook's corresponding analysis. |
| `logical_margin_threshold` | `0.25` | Decision/diagnostic cutoff; scope and inequalities are specified by the initializer and analyses. |
| `milp_time_limit_seconds` | `60.0` | Stage-specific milp time limit seconds; used in the notebook's corresponding analysis. |
| `logical_cost_case_limit_per_group` | `100` | Stage-specific logical cost case limit per group; used in the notebook's corresponding analysis. |
| `bootstrap_repetitions` | `1000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `development_fraction` | `0.5` | Stage-specific development fraction; used in the notebook's corresponding analysis. |
| `final_case_limit_per_point` | `-1` | Stage-specific final case limit per point; used in the notebook's corresponding analysis. |
| `random_seed` | `20260912` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |
| `tesseract_timing_shots_per_batch` | `100` | Population/selection budget; a target is not proof that enough eligible events were obtained. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/gap_mechanisms/stages/01_mechanism_attribution/plots/principle_fixes_v1` |
| `outputs.config_json` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/config/config.json` |
| `outputs.paired_population_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/paired_population.csv` |
| `outputs.gap_bp_trajectories_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/gap_bp_trajectories.csv` |
| `outputs.population_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/population_summary.csv` |
| `outputs.experiment_1_early_stop_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_1_early_stop.csv` |
| `outputs.experiment_2_osd_orders_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_2_osd_orders.csv` |
| `outputs.osd_order_decomposition_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/osd_order_decomposition.csv` |
| `outputs.early_stop_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/early_stop_summary.csv` |
| `outputs.pcm_view_success_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/pcm_view_success_summary.csv` |
| `outputs.paired_pcm_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/paired_pcm_summary.csv` |
| `outputs.restart_basin_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/restart_basin_summary.csv` |
| `outputs.update_rule_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/update_rule_summary.csv` |
| `outputs.damping_scaling_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/damping_scaling_summary.csv` |
| `outputs.logical_margin_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/logical_margin_summary.csv` |
| `outputs.prior_sensitivity_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/prior_sensitivity_summary.csv` |
| `outputs.targeted_row_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/targeted_row_summary.csv` |
| `outputs.paired_trapping_subgraph_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/paired_trapping_subgraph_summary.csv` |
| `outputs.mechanism_frequency_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/mechanism_frequency.csv` |
| `outputs.mechanism_overlap_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/mechanism_overlap.csv` |
| `outputs.mechanism_conditional_repair_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/mechanism_conditional_repair.csv` |
| `outputs.mechanism_aware_policy_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/mechanism_aware_policy.csv` |
| `outputs.targeted_repair_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/targeted_repair_summary.csv` |
| `outputs.final_benchmark_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/final_benchmark.csv` |
| `outputs.final_paired_bootstrap_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/final_paired_bootstrap.csv` |
| `outputs.experiment_status_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats/experiment_status.csv` |
| `outputs.experiment_4_paired_pcm_differences_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_4_paired_pcm_differences.csv` |
| `outputs.experiment_9_prior_sensitivity_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_9_prior_sensitivity.csv` |
| `outputs.experiment_10_targeted_rows_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_10_targeted_rows.csv` |
| `outputs.experiment_5_restarts_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_5_restarts.csv` |
| `outputs.experiment_8_logical_costs_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_8_logical_costs.csv` |
| `outputs.experiment_6_update_rules_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_6_update_rules.csv` |
| `outputs.experiment_4_pcm_pairs_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_4_pcm_pairs.csv` |
| `outputs.experiment_5_restart_case_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_5_restart_case_summary.csv` |
| `outputs.experiment_12_mechanism_matrix_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_12_mechanism_matrix.csv` |
| `outputs.experiment_3_pcm_case_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_3_pcm_case_summary.csv` |
| `outputs.experiment_11_paired_subgraph_differences_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_11_paired_subgraph_differences.csv` |
| `outputs.experiment_3_pcm_views_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_3_pcm_views.csv` |
| `outputs.experiment_11_trapping_subgraphs_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_11_trapping_subgraphs.csv` |
| `outputs.final_benchmark_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/final_benchmark_trials.csv` |
| `outputs.experiment_7_repair_classes_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_7_repair_classes.csv` |
| `outputs.targeted_repair_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/targeted_repair_trials.csv` |
| `outputs.experiment_7_damping_scaling_csv` | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace/experiment_7_damping_scaling.csv` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- p_values, validation_p_values: each item is in (0, 1).
- Expected relationship: minimum_shots_per_point <= maximum_shots_per_point.
- Expected relationship: continued_iterations > bp_standard_iterations.
- pcm_views: positive.
- Expected relationship: selected_pcm_views >= pcm_views.
- sampling_batch_size: positive.
- target_gap_cases_per_point: positive.
- analysis_gap_limit_per_point: positive.
- bootstrap_repetitions: positive.
- restart_count: positive.
- ambiguity_top_k: positive.
- restart_message_sigma: non-negative.
- restart_prior_sigma: non-negative.
- prior_variable_sigma: non-negative.
- nms_alpha_values, scaling_values: each item is positive.
- oms_beta_values: each item is non-negative.
- configured_full_osd_order_by_distance: keys exactly {'5', '9'}.
- values of configured_full_osd_order_by_distance: each item is at least 2.
- tesseract_beam: positive.
- tesseract_pq_limit: positive.
- tesseract_timing_shots_per_batch: positive.
- final_case_limit_per_point: -1 or positive.
- development_fraction: in (0, 1).

Environment variables read literally in these cells: `BPOSD_GAP_MECHANISM_REBUILD`, `BPOSD_GAP_MECHANISM_SMOKE`. Initialization code defines their precedence and parsing.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

`algorithm_versions.osd` is computed as `cpp_cs_v1`. Candidate costs use original channel priors even when an experiment deliberately supplies a different ordering score. For `k = n - rank(H)` and order `q`, the tested candidate budget is `1` for `q=0`, otherwise `1 + k + choose(min(q,k),2)`. Negative orders, non-finite scores and priors outside `(0,1)` are rejected. The postprocessor supports the benchmark `osd_cs` method; unsupported methods are rejected rather than silently mapped to another search.
