# Configuration for 02 pcm fixed points

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [next_gap_mechanism_pcm_fixed_point_experiments.json](next_gap_mechanism_pcm_fixed_point_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## next_gap_mechanism_pcm_fixed_point_experiments.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `schema_version` | `8` | Supported schema/version identifier; preserve the version accepted by the reader. |
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
| `analysis_gap_limit_per_point` | `200` | Stage-specific analysis gap limit per point; used in the notebook's corresponding analysis. |
| `development_fraction` | `0.5` | Stage-specific development fraction; used in the notebook's corresponding analysis. |
| `bp_standard_iterations` | `60` | Iteration budget/window; retain the ordering relationships listed below. |
| `continued_iterations` | `100` | Iteration budget/window; retain the ordering relationships listed below. |
| `configured_full_osd_order_by_distance` | `{"5":3,"9":3}` | Code-family levels/distances evaluated or selected by this stage. |
| `osd_orders` | `[0,1,2]` | Stage-specific osd orders; used in the notebook's corresponding analysis. |
| `basis_views` | `50` | Stage-specific basis views; used in the notebook's corresponding analysis. |
| `redundant_views` | `50` | Stage-specific redundant views; used in the notebook's corresponding analysis. |
| `redundant_rows_per_view` | `[1,2,3]` | Stage-specific redundant rows per view; used in the notebook's corresponding analysis. |
| `candidate_max_row_weight` | `10` | Error/row/feature weighting setting; keep the named role with the corresponding experiment. |
| `restart_count` | `100` | Stage-specific restart count; used in the notebook's corresponding analysis. |
| `basin_case_limit` | `50` | Stage-specific basin case limit; used in the notebook's corresponding analysis. |
| `restart_message_sigma` | `0.05` | Stage-specific restart message sigma; used in the notebook's corresponding analysis. |
| `fixed_point_tail_window` | `3` | Stage-specific fixed point tail window; used in the notebook's corresponding analysis. |
| `fixed_point_belief_tolerance` | `1e-06` | Stage-specific fixed point belief tolerance; used in the notebook's corresponding analysis. |
| `logical_margin_control_limit_per_group` | `100` | Stage-specific logical margin control limit per group; used in the notebook's corresponding analysis. |
| `logical_margin_quartiles` | `4` | Stage-specific logical margin quartiles; used in the notebook's corresponding analysis. |
| `ambiguity_top_k` | `8` | Stage-specific ambiguity top k; used in the notebook's corresponding analysis. |
| `pcm_budgets` | `[1,2,4,8,16]` | Stage-specific pcm budgets; used in the notebook's corresponding analysis. |
| `targeting_strategies` | `["random","instability","ambiguity","conflict","hybrid"]` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `targeting_constructions` | `["top_m_distinct","progressive","diversity_constrained"]` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `diversity_mu_values` | `[0.0,0.25,0.5,1.0]` | Stage-specific diversity mu values; used in the notebook's corresponding analysis. |
| `selector_eta_values` | `[0.0,0.25,0.5,1.0,2.0,4.0]` | Stage-specific selector eta values; used in the notebook's corresponding analysis. |
| `selector_temperature_values` | `[0.25,0.5,1.0,2.0]` | Stage-specific selector temperature values; used in the notebook's corresponding analysis. |
| `selector_confidence_epsilon` | `1e-09` | Stage-specific selector confidence epsilon; used in the notebook's corresponding analysis. |
| `bootstrap_repetitions` | `1000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `final_pcm_views` | `16` | Stage-specific final pcm views; used in the notebook's corresponding analysis. |
| `final_case_limit_per_point` | `-1` | Stage-specific final case limit per point; used in the notebook's corresponding analysis. |
| `e1_stream_chunk_shots` | `100` | Population/selection budget; a target is not proof that enough eligible events were obtained. |
| `milp_time_limit_seconds` | `60.0` | Stage-specific milp time limit seconds; used in the notebook's corresponding analysis. |
| `random_seed` | `20260913` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |
| `tesseract_timing_shots_per_batch` | `100` | Population/selection budget; a target is not proof that enough eligible events were obtained. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1` |
| `outputs.stats` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats` |
| `outputs.trace` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace` |
| `outputs.config` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/config` |
| `outputs.plots` | `experiment/research/gap_mechanisms/stages/02_pcm_fixed_points/plots/principle_fixes_v1` |
| `outputs.config_json` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/config/config.json` |
| `outputs.population_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/population.csv` |
| `outputs.population_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/population_summary.csv` |
| `outputs.A1_basis_vs_redundant_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/A1_basis_vs_redundant.csv` |
| `outputs.A2_paired_global_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/A2_paired_global.csv` |
| `outputs.A2_paired_local_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/A2_paired_local.csv` |
| `outputs.B1_B2_margin_correlations_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/B1_B2_margin_correlations.csv` |
| `outputs.B1_B2_margin_bins_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/B1_B2_margin_bins.csv` |
| `outputs.B1_cohort_margins_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/B1_cohort_margins.csv` |
| `outputs.B3_basin_probability_shift_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/B3_basin_probability_shift.csv` |
| `outputs.B4_margin_oscillation_correlation_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/B4_margin_oscillation_correlation.csv` |
| `outputs.B4_grouped_oscillation_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/B4_grouped_oscillation.csv` |
| `outputs.C1_C2_distinct_targeted_scaling_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/C1_C2_distinct_targeted_scaling.csv` |
| `outputs.C2_stratified_generation_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/C2_stratified_generation.csv` |
| `outputs.C3_row_coverage_repair_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/C3_row_coverage_repair.csv` |
| `outputs.D1_D2_generation_diversity_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/D1_D2_generation_diversity.csv` |
| `outputs.D2_margin_sector_entropy_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/D2_margin_sector_entropy.csv` |
| `outputs.D3_selector_tuning_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/D3_selector_tuning.csv` |
| `outputs.D3_frozen_selector_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/D3_frozen_selector.csv` |
| `outputs.D5_generation_selection_gap_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/D5_generation_selection_gap.csv` |
| `outputs.D4_selector_margin_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/D4_selector_margin.csv` |
| `outputs.E_frozen_policy_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/E_frozen_policy.csv` |
| `outputs.E1_pcm_ensemble_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/E1_pcm_ensemble_summary.csv` |
| `outputs.E2_qmin_improvement_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/E2_qmin_improvement.csv` |
| `outputs.E3_mechanism_conditional_repair_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/E3_mechanism_conditional_repair.csv` |
| `outputs.E4_final_benchmark_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/E4_final_benchmark.csv` |
| `outputs.E4_paired_bootstrap_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/E4_paired_bootstrap.csv` |
| `outputs.experiment_status_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/experiment_status.csv` |
| `outputs.C1_C2_targeted_candidates_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/C1_C2_targeted_candidates.csv` |
| `outputs.A2_paired_views_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/A2_paired_views.csv` |
| `outputs.D1_D2_candidate_prefixes_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/D1_D2_candidate_prefixes.csv` |
| `outputs.C3_row_candidates_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/C3_row_candidates.csv` |
| `outputs.A1_pcm_views_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/A1_pcm_views.csv` |
| `outputs.B2_B3_restart_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/B2_B3_restart_trials.csv` |
| `outputs.logical_sector_costs_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/logical_sector_costs.csv` |
| `outputs.E1_pcm_ensemble_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/E1_pcm_ensemble_trials.csv` |
| `outputs.logical_margins_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/logical_margins.csv` |
| `outputs.C2_prefix_outcomes_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/C2_prefix_outcomes.csv` |
| `outputs.A2_paired_local_graph_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/A2_paired_local_graph.csv` |
| `outputs.D3_D5_selector_trials_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/D3_D5_selector_trials.csv` |
| `outputs.E1_candidate_details_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/E1_candidate_details.csv` |
| `outputs.B1_B2_analysis_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/B1_B2_analysis.csv` |
| `outputs.B2_B3_basin_summary_csv` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace/B2_B3_basin_summary.csv` |
| `templates.name` | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats/{name}` |

Placeholder arguments: `name`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- p_values, validation_p_values: each item is in (0, 1).
- Expected relationship: 0 < minimum_shots_per_point <= maximum_shots_per_point.
- sampling_batch_size: positive.
- target_gap_cases_per_point: positive.
- Expected relationship: continued_iterations > bp_standard_iterations > 0.
- configured_full_osd_order_by_distance: keys exactly {'5', '9'}.
- values of configured_full_osd_order_by_distance: each item is at least 2.
- basis_views: positive.
- redundant_views: positive.
- analysis_gap_limit_per_point: positive.
- candidate_max_row_weight: positive.
- redundant_rows_per_view: non-empty; each item is positive.
- restart_count: positive.
- basin_case_limit: positive.
- restart_message_sigma: non-negative.
- fixed_point_tail_window: at least 2.
- fixed_point_belief_tolerance: positive.
- logical_margin_control_limit_per_group: positive.
- logical_margin_quartiles: at least 2.
- ambiguity_top_k: positive.
- selector_confidence_epsilon: positive.
- Expected relationship: minimum of pcm_budgets >= 1.
- Expected relationship: final_pcm_views >= maximum of pcm_budgets.
- targeting_strategies: values exactly {'random', 'instability', 'ambiguity', 'conflict',
- targeting_constructions: values exactly {'top_m_distinct', 'progressive',
- diversity_mu_values: each item is non-negative.
- selector_eta_values: each item is non-negative.
- selector_temperature_values: each item is positive.
- bootstrap_repetitions: positive.
- development_fraction: in (0, 1).
- final_case_limit_per_point: -1 or positive.
- e1_stream_chunk_shots: positive.
- milp_time_limit_seconds: positive.
- tesseract_beam: positive.
- tesseract_pq_limit: positive.
- tesseract_timing_shots_per_batch: positive.

Environment variables read literally in these cells: `NEXT_GAP_PCM_REBUILD`, `NEXT_GAP_PCM_SMOKE`. Initialization code defines their precedence and parsing.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.

## Algorithm-version and output compatibility

The current output run is `principle_fixes_v1`; the previous `main` files are retained as historical results. Inputs keep their explicitly selected paths. Before any result is reused, `check_algorithm_cache` checks `config/algorithm_versions.json`. A missing/mismatched marker in a nonempty result directory raises an assertion; it never silently accepts, deletes or overwrites older data. Fresh runs record the version and include the generated `algorithm_versions` object in the effective configuration. This object is implementation-owned, not a user-selectable algorithm switch.

`algorithm_versions.osd` is computed as `cpp_cs_v1`. Candidate costs use original channel priors even when an experiment deliberately supplies a different ordering score. For `k = n - rank(H)` and order `q`, the tested candidate budget is `1` for `q=0`, otherwise `1 + k + choose(min(q,k),2)`. Negative orders, non-finite scores and priors outside `(0,1)` are rejected. The postprocessor supports the benchmark `osd_cs` method; unsupported methods are rejected rather than silently mapped to another search.

| Key | Allowed values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `fixed_point_message_tolerance` | Finite positive number | `1e-6` | Maximum absolute consecutive C2V/V2C message change allowed over the existing tail window. Both directions must pass, alongside the existing belief and hard-decision checks. |
| `fixed_point_state_tolerance` | Finite positive number | `1e-6` | Maximum absolute full terminal-state difference within an observed endpoint cluster. Cluster members must also have identical hard decisions. |
| `algorithm_versions.fixed_point` | `message_state_v1` | Computed | Distinguishes message-state convergence and endpoint clustering from the earlier belief-only/sector-count interpretation. |

`fixed_point_tail_window` must be at least 2. Restart output adds `fixed_point_id` (`-1` for an ineligible endpoint) and `terminal_message_step`. Aggregate `fixed_point_count` counts eligible endpoint clusters, `basin_entropy` uses their restart frequencies, `restart_sector_count` counts logical sectors and `sector_entropy` uses sector frequencies. `P_true` remains conditional on eligible arrivals; read it with `valid_basin_arrival_rate`.
