# Configuration for 03 logical degeneracy

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [instability_degeneracy_experiments.json](instability_degeneracy_experiments.json) | Named directly by the current notebook code. |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## instability_degeneracy_experiments.json

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
| `stable_matched_pairs` | `100` | Stage-specific stable matched pairs; used in the notebook's corresponding analysis. |
| `unstable_matched_pairs` | `100` | Stage-specific unstable matched pairs; used in the notebook's corresponding analysis. |
| `exact_enumeration_chunk_size` | `8192` | Stage-specific exact enumeration chunk size; used in the notebook's corresponding analysis. |
| `exact_max_nullity` | `20` | Stage-specific exact max nullity; used in the notebook's corresponding analysis. |
| `susceptibility_k_values` | `[4,8,12,16]` | Stage-specific susceptibility k values; used in the notebook's corresponding analysis. |
| `perturbation_eta_values` | `[0.01,0.02,0.05,0.1]` | Stage-specific perturbation eta values; used in the notebook's corresponding analysis. |
| `num_perturbations` | `50` | Stage-specific num perturbations; used in the notebook's corresponding analysis. |
| `random_subset_replicates` | `5` | Stage-specific random subset replicates; used in the notebook's corresponding analysis. |
| `probability_clip` | `1e-12` | Stage-specific probability clip; used in the notebook's corresponding analysis. |
| `cv_folds` | `5` | Cross-validation split setting; preserve grouping and development/test separation. |
| `bootstrap_repetitions` | `2000` | Bootstrap resampling setting for the reported uncertainty calculation. |
| `random_seed` | `20260907` | Seed for the configured random support/decoder experiment; retain it with the result. |

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `outputs.run` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main` |
| `outputs.stats` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats` |
| `outputs.trace` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/trace` |
| `outputs.config` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/config` |
| `outputs.plots` | `experiment/research/bp_failure_dynamics/stages/03_logical_degeneracy/plots/main` |
| `outputs.instability_degeneracy_input_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/instability_degeneracy_input_audit.csv` |
| `outputs.instability_degeneracy_population_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/trace/instability_degeneracy_population_cases.csv` |
| `outputs.instability_degeneracy_population_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/instability_degeneracy_population_summary.csv` |
| `outputs.instability_degeneracy_stability_sensitivity_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/instability_degeneracy_stability_sensitivity.csv` |
| `outputs.logical_ambiguity_selection_audit_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_selection_audit.csv` |
| `outputs.exact_logical_mass_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/trace/exact_logical_mass_cases.csv` |
| `outputs.exact_logical_mass_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/exact_logical_mass_summary.csv` |
| `outputs.logical_ambiguity_group_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_group_summary.csv` |
| `outputs.logical_ambiguity_paired_comparison_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_paired_comparison.csv` |
| `outputs.logical_ambiguity_auc_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_auc.csv` |
| `outputs.logical_ambiguity_conditional_curves_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_conditional_curves.csv` |
| `outputs.logical_perturbation_trials_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/trace/logical_perturbation_trials.csv` |
| `outputs.logical_susceptibility_cases_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/trace/logical_susceptibility_cases.csv` |
| `outputs.logical_susceptibility_summary_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_susceptibility_summary.csv` |
| `outputs.logical_susceptibility_group_comparison_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_susceptibility_group_comparison.csv` |
| `outputs.susceptibility_vs_logical_gap_csv` | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/susceptibility_vs_logical_gap.csv` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- Expected relationship: long_iterations >= benchmark_iterations >= diagnostic_iterations >= 3.
- ms_scaling_factor: finite and positive.
- late_fraction: in (0, 1].
- message_zero_tolerance: non-negative.
- primary_stability_threshold: chosen from stability_thresholds.
- stability_thresholds: each item is finite and positive.
- stable_matched_pairs: positive.
- unstable_matched_pairs: positive.
- exact_enumeration_chunk_size: positive.
- exact_max_nullity: non-negative.
- susceptibility_k_values: each item is positive.
- perturbation_eta_values: each item is finite and positive.
- Expected relationship: num_perturbations >= random_subset_replicates > 1.
- probability_clip: in (0, 0.5).
- cv_folds: at least 2.
- bootstrap_repetitions: positive.

Native BP results use the recorded backend marker. Non-empty historical outputs without a compatible marker are not silently mixed with native results. Reuse follows the configured cache flag, file availability and any backend/algorithm version checks; configuration hashes do not reject caches. Timing remains backend-specific.
