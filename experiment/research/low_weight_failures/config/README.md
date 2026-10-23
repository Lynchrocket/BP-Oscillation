# Configuration for low weight failures

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [paths.json](paths.json) | Named directly by the current notebook code. |
| [prior_scaling.json](prior_scaling.json) | Named directly by the current notebook code. |
| [prior_scaling_paths.json](prior_scaling_paths.json) | Named directly by the current notebook code. |
| [runtime.json](runtime.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.failure_cases` | `experiment/data/input_snapshots/low_weight_failures/root/failure_cases` |
| `inputs.bfts_cases_csv` | `experiment/data/input_snapshots/low_weight_failures/root/failure_cases/bfts_cases.csv` |
| `outputs.run` | `experiment/data/failure_analysis/low_weight_failures/root/main` |
| `outputs.stats` | `experiment/data/failure_analysis/low_weight_failures/root/main/stats` |
| `outputs.trace` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace` |
| `outputs.config` | `experiment/data/failure_analysis/low_weight_failures/root/main/config` |
| `outputs.plots` | `experiment/research/low_weight_failures/plots/main` |
| `outputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `outputs.fixed_weight_discrepancy_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/stats/fixed_weight_discrepancy.csv` |
| `outputs.weight_2_p=np_geomspace(5e-2, 1e-1, 30)_alphas=np_linspace(0_5, 2_0, 10)_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(5e-2, 1e-1, 30)_alphas=np.linspace(0.5, 2.0, 10)_data.csv` |
| `outputs.weight_4_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_4_data.csv` |
| `outputs.weight_2_p=np_geomspace(1e-8, 5e-2, 30)_alpha=0_8_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(1e-8, 5e-2, 30)_alpha=0.8_data.csv` |
| `outputs.weight_2_p=np_geomspace(5e-2, 1e-1, 30)_alpha=1_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(5e-2, 1e-1, 30)_alpha=1.2_data.csv` |
| `outputs.weight_2_p=np_geomspace(5e-2, 1e-1, 30)_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(5e-2, 1e-1, 30)_data.csv` |
| `outputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `outputs.weight_2_tag_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_tag_data.csv` |
| `outputs.weight_2_p=np_geomspace(5e-2, 1e-1, 30)_alpha=3_0_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(5e-2, 1e-1, 30)_alpha=3.0_data.csv` |
| `outputs.weight_2_p=np_geomspace(5e-2, 1e-1, 30)_alphas=np_linspace(0_5, 2_0, 6)_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(5e-2, 1e-1, 30)_alphas=np.linspace(0.5, 2.0, 6)_data.csv` |
| `outputs.weight_2_p=np_geomspace(1e-8, 1e-1, 30)_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(1e-8, 1e-1, 30)_data.csv` |
| `outputs.weight_2_p=np_geomspace(5e-2, 1e-1, 30)_alpha=0_8_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_p=np.geomspace(5e-2, 1e-1, 30)_alpha=0.8_data.csv` |
| `outputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `templates.weight__low_weight_failure_analysis_weight__data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_{low_weight_failure_analysis_weight}_data.csv` |
| `templates.d_parameter_1__weight__error_weight_prob_weight__data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/d={parameter_1}_weight_{error_weight_prob_weight}_data.csv` |
| `templates.weight__error_weight_prob_weight__p_np_geomspace_5e_2_1e_1_30__data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_{error_weight_prob_weight}_p=np.geomspace(5e-2, 1e-1, 30)_data.csv` |
| `templates.weight__error_weight_prob_weight__p_np_geomspace_5e_2_1e_1_30__alpha_alpha__data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_{error_weight_prob_weight}_p=np.geomspace(5e-2, 1e-1, 30)_alpha={alpha}_data.csv` |
| `templates.weight__error_weight_prob_weight__p_np_geomspace_1e_8_5e_2_30__alpha_alpha__data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_{error_weight_prob_weight}_p=np.geomspace(1e-8, 5e-2, 30)_alpha={alpha}_data.csv` |
| `templates.weight__error_weight_prob_weight__p_np_geomspace_5e_2_1e_1_30__alphas_np_linspace_0_5_2_0_10__data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_{error_weight_prob_weight}_p=np.geomspace(5e-2, 1e-1, 30)_alphas=np.linspace(0.5, 2.0, 10)_data.csv` |

Placeholder arguments: `alpha`, `error_weight_prob_weight`, `low_weight_failure_analysis_weight`, `parameter_1`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## prior_scaling.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `schema_version` | `1` | Supported schema/version identifier; preserve the version accepted by the reader. |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |
| `reuse_existing_data` | `true` | Reuse available cached outputs when true; backend/algorithm version checks still apply where implemented. Configuration hashes are not checked. |
| `random_seed` | `20261022` | Seed for the configured random support/decoder experiment; retain it with the result. |
| `cases` | `[{"distance":5,"weight":2,"sampling":"enumerate","sample_count":null},{"distance":9,"weight":2,"sampling":"enumerate","sample_count":null},{"distance":9,"weight":4,"sampling":"sample","sample_count":5000}]` | Distance/weight/sampling cohorts; keep their populations separate. |
| `p_grids` | `{"high_p":{"start":0.05,"stop":0.1,"count":30},"low_p":{"start":1e-08,"stop":0.05,"count":30}}` | Named physical-error-probability grids; the grid name remains part of the result key. |
| `alpha_grid` | `{"start":0.5,"stop":2.0,"count":10,"extra":[0.8,1.0,1.2,3.0]}` | Minimum-sum scaling scan; do not assume monotonic decoder behavior. |
| `low_p_alphas` | `[0.8,1.0]` | Scaling values for the low-p grid, separate from the high-p scan. |
| `max_iter` | `60` | Decoder iteration limit; interpretation follows the selected backend/section. |
| `d_rounds` | `1` | Number of circuit measurement rounds. |
| `p_meas` | `0.0` | Measurement-readout error probability; keep explicit zeros distinct from inherited p. |
| `osd_order` | `10` | OSD search order; compare with the configured OSD method and decoder budget. |
| `osd_method` | `"osd_cs"` | OSD search method string passed to the decoder. |
| `tesseract_beam` | `8` | Tesseract beam limit passed by the helper. |
| `tesseract_pq_limit` | `200000` | Tesseract priority-queue limit passed by the helper. |

## prior_scaling_paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs` | `{}` (no named entries) |
| `outputs.run` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main` |
| `outputs.stats` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/stats` |
| `outputs.trace` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/trace` |
| `outputs.config` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/config` |
| `outputs.plots` | `experiment/research/low_weight_failures/plots/prior_scaling` |
| `outputs.summary_csv` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/stats/parameter_summary.csv` |
| `outputs.boundary_summary_csv` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/stats/critical_alpha_distribution.csv` |
| `templates.supports_csv` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/trace/distance_{distance}/weight_{weight}/supports.csv` |
| `templates.point_cases_csv` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/trace/distance_{distance}/weight_{weight}/{grid}/p_{p_index:02d}_alpha_{alpha_index:02d}.csv` |
| `templates.boundary_cases_csv` | `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/trace/distance_{distance}/weight_{weight}/critical_alpha/{grid}/p_{p_index:02d}.csv` |
| `templates.phase_plot` | `experiment/research/low_weight_failures/plots/prior_scaling/distance_{distance}_weight_{weight}_{grid}.png` |
| `templates.example_plot` | `experiment/research/low_weight_failures/plots/prior_scaling/distance_{distance}_weight_{weight}_example.png` |

Placeholder arguments: `alpha_index`, `distance`, `grid`, `p_index`, `weight`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## runtime.json

These are the values stored in the local JSON file. The notebook reads this configuration directly; where `smoke_mode` is enabled, its documented reduced-budget settings are applied after loading.

| Key | Current value | Role / dependency |
| --- | --- | --- |
| `print_level` | `2` | Logging: 1 silent configured output, 2 key information, 3 detailed/debug output; does not change scientific populations. |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.

- Expected configuration values:
- print_level: integer; 1 = silent, 2 = key information, 3 = debug.
- All alpha values are positive; every p value is between 0 and 0.5.

## Prior-scaling field domains and dependencies

[prior_scaling.json](prior_scaling.json) is read directly by the notebook and validated by [prior_scaling.py](../../../utils/prior_scaling.py). Invalid values assert rather than falling back. Its fields are:

| Key | Default | Allowed values and dependencies |
| --- | --- | --- |
| `schema_version` | `1` | This configuration schema supports only `1`. |
| `print_level` | `2` | `1`: no logging or progress bars; `2`: cohort/parameter progress and final paths; `3`: per-support progress and parameter diagnostics. Scientific displays and figures remain notebook outputs. |
| `reuse_existing_data` | `true` | Boolean. `true` reuses existing completed parameter-point CSVs; `false` recomputes points. Configuration/source hashes are not compared. |
| `random_seed` | `20261022` | Non-negative integer. Controls the sampled support cohort and the Tesseract seed. |
| `cases` | d=5/w=2 census; d=9/w=2 census; d=9/w=4 sample | Non-empty list. Each `(distance, weight)` pair must be unique. |
| `cases[].distance` | `5` or `9` as listed | `5` selects `CA_9_2_5`; `9` selects `TriangularCode(3)`. |
| `cases[].weight` | `2`, `2`, `4` | Positive integer no larger than the DEM error count. |
| `cases[].sampling` | `enumerate`, `enumerate`, `sample` | `enumerate` visits every support; `sample` draws distinct supports uniformly without replacement. |
| `cases[].sample_count` | `null`, `null`, `5000` | Must be `null` for enumeration; a positive integer at most `comb(num_errors, weight)` for sampling. |
| `p_grids` | `high_p` and `low_p` | Exactly these two named geometric grids. Each contains `start`, `stop`, and `count`. |
| `p_grids.high_p` | `start=0.05`, `stop=0.1`, `count=30` | Finite `0 < start < stop < 0.5`; integer `count >= 2`. Uses the full alpha scan. |
| `p_grids.low_p` | `start=1e-8`, `stop=0.05`, `count=30` | Same domain restrictions. Uses `low_p_alphas`. The shared p=0.05 boundary belongs to both named grids and is never pooled across them. |
| `alpha_grid` | `start=0.5`, `stop=2.0`, `count=10`, `extra=[0.8,1.0,1.2,3.0]` | Positive finite bounds with `start < stop`, integer `count >= 2`, and positive finite extra values. The sorted unique union of the linear grid and extras gives 13 default values. |
| `low_p_alphas` | `[0.8,1.0]` | Non-empty, distinct positive finite values. Independent of the high-p alpha grid. |
| `max_iter` | `60` | Positive integer; BP-only and BPOSD iteration limit. |
| `d_rounds` | `1` | Positive integer used by the circuit builder. Changing it changes matrices; set `reuse_existing_data=false` to recompute affected results. |
| `p_meas` | `0.0` | Measurement error probability in `[0,0.5)`; other error rates inherit p. Set `reuse_existing_data=false` to recompute results after changing it. |
| `osd_order` | `10` | Non-negative integer; applies to BPOSD only. |
| `osd_method` | `osd_cs` | Only this method is supported in this experiment. |
| `tesseract_beam` | `8` | Positive integer; Tesseract beam setting. |
| `tesseract_pq_limit` | `200000` | Positive integer; Tesseract queue limit. |

The generated `effective_config.json` contains a `config` object with the effective input fields and a `backends` mapping recording the BP version, accumulation order and BPOSD/Tesseract wrapper names. It contains no new source/binary hashes or fingerprint. These are output metadata, not user overrides. An existing snapshot is retained when reusing data; otherwise the current configuration is recorded. Legacy snapshots with hash fields remain readable without migration or comparison. Point CSVs are replaced atomically after a complete point, allowing interrupted runs to resume completed work. Changing scientific settings requires explicitly choosing whether to reuse or recompute available results.
