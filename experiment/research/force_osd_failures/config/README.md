# Configuration for force osd failures

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [paths.json](paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.benchmark` | `experiment/data/input_snapshots/force_osd_failures/root/benchmark` |
| `inputs.2026_05_28_circuit_level_1_perfect_meas_shots_1000000000_csv` | `experiment/data/input_snapshots/force_osd_failures/root/benchmark/2026-05-28,circuit_level-1,perfect meas,shots=1000000000.csv` |
| `outputs.run` | `experiment/data/failure_analysis/force_osd_failures/root/main` |
| `outputs.stats` | `experiment/data/failure_analysis/force_osd_failures/root/main/stats` |
| `outputs.trace` | `experiment/data/failure_analysis/force_osd_failures/root/main/trace` |
| `outputs.config` | `experiment/data/failure_analysis/force_osd_failures/root/main/config` |
| `outputs.plots` | `experiment/research/force_osd_failures/plots/main` |
| `outputs.force_OSD` | `experiment/data/failure_analysis/force_osd_failures/root/main/trace/force_OSD` |
| `outputs.bfts_cases_csv` | `experiment/data/failure_analysis/force_osd_failures/root/main/trace/force_OSD/bfts_cases.csv` |
| `outputs.bs_cases_csv` | `experiment/data/failure_analysis/force_osd_failures/root/main/trace/force_OSD/bs_cases.csv` |
| `outputs.bstf_cases_csv` | `experiment/data/failure_analysis/force_osd_failures/root/main/trace/force_OSD/bstf_cases.csv` |
| `outputs.bf_cases_csv` | `experiment/data/failure_analysis/force_osd_failures/root/main/trace/force_OSD/bf_cases.csv` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.
