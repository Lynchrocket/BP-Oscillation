# Configuration for quotient

This is the local configuration reference for [the experiment workflow and conclusions](../README.md). Use conda `qec` and working directory `Decoders/`. Paths and scientific settings have different roles; loading a file does not itself run an experiment.

## Which files are used

| JSON file | Selection in the current notebooks |
| --- | --- |
| [quotient_bp_paths.json](quotient_bp_paths.json) | Named directly by the current notebook code. |
| [quotient_split_paths.json](quotient_split_paths.json) | Named directly by the current notebook code. |

## Path-manifest contract

`schema_version` is integer **1**. `inputs`, `outputs` and `templates` are required mappings from non-empty keys to non-empty repository-relative paths. Absolute paths, parent traversal and paths escaping `Decoders/` are rejected. A named input must exist when read. Template arguments must exactly match the selected placeholder names; there is no automatic version selection or missing-input fallback.

The output meanings are `run` = main dataset root; `stats` = aggregate tables; `trace` = per-case/iteration records and caches; `config` = effective settings and version metadata; `plots` = local figure destination. A path reservation does not mean a result has been generated. Directory names in historical shared-data producers remain for compatibility.

## quotient_bp_paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `outputs.run` | `experiment/data/quotient_decoding/quotient_bp/root/main` |
| `outputs.stats` | `experiment/data/quotient_decoding/quotient_bp/root/main/stats` |
| `outputs.trace` | `experiment/data/quotient_decoding/quotient_bp/root/main/trace` |
| `outputs.config` | `experiment/data/quotient_decoding/quotient_bp/root/main/config` |
| `outputs.plots` | `experiment/research/quotient` |
| `templates` | `{}` (no named entries) |

## quotient_split_paths.json

This file selects locations; it does not set decoder/scientific parameters. All current path keys are listed below.

| Section/key | Current repository-relative value |
| --- | --- |
| `schema_version` | `1` |
| `inputs.error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `inputs.weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `inputs.weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `outputs.run` | `experiment/data/quotient_decoding/quotient_split/root/main` |
| `outputs.stats` | `experiment/data/quotient_decoding/quotient_split/root/main/stats` |
| `outputs.trace` | `experiment/data/quotient_decoding/quotient_split/root/main/trace` |
| `outputs.config` | `experiment/data/quotient_decoding/quotient_split/root/main/config` |
| `outputs.plots` | `experiment/research/quotient` |
| `templates` | `{}` (no named entries) |

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.
