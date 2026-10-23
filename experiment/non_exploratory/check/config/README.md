# Check path manifests

Each notebook selects its own manifest: `adaptive.ipynb` uses [adaptive_bp_repair_paths.json](adaptive_bp_repair_paths.json), `matching_graph.ipynb` uses [matching_graph_paths.json](matching_graph_paths.json), and `check_noise_model.ipynb` uses [noise_model_sanity_check_paths.json](noise_model_sanity_check_paths.json). Scientific parameters and output interpretation are described in [the check README](../README.md). These manifests provide paths only; there is no `print_level` key for these checks.

All path manifests contain the following required fields; no field has an implicit fallback:

| Field | Current/default value | Meaning and constraints |
| --- | --- | --- |
| `schema_version` | `1` | Integer schema version; only 1 is accepted. |
| `inputs` | Mapping listed below | Named existing inputs; `{}` when no dataset is read. |
| `outputs` | Mapping listed below | Explicit output files/directories; includes `run`, `stats`, `trace`, `config` and `plots`. |
| `templates` | Mapping listed below | Parameterized output filenames; `{}` when unused. |

Every mapping key and value must be a non-empty string. Values are paths relative to `Decoders/`; absolute paths, parent traversal and paths escaping the project are rejected. Named inputs must exist when read. Outputs are resolved explicitly; these path fields do not create scientific settings or enable sampling. `run` is the default main dataset root, `stats` is for aggregate results, `trace` for per-case/iteration data, `config` for effective settings, and `plots` for local figures. Supplying another run name uses the path helper's run-specific directory; the notebooks shown here use their stated main paths. Do not run concurrent writers into the same result file.

## adaptive_bp_repair_paths.json

### `inputs`

| Key | Repository-relative value |
| --- | --- |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |

### `outputs`

| Key | Repository-relative value |
| --- | --- |
| `run` | `experiment/data/decoder_intervention/adaptive_bp_repair/root/main` |
| `stats` | `experiment/data/decoder_intervention/adaptive_bp_repair/root/main/stats` |
| `trace` | `experiment/data/decoder_intervention/adaptive_bp_repair/root/main/trace` |
| `config` | `experiment/data/decoder_intervention/adaptive_bp_repair/root/main/config` |
| `plots` | `experiment/non_exploratory/check/plots/adaptive_bp_repair/main` |

### `templates`

Empty for this manifest.

## matching_graph_paths.json

### `inputs`

Empty for this manifest.

### `outputs`

| Key | Repository-relative value |
| --- | --- |
| `run` | `experiment/data/graph_analysis/ca_matching_graphs/root/main` |
| `stats` | `experiment/data/graph_analysis/ca_matching_graphs/root/main/stats` |
| `trace` | `experiment/data/graph_analysis/ca_matching_graphs/root/main/trace` |
| `config` | `experiment/data/graph_analysis/ca_matching_graphs/root/main/config` |
| `plots` | `experiment/non_exploratory/check/plots/matching_graph` |
| `circ_code1_svg` | `experiment/non_exploratory/check/plots/matching_graph/mg_test/circ_code1.svg` |
| `mg_code1_svg` | `experiment/non_exploratory/check/plots/matching_graph/mg_test/mg_code1.svg` |
| `circ_code2_svg` | `experiment/non_exploratory/check/plots/matching_graph/mg_test/circ_code2.svg` |
| `mg_code2_svg` | `experiment/non_exploratory/check/plots/matching_graph/mg_test/mg_code2.svg` |
| `circuit_level_round_1_svg` | `experiment/non_exploratory/check/plots/matching_graph/mg_test/circuit-level,round=1.svg` |

### `templates`

Empty for this manifest.

## noise_model_sanity_check_paths.json

### `inputs`

| Key | Repository-relative value |
| --- | --- |
| `reference` | `experiment/data/input_snapshots/noise_model_sanity_check/root/reference` |
| `2026_04_14_data_err_rates_c_ldpc_dec_bposd_dist_5_r_3_shots_1_0e_09_txt` | `experiment/data/input_snapshots/noise_model_sanity_check/root/reference/2026-04-14,data=err_rates,c=ldpc,dec=bposd,dist=5,r=3,shots=1.0e+09.txt` |
| `2026_04_15_data_err_rates_c_ldpc_dec_tesseract_dist_5_r_3_shots_1_0e_09_txt` | `experiment/data/input_snapshots/noise_model_sanity_check/root/reference/2026-04-15,data=err_rates,c=ldpc,dec=tesseract,dist=5,r=3,shots=1.0e+09.txt` |

### `outputs`

| Key | Repository-relative value |
| --- | --- |
| `run` | `experiment/data/noise_validation/noise_model_sanity_check/root/main` |
| `stats` | `experiment/data/noise_validation/noise_model_sanity_check/root/main/stats` |
| `trace` | `experiment/data/noise_validation/noise_model_sanity_check/root/main/trace` |
| `config` | `experiment/data/noise_validation/noise_model_sanity_check/root/main/config` |
| `plots` | `experiment/non_exploratory/check/plots/noise_model_sanity_check/main` |
| `bposd_tesseract_d_5_r_3_csv` | `experiment/data/noise_validation/noise_model_sanity_check/root/main/stats/bposd_tesseract_d=5_r=3.csv` |

### `templates`

Empty for this manifest.

The matching-graph shared-data root retains its historical producer name for compatibility. Active notebook/configuration/figure names use matching graph. Its SVG keys identify the local example code (code1 = 5-bit, code2 = 9-bit) and diagram kind (`circ` = timeline, `mg` = detector matching graph). The two long noise-reference input keys select BPOSD and Tesseract reference tables respectively; the reference directory itself is also named explicitly.
