# Configuration for failure regions

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
| `inputs.benchmark` | `experiment/data/input_snapshots/failure_regions/root/benchmark` |
| `inputs.trapping_structure` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis` |
| `inputs.phenomenological_trapping_structure` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis` |
| `inputs.event_graph_features_csv` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/circuit_level/graph_analysis/event_graph_features.csv` |
| `inputs.event_graph_features_csv_7` | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace/phenomenological/graph_analysis/event_graph_features.csv` |
| `inputs.2026_05_28_circuit_level_1_perfect_meas_shots_1000000000_csv` | `experiment/data/input_snapshots/failure_regions/root/benchmark/2026-05-28,circuit_level-1,perfect meas,shots=1000000000.csv` |
| `outputs.run` | `experiment/data/failure_analysis/failure_regions/root/main` |
| `outputs.stats` | `experiment/data/failure_analysis/failure_regions/root/main/stats` |
| `outputs.trace` | `experiment/data/failure_analysis/failure_regions/root/main/trace` |
| `outputs.config` | `experiment/data/failure_analysis/failure_regions/root/main/config` |
| `outputs.plots` | `experiment/research/failure_regions/plots` |
| `outputs.code_capacity` | `experiment/data/failure_analysis/failure_regions/root/main/trace/code_capacity` |
| `outputs.cases` | `experiment/data/failure_analysis/failure_regions/root/main/trace/code_capacity/cases` |
| `outputs.ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/code_capacity/ca_bfts_cases.csv` |
| `outputs.phenomenological` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological` |
| `outputs.cases_1` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases` |
| `outputs.ca_bfts_cases_csv_2` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases/ca_bfts_cases.csv` |
| `outputs.ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases/ca_bs_cases.csv` |
| `outputs.ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases/ca_bstf_cases.csv` |
| `outputs.ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases/ca_bf_cases.csv` |
| `outputs.ca_bfts_cases_csv_3` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases/ca_bfts_cases.csv` |
| `outputs.circuit_level` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level` |
| `outputs.cases_4` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases` |
| `outputs.ca_bfts_cases_csv_5` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bfts_cases.csv` |
| `outputs.ca_bs_cases_csv_6` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bs_cases.csv` |
| `outputs.ca_bstf_cases_csv_7` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bstf_cases.csv` |
| `outputs.ca_bf_cases_csv_8` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bf_cases.csv` |
| `outputs.ca_bf_cases_csv_9` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bf_cases.csv` |
| `outputs.ca_bs_cases_csv_10` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bs_cases.csv` |
| `outputs.ca_bfts_cases_csv_11` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bfts_cases.csv` |
| `outputs.ca_bstf_cases_csv_12` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bstf_cases.csv` |
| `outputs.belief_track` | `experiment/research/failure_regions/plots/migrated/trace/circuit_level/belief_track` |
| `outputs.belief_history_01_txt` | `experiment/research/failure_regions/plots/migrated/trace/circuit_level/belief_track/belief_history_01.txt` |
| `outputs.hard_decision_evolution_gif` | `experiment/research/failure_regions/plots/migrated/trace/circuit_level/belief_track/hard_decision_evolution.gif` |
| `outputs.ca_bfts_cases_csv_1` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bfts_cases.csv` |
| `outputs.graph_vis` | `experiment/research/failure_regions/plots/graph_vis` |
| `outputs.main_cases` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases` |
| `outputs.dem_terms_expanded_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/code_capacity/dem_structure_analysis/dem_terms_expanded.csv` |
| `outputs.dem_node_degrees_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/code_capacity/dem_structure_analysis/dem_node_degrees.csv` |
| `outputs.cases_ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/code_capacity/cases/ca_bfts_cases.csv` |
| `outputs.bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/circuit_level-1,perfect meas,shots=1000000000/bfts_cases.csv` |
| `outputs.bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/circuit_level-1,perfect meas,shots=1000000000/bs_cases.csv` |
| `outputs.bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/circuit_level-1,perfect meas,shots=1000000000/bf_cases.csv` |
| `outputs.bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/circuit_level-1,perfect meas,shots=1000000000/bstf_cases.csv` |
| `outputs.CA_9_2_5,p=0_1,belief_history_01_txt` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/belief_track/CA_9_2_5,p=0.1,belief_history_01.txt` |
| `outputs.belief_track_belief_history_01_txt` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/belief_track/belief_history_01.txt` |
| `outputs.circuit_level-1,perfect meas,shots=1000000_ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/ca_bs_cases.csv` |
| `outputs.circuit_level-1,perfect meas,shots=1000000_ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/ca_bfts_cases.csv` |
| `outputs.circuit_level-1,perfect meas,shots=1000000_ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/ca_bf_cases.csv` |
| `outputs.circuit_level-1,perfect meas,shots=1000000_ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/ca_bstf_cases.csv` |
| `outputs.with_errors_ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bs_cases.csv` |
| `outputs.with_errors_ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bf_cases.csv` |
| `outputs.with_errors_ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bstf_cases.csv` |
| `outputs.ca_sinter_bposd_tesseract_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/phenomenological/benchmark/ca_sinter_bposd_tesseract.csv` |
| `outputs.cnt_cases,phenomenological-3,shots=5000_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/phenomenological/cases/cnt_cases,phenomenological-3,shots=5000.csv` |
| `outputs.cnt_cases,phenomenological-3,shots=1000000_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/phenomenological/cases/cnt_cases,phenomenological-3,shots=1000000.csv` |
| `outputs.dem_structure_summary_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/code_capacity/dem_structure_analysis/dem_structure_summary.csv` |
| `outputs.benchmark_ca_sinter_bposd_tesseract_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/circuit_level/benchmark/ca_sinter_bposd_tesseract.csv` |
| `outputs.cnt_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/main_cases/trace/circuit_level-1,perfect meas,shots=1000000000/cnt_cases.csv` |
| `outputs.circuit_level-1,perfect meas,shots=1000000000_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/circuit_level/benchmark/circuit_level-1,perfect meas,shots=1000000000.csv` |
| `outputs.cnt_cases,circuit_level-3,shots=1000000_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/circuit_level/cases/cnt_cases,circuit_level-3,shots=1000000.csv` |
| `outputs.circuit_level-1,perfect meas,shots=1000000_cnt_cases,circuit_level-3,shots=1000000_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/cnt_cases,circuit_level-3,shots=1000000.csv` |
| `templates.ca_sinter_bposd_tesseract_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/{path_component}/benchmark/ca_sinter_bposd_tesseract.csv` |
| `templates.cnt_cases_noise_model_d_rounds_shots_shots_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/phenomenological/cases/cnt_cases,{noise_model}-{d_rounds},shots={shots}.csv` |
| `templates.cnt_cases_noise_model_d_rounds_shots_shots_csv_13` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/cnt_cases,{noise_model}-{d_rounds},shots={shots}.csv` |
| `templates.noise_model_d_rounds_perfect_meas_shots_max_shots_csv` | `experiment/data/failure_analysis/failure_regions/root/main/stats/{path_component}/benchmark/{noise_model}-{d_rounds},perfect meas,shots={max_shots}.csv` |
| `templates.noise_model_d_rounds_perfect_meas_shots_shots` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}` |
| `templates.ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/ca_bfts_cases.csv` |
| `templates.ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/ca_bs_cases.csv` |
| `templates.ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/ca_bstf_cases.csv` |
| `templates.ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/ca_bf_cases.csv` |
| `templates.cnt_cases_noise_model_d_rounds_shots_shots_csv_14` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/cnt_cases,{noise_model}-{d_rounds},shots={shots}.csv` |
| `templates.with_errors` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors` |
| `templates.ca_bfts_cases_csv_15` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bfts_cases.csv` |
| `templates.ca_bs_cases_csv_16` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bs_cases.csv` |
| `templates.ca_bstf_cases_csv_17` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bstf_cases.csv` |
| `templates.ca_bf_cases_csv_18` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bf_cases.csv` |
| `templates.ca_bf_cases_csv_19` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bf_cases.csv` |
| `templates.ca_bs_cases_csv_20` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bs_cases.csv` |
| `templates.ca_bfts_cases_csv_21` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bfts_cases.csv` |
| `templates.ca_bstf_cases_csv_22` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/{noise_model}-{d_rounds},perfect meas,shots={shots}/with_errors/ca_bstf_cases.csv` |
| `templates.parameter_1_p_parameter_2_belief_history_01_txt` | `experiment/research/failure_regions/plots/migrated/trace/circuit_level/belief_track/{parameter_1},p={parameter_2},belief_history_01.txt` |
| `templates.cases` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases` |
| `templates.ca_bfts_cases_csv_1` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases/ca_bfts_cases.csv` |
| `templates.mg_bfts_phenomenological_CA_9_2_5_p_p_det_parameter_2_gif` | `experiment/research/failure_regions/plots/graph_vis/mg,bfts,phenomenological,CA_9_2_5(p={p}),det={parameter_2}.gif` |
| `templates.tg_bfts_phenomenological_CA_9_2_5_p_p_det_parameter_2_gif` | `experiment/research/failure_regions/plots/graph_vis/tg,bfts,phenomenological,CA_9_2_5(p={p}),det={parameter_2}.gif` |
| `templates.ca_bstf_cases_csv_2` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases/ca_bstf_cases.csv` |
| `templates.mg_bstf_phenomenological_CA_9_2_5_p_p_det_parameter_2_gif` | `experiment/research/failure_regions/plots/graph_vis/mg,bstf,phenomenological,CA_9_2_5(p={p}),det={parameter_2}.gif` |
| `templates.tg_bstf_phenomenological_CA_9_2_5_p_p_det_parameter_2_gif` | `experiment/research/failure_regions/plots/graph_vis/tg,bstf,phenomenological,CA_9_2_5(p={p}),det={parameter_2}.gif` |
| `templates.ca_bfts_cases_csv_3` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases/ca_bfts_cases.csv` |
| `templates.ca_bstf_cases_csv_4` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases/ca_bstf_cases.csv` |
| `templates.ca_bs_cases_csv_5` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases/ca_bs_cases.csv` |
| `templates.ca_bf_cases_csv_6` | `experiment/data/failure_analysis/failure_regions/root/main/trace/{noise_model}/cases/ca_bf_cases.csv` |
| `templates.noise_model_d_rounds_perfect_meas_shots_max_shots` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/{noise_model}-{d_rounds},perfect meas,shots={max_shots}` |
| `templates.cnt_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/{noise_model}-{d_rounds},perfect meas,shots={max_shots}/cnt_cases.csv` |
| `templates.bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/{noise_model}-{d_rounds},perfect meas,shots={max_shots}/bfts_cases.csv` |
| `templates.bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/{noise_model}-{d_rounds},perfect meas,shots={max_shots}/bs_cases.csv` |
| `templates.bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/{noise_model}-{d_rounds},perfect meas,shots={max_shots}/bstf_cases.csv` |
| `templates.bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/main_cases/trace/{noise_model}-{d_rounds},perfect meas,shots={max_shots}/bf_cases.csv` |

Placeholder arguments: `d_rounds`, `max_shots`, `noise_model`, `p`, `parameter_1`, `parameter_2`, `path_component`, `shots`. Values come from the corresponding notebook variables. Match each template's exact names rather than passing every placeholder to every call.

## Validation, overrides and result compatibility

Configuration loading and validation follow the notebook/helper code. The constraints below are taken from its initialization comments; they are not a claim that every condition is dynamically asserted. Invalid required paths/schema/template arguments raise errors. Keep iteration budgets, cohort definitions, seeds and calibration/test separation consistent with any resumed result.
