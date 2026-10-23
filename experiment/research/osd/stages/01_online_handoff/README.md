# Causal BP-to-OSD handoff

Study when BP first exposes useful soft information and whether an online policy can retain that information without seeing future iterations or true logical outcomes.

## Directory and entry

```text
01_online_handoff/
  README.md
  notebooks/online_handoff.ipynb
  config/
```

Run [online_handoff.ipynb](notebooks/online_handoff.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

This stage generates its own sampled populations; its path manifest has no external dataset input. Related earlier stages provide scientific motivation, not hidden in-memory variables.

1. Generate separate mechanism and unbiased populations with explicit cohort labels.
2. Measure q-min trajectories, feature availability, useful-snapshot persistence and endpoint/prefix/full-oracle ceilings.
3. Fit held-out causal handoff policies and temporal-shuffle controls; test prefix buffers and adaptive OSD order.
4. Scan BP soft-generator settings and run the frozen paired accuracy/cost benchmark.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

The central question is whether a causal policy can hand a useful BP snapshot to OSD before the opportunity disappears, while choosing an appropriate search order. Define the BP endpoint as its first syndrome-consistent iterate, or the configured iteration limit if no such iterate occurs. Forced endpoint OSD and OSD after a fixed maximum-length trace are different references and must remain distinct.

The diagnostic `q_min(t)` is the smallest tested OSD order that succeeds at snapshot t. It uses logical truth and labels the opportunity offline. A runtime feature or prefix buffer at time t can use only observations available through t; a full-trajectory oracle can look ahead and is a separate ceiling.

| Experiment block | Controlled comparison | Purpose |
| --- | --- | --- |
| Temporal map and oracle ceilings | Endpoint, prefix and full-trajectory success under tested orders | Locate when lower-order search is sufficient |
| Feature/model evaluation | Causal reliability and dynamics features on held-out inputs | Test whether a runtime policy can recognize that opportunity |
| Temporal shuffle | Disrupt temporal order while retaining snapshot information | Test whether sequence order contributes beyond snapshot quality |
| Prefix buffer and adaptive order | Select among already observed snapshots and configured orders | Trade retained history/search work against repair quality |
| Final paired benchmark | Frozen policy versus audited endpoint/fixed-budget references | Measure LER, regressions and computational cost on the same shots |

Do not use future snapshots, the true label or an oracle's chosen order as online inputs. Conditional failure-cohort ceilings should be reported with event counts, while the final policy must also be evaluated on baseline successes. A small prefix-oracle example demonstrates an opportunity in those cases; it does not establish a calibrated stopping rule or population-level improvement.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

## Current results and conclusions

The saved temporal and prefix-oracle tables contain examples where earlier soft information offers a successful OSD correction that the endpoint misses. This establishes a diagnostic opportunity in a selected historical population, not a causal online policy result.

The policy-fitting, held-out repair/regression and final cost evaluation are incomplete: 该部分实验未进行，可参考data_packedup自行分析

Evidence: [qmin_temporal_summary.csv](../../../../data/osd/osd_handoff/01_online_handoff/main/stats/qmin_temporal_summary.csv), [prefix_oracle_ceiling.csv](../../../../data/osd/osd_handoff/01_online_handoff/main/stats/prefix_oracle_ceiling.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/stats` |
| trace | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/trace` |
| config | `experiment/data/osd/osd_handoff/01_online_handoff/principle_fixes_v1/config` |
| plots | `experiment/research/osd/stages/01_online_handoff/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

An oracle sees outcomes a causal handoff rule cannot. Small or enriched mechanism cohorts cannot establish a physical-noise LER improvement, and early usefulness alone does not guarantee a safe stopping decision.

Intervention gating separates whether extra work is justified from which snapshot/order to use.

Return to [the research-theme workflow](../../README.md).

## JSON execution switches

Execution switches are read from the stage's JSON configuration; environment variables do not override them.

| Key | Default | Values and behavior |
| --- | --- | --- |
| `smoke_mode` | `false` | `false`: use the configured experiment budget. `true`: apply the notebook's existing reduced sampling, iteration and analysis settings and redirect generated outputs to a temporary directory. |
| `reuse_existing_data` | `true` | `true`: reuse available compatible cached results. `false`: recompute the corresponding results. Existing backend and algorithm version checks still apply. |

Set these booleans before running the notebook's configuration cell. Smoke mode does not implicitly disable cache reuse. Temporary smoke outputs are isolated from the configured research outputs.

Configuration: [bp_dynamics_online_osd_handoff_experiments.json](config/bp_dynamics_online_osd_handoff_experiments.json).
