# Two-stage intervention gating and snapshot selection

Separate the decision to invoke extra OSD work from the conditional decision about snapshot and order, under calibrated cost and held-out evaluation.

## Directory and entry

```text
02_intervention_gating/
  README.md
  notebooks/intervention_gating.ipynb
  config/
```

Run [intervention_gating.ipynb](notebooks/intervention_gating.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

This stage generates its own sampled populations; its path manifest has no external dataset input. Related earlier stages provide scientific motivation, not hidden in-memory variables.

1. Reconcile trace/native/reference OSD definitions and build an unbiased intervention population.
2. Fit and calibrate the Stage-1 intervention gate using endpoint-available features and cost-sensitive thresholds.
3. Condition Stage-2 snapshot selection and adaptive order on the gated population.
4. Freeze both stages and evaluate the paired held-out LER, repairs/regressions, invocation rate and runtime.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

A snapshot intervention can repair a baseline failure and also replace a correct baseline output with a wrong one. This stage separates two decisions: whether to override standard BPOSD at all, and which stored snapshot/order to use after intervention is allowed. Optimizing conditional recovery alone ignores the first decision's regression cost.

Before learning either decision, reconcile the forced-OSD references. The audit checks the decoding input, BP endpoint, signed soft information versus ordering magnitudes, configured method/order/search budget and ordering/tie behavior. A baseline mismatch can otherwise appear as an intervention gain without any new decision rule.

| Phase | Input and decision | Evaluation role |
| --- | --- | --- |
| Reference reconciliation | Identical PCM, syndrome, priors and endpoint soft information | Establish comparable decoder semantics before attributing gains |
| Stage-1 gate | Available endpoint/prefix features decide whether to intervene | Balance recoverable failures against harmful overrides |
| Stage-2 recovery | Select from permitted stored snapshots and tested orders | Recover gated cases without lookahead or logical-truth features |
| Calibration and ablation | Development-only thresholds, feature groups and order policies | Identify the source of benefit and freeze the operating point |
| Held-out paired benchmark | Apply the frozen gate and recovery policy to all shots | Measure total LER, repairs, regressions, invocation and cost |

Training labels may use truth offline; runtime features and the evaluated policy cannot. Calibrate cost-sensitive thresholds on development data rather than the final benchmark. Gate precision, recovery success conditional on intervention and whole-population LER have different denominators and must be read together. No saved aggregate currently establishes the performance of this design.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

## Current results and conclusions

该部分实验未进行，可参考data_packedup自行分析

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/stats` |
| trace | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/trace` |
| config | `experiment/data/osd/osd_handoff/02_intervention_gating/principle_fixes_v1/config` |
| plots | `experiment/research/osd/stages/02_intervention_gating/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

Do not infer selector performance from an oracle or from training data, and do not equate low invocation rate with unchanged LER. This stage generates its own population; the schedule follows online handoff, but no in-memory state is shared.

Use the eventual final paired benchmark to decide whether the complete two-stage policy is justified.

Return to [the research-theme workflow](../../README.md).

## JSON execution switches

Execution switches are read from the stage's JSON configuration; environment variables do not override them.

| Key | Default | Values and behavior |
| --- | --- | --- |
| `smoke_mode` | `false` | `false`: use the configured experiment budget. `true`: apply the notebook's existing reduced sampling, iteration and analysis settings and redirect generated outputs to a temporary directory. |
| `reuse_existing_data` | `true` | `true`: reuse available compatible cached results. `false`: recompute the corresponding results. Existing backend and algorithm version checks still apply. |

Set these booleans before running the notebook's configuration cell. Smoke mode does not implicitly disable cache reuse. Temporary smoke outputs are isolated from the configured research outputs.

Configuration: [two_stage_bposd_intervention_snapshot_experiments.json](config/two_stage_bposd_intervention_snapshot_experiments.json).
