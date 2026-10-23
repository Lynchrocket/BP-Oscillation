# BP trajectories as OSD reliability information

Test whether temporal BP reliability, instability penalties, conflict or sign consistency improve OSD variable ranking and snapshot selection.

## Directory and entry

```text
01_reliability_dynamics/
  README.md
  notebooks/reliability_dynamics.ipynb
  config/
```

Run [reliability_dynamics.ipynb](notebooks/reliability_dynamics.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

This stage generates its own sampled populations; its path manifest has no external dataset input. Related earlier stages provide scientific motivation, not hidden in-memory variables.

1. Construct native decoder contexts and fixed-budget trajectories.
2. Screen event-driven cohorts and measure error-variable ranks, rank changes and reliability/instability associations.
3. Compare temporal mean/median, sign consistency, SNR, conflict and combined reliability; select parameters on development data.
4. Evaluate handoff/snapshot rules, BP parameter scans and the unbiased BPOSD–Tesseract gap benchmark.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

BP can produce useful soft information while its hard decision is still wrong. Conversely, a stable or large-magnitude LLR can put a required repair variable in an unfavorable OSD position. This stage asks which parts of the trajectory improve reliability ordering, rather than assuming that BP convergence is the best handoff signal.

| Experiment block | Comparison | Interpretation |
| --- | --- | --- |
| Snapshots and rank dynamics | Run the configured OSD variants on multiple snapshots of the same trajectory | Separate soft-information quality from hard-decision correctness |
| Relevant-variable analysis | Compare ranks and motion for variables involved in the diagnostic repair | Test whether instability exposes misleading confidence |
| Dynamics-aware reliability | Compare temporal aggregation, consistency/conflict measures and baseline absolute LLR | Determine whether trajectory features improve the OSD input |
| Handoff and snapshot choice | Compare runtime rules with truth-based best-snapshot ceilings | Separate an available good snapshot from the ability to choose it |
| BP parameter scan | Evaluate BP as a soft-information generator for OSD | Avoid optimizing only stand-alone BP convergence |
| Unbiased Monte Carlo | Freeze a method and measure paired LER, invocation and cost | Check whether conditional gains survive the full population |

The event-driven cohort deliberately searches for standard-failure/forced-OSD-success examples. Its target count, achieved count and fixed-weight selection are part of the result; it is not an unbiased sample of physical noise. Keep distance-5 and distance-9 event counts separate. Snapshot and reliability comparisons must also retain the same OSD order/search budget, since stronger search can mask a worse ordering.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

## Current results and conclusions

The saved screening status reports insufficient gap events, with especially sparse evidence for the larger-distance study. The benchmark tables contain settings without target directional gaps, so gap-closure metrics there are not informative.

The results justify studying trajectory-dependent reliability rankings but do not establish robust improvement across distances. Historical ordering comparisons must be regenerated using the corrected benchmark-compatible OSD semantics.

Evidence: [experiment_status.csv](../../../../data/osd/osd_reliability/01_reliability_dynamics/main/stats/experiment_status.csv), [reliability_method_osd_benchmark.csv](../../../../data/osd/osd_reliability/01_reliability_dynamics/main/stats/reliability_method_osd_benchmark.csv), [bposd_tesseract_gap_closure.csv](../../../../data/osd/osd_reliability/01_reliability_dynamics/main/stats/bposd_tesseract_gap_closure.csv), [selected_reliability_parameters.csv](../../../../data/osd/osd_reliability/01_reliability_dynamics/main/stats/selected_reliability_parameters.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/stats` |
| trace | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/trace` |
| config | `experiment/data/osd/osd_reliability/01_reliability_dynamics/principle_fixes_v1/config` |
| plots | `experiment/research/osd/stages/01_reliability_dynamics/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

Selected reliability parameters are development choices, not standalone test evidence. Conditional OSD success on screened events is not overall LER; inspect denominator and event status for each row.

The rank follow-up tests repair-variable coverage, oracle ceilings, and whether practical corrections retain the apparent opportunity.

Return to [the research-theme workflow](../../README.md).

## JSON execution switches

Execution switches are read from the stage's JSON configuration; environment variables do not override them.

| Key | Default | Values and behavior |
| --- | --- | --- |
| `smoke_mode` | `false` | `false`: use the configured experiment budget. `true`: apply the notebook's existing reduced sampling, iteration and analysis settings and redirect generated outputs to a temporary directory. |
| `reuse_existing_data` | `true` | `true`: reuse available compatible cached results. `false`: recompute the corresponding results. Existing backend and algorithm version checks still apply. |

Set these booleans before running the notebook's configuration cell. Smoke mode does not implicitly disable cache reuse. Temporary smoke outputs are isolated from the configured research outputs.

Configuration: [bp_message_dynamics_for_osd_reliability_experiments.json](config/bp_message_dynamics_for_osd_reliability_experiments.json).
