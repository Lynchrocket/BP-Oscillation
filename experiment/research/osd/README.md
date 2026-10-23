# OSD reliability and handoff studies

Study two distinct decisions: how BP trajectories should rank OSD repair variables, and when/how extra OSD work should be invoked.

## Folder structure

```text
osd/
  README.md
  stages/
    01_reliability_dynamics/
      README.md
      notebooks/
      config/
    02_rank_followup/
      README.md
      notebooks/
      config/
    01_online_handoff/
      README.md
      notebooks/
      config/
    02_intervention_gating/
      README.md
      notebooks/
      config/
```

Figures are saved under each stage’s configured local output, while shared data remain under `experiment/data/`. Each stage README identifies its input versions, workflow, actual evidence and limitations.

## Reading and execution order

Branch A: reliability dynamics → repair-rank/snapshot follow-up. Branch B: causal online handoff → intervention gate and conditional snapshot/order selection. Numeric stage prefixes belong to their branches, not to one linear four-stage execution chain.

| Stage | Question and entry |
| --- | --- |
| [01_reliability_dynamics](stages/01_reliability_dynamics/README.md) | Test whether temporal BP reliability, instability penalties, conflict or sign consistency improve OSD variable ranking and snapshot selection. |
| [02_rank_followup](stages/02_rank_followup/README.md) | Explain which repair variables are misplaced in the reliability ordering and determine whether practical rank/snapshot rules can capture the oracle opportunity. |
| [01_online_handoff](stages/01_online_handoff/README.md) | Study when BP first exposes useful soft information and whether an online policy can retain that information without seeing future iterations or true logical outcomes. |
| [02_intervention_gating](stages/02_intervention_gating/README.md) | Separate the decision to invoke extra OSD work from the conditional decision about snapshot and order, under calibrated cost and held-out evaluation. |

Use conda `qec` with working directory `Decoders/`. Initialize the selected stage independently, confirm its named input datasets, then run acquisition/diagnostics before aggregate tables and plots. Follow actual data dependencies rather than executing all notebooks in one shared namespace. A matching completed cache may be reused where the stage’s own configuration allows it.

## How the two research branches connect

The reliability branch first asks whether BP trajectories contain better OSD orderings than their final hard decision suggests, then studies repair-variable ranks, oracle ceilings and practical ranking rules. The handoff branch asks when a causal policy can use that information, then separates permission to intervene from conditional snapshot/order selection. These branches share a question about useful soft information, not a four-stage mandatory execution chain.

Across both branches, distinguish an endpoint forced-OSD baseline, a fixed-budget snapshot, a truth-based oracle and a causal policy. Search order and candidate budget remain part of each comparison. A high oracle ceiling motivates selection research; it does not validate the existing selector. Final conclusions require enough events at each distance and whole-population evaluation including regressions on baseline successes.

## What the current results support

Historical reliability and snapshot studies expose opportunities that the practical dynamic-rank rule does not reliably exploit. Evidence for the larger-distance setting is sparse. Online-handoff evaluation is incomplete; the intervention-gating experiment has no current results. Historical OSD comparisons cannot validate the corrected implementation.

## Configuration and results

Each `config/README.md` documents the local files, values, path contract and code-level constraints. No root-level configuration document is required. `stats/` is for aggregates, `trace/` for per-case/iteration data, and figures are local. Native BP research traces use the existing Cython adapters; historical Python timings and current native timings are not interchangeable. Existing producer IDs and dataset paths have not changed.

Return to [the research index](../README.md).
