# Instability-guided BP

Test whether an instability diagnostic is useful for intervention selection, and whether candidate-level gains survive held-out feedback evaluation.

## Folder structure

```text
instability_guided_bp/
  README.md
  stages/
    01_graph_intervention/
      README.md
      notebooks/
      config/
    02_dynamics_feedback/
      README.md
      notebooks/
      config/
```

Figures are saved under each stage’s configured local output, while shared data remain under `experiment/data/`. Each stage README identifies its input versions, workflow, actual evidence and limitations.

## Reading and execution order

Localize/probe/compare interventions → audit the signal and candidate pool → freeze feedback and evaluate stochastic LER. Stage 02 reads a named historical stage-01 snapshot.

| Stage | Question and entry |
| --- | --- |
| [01_graph_intervention](stages/01_graph_intervention/README.md) | Test whether local message-instability information helps choose redundant checks, damping or scaling that repair BP failures, rather than merely reducing a diagnostic residual. |
| [02_dynamics_feedback](stages/02_dynamics_feedback/README.md) | Determine whether short-probe instability reduction predicts repair and whether a frozen feedback policy improves held-out stochastic decoding. |

Use conda `qec` with working directory `Decoders/`. Initialize the selected stage independently, confirm its named input datasets, then run acquisition/diagnostics before aggregate tables and plots. Follow actual data dependencies rather than executing all notebooks in one shared namespace. A matching completed cache may be reused where the stage’s own configuration allows it.

## How the scientific argument progresses

Stage 01 asks whether spatial instability and short-probe response can locate useful actions even when instability is an imperfect global error detector. Stage 02 tests alternative explanations for that apparent targeting value: score normalization, biased candidate pools, simple stabilization and extra search budget. It then evaluates repeated feedback, stopping and gating as a complete decoder.

The intervention unit is an action on a fixed decoding input; the deployment unit is a whole physical shot. Compare local repair and regression results with the final stochastic benchmark before concluding that feedback helps. In particular, the current final negative comparison remains part of the conclusion even where a selected-cohort targeting ablation is favorable.

## What the current results support

Saved selected-cohort ablations favor some targeted interventions, but instability reduction is only weakly predictive of repair. The final feedback comparison does not establish a LER advantage over BPOSD in the studied regime. Interpret targeting gains together with regressions and whole-population evaluation.

## Configuration and results

Each `config/README.md` documents the local files, values, path contract and code-level constraints. No root-level configuration document is required. `stats/` is for aggregates, `trace/` for per-case/iteration data, and figures are local. Native BP research traces use the existing Cython adapters; historical Python timings and current native timings are not interchangeable. Existing producer IDs and dataset paths have not changed.

Return to [the research index](../README.md).
