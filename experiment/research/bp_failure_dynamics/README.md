# BP failure dynamics

Determine what BP instability diagnoses, classify failure trajectories, and investigate why some logically wrong endpoints remain stable.

## Folder structure

```text
bp_failure_dynamics/
  README.md
  stages/
    01_message_dynamics/
      README.md
      notebooks/
      config/
    02_failure_mechanisms/
      README.md
      notebooks/
      config/
    03_logical_degeneracy/
      README.md
      notebooks/
      config/
```

Figures are saved under each stage’s configured local output, while shared data remain under `experiment/data/`. Each stage README identifies its input versions, workflow, actual evidence and limitations.

## Reading and execution order

Message diagnostics → failure taxonomy and controlled repair → logical-sector mechanism study. The final branch investigates stable logical failures; it is not another benchmark that should be pooled with the others.

| Stage | Question and entry |
| --- | --- |
| [01_message_dynamics](stages/01_message_dynamics/README.md) | Determine whether BP message instability diagnoses failure of the syndrome constraint, logical decoding failure, or both. Separate syndrome consistency, logical correctness and exact recovery of the sampled physical error. |
| [02_failure_mechanisms](stages/02_failure_mechanisms/README.md) | Separate dynamical syndrome failures, premature first-valid failures and persistent logical failures, then test which can be repaired by stopping changes or Adaptive BP. |
| [03_logical_degeneracy](stages/03_logical_degeneracy/README.md) | Test whether stable logical failures are associated with competition between logical sectors, and whether targeted perturbations reveal useful sector changes. |

Use conda `qec` with working directory `Decoders/`. Initialize the selected stage independently, confirm its named input datasets, then run acquisition/diagnostics before aggregate tables and plots. Follow actual data dependencies rather than executing all notebooks in one shared namespace. A matching completed cache may be reused where the stage’s own configuration allows it.

## How the scientific argument progresses

Stage 01 tests the original instability explanation and establishes that satisfying constraints and selecting the correct logical sector are different problems. Stage 02 asks which failure trajectories respond to a different stopping time, more computation or a different graph, with equal-budget controls. Stage 03 follows the stable-wrong branch using exact logical-sector ambiguity and perturbation susceptibility.

Read the three conclusions together: residuals diagnose constraint failure more strongly than conditional logical failure; graph diversity repairs many selected failures; logical ambiguity provides a separate explanation worth testing for stable wrong endpoints. None of these observations alone establishes an online logical-error detector or an unbiased LER improvement.

## What the current results support

Saved diagnostics separate syndrome consistency from logical correctness: residuals are stronger indicators of unresolved constraints, while sector ambiguity helps explain selected stable wrong endpoints. Historical redundant-view and Adaptive BP comparisons suggest that graph diversity matters beyond additional plain-BP iterations. The corrected randomized restart comparison has not been regenerated. Consult stage-level evidence and implementation boundaries before combining these observations.

## Configuration and results

Each `config/README.md` documents the local files, values, path contract and code-level constraints. No root-level configuration document is required. `stats/` is for aggregates, `trace/` for per-case/iteration data, and figures are local. Native BP research traces use the existing Cython adapters; historical Python timings and current native timings are not interchangeable. Existing producer IDs and dataset paths have not changed.

Return to [the research index](../README.md).
