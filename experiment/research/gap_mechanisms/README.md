# BPOSD–Tesseract gap mechanisms

Locate a real decoder gap, attribute sensitivity on enriched cases, then test whether PCM candidate generation can be converted into useful non-oracle selection.

## Folder structure

```text
gap_mechanisms/
  README.md
  stages/
    00_gap_exploration/
      README.md
      notebooks/
      config/
    01_mechanism_attribution/
      README.md
      notebooks/
      config/
    02_pcm_fixed_points/
      README.md
      notebooks/
      config/
```

Figures are saved under each stage’s configured local output, while shared data remain under `experiment/data/`. Each stage README identifies its input versions, workflow, actual evidence and limitations.

## Reading and execution order

Gap detection → paired mechanism interventions → PCM generation/selection and frozen final LER. Each stage generates its own populations under its own configuration.

| Stage | Question and entry |
| --- | --- |
| [00_gap_exploration](stages/00_gap_exploration/README.md) | Establish paired decoder outcomes before optimizing stopping, OSD order or handoff. Determine whether the sampled regime actually contains the target BPOSD-fail/Tesseract-success events. |
| [01_mechanism_attribution](stages/01_mechanism_attribution/README.md) | Distinguish premature stopping, PCM sensitivity, restart basins, update-rule sensitivity and logical ambiguity using paired interventions on actual gap cases. |
| [02_pcm_fixed_points](stages/02_pcm_fixed_points/README.md) | Test whether equivalent PCM changes expose alternative BP fixed points, how this relates to logical near-degeneracy, and whether non-oracle selection can turn generated candidates into lower LER. |

Use conda `qec` with working directory `Decoders/`. Initialize the selected stage independently, confirm its named input datasets, then run acquisition/diagnostics before aggregate tables and plots. Follow actual data dependencies rather than executing all notebooks in one shared namespace. A matching completed cache may be reused where the stage’s own configuration allows it.

## How the scientific argument progresses

Stage 00 screens the early-acceptance/OSD-handoff hypothesis on paired shots. Its current low-gap population is a screening result, not evidence that every later conditional mechanism is absent. Stage 01 deliberately enriches disagreements and probes stopping, graph representation, update dynamics and alternative basins; mechanism labels can overlap. Stage 02 refines the graph explanation by separating basis changes from redundancy and candidate generation from candidate selection.

The final scientific test is a frozen decoder on an unbiased population. Conditional repair opportunity, oracle ensemble coverage and final LER are distinct quantities. A favorable result in an enriched or oracle study cannot replace the final repair/regression comparison, and results from differently sampled stages should not be combined as one denominator.

## What the current results support

The low-noise screening population lacks target gap events, while enriched historical populations reveal PCM and trajectory sensitivity. Partial PCM studies expose a difference between candidate-generation opportunity and practical selection. The final-policy evaluation remains incomplete, and the corrected OSD and message-state definitions have no regenerated evidence.

## Configuration and results

Each `config/README.md` documents the local files, values, path contract and code-level constraints. No root-level configuration document is required. `stats/` is for aggregates, `trace/` for per-case/iteration data, and figures are local. Native BP research traces use the existing Cython adapters; historical Python timings and current native timings are not interchangeable. Existing producer IDs and dataset paths have not changed.

Return to [the research index](../README.md).
