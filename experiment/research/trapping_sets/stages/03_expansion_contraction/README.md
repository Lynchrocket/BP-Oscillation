# Expansion and contraction of failure subgraphs

Measure how seeded local structures expand/contract and how generated-event decoding changes with depth, weight and physical error rate.

## Directory and entry

```text
03_expansion_contraction/
  README.md
  notebooks/expansion_contraction.ipynb
  config/
  plots/                    # generated when configured figures are saved
```

Run [expansion_contraction.ipynb](notebooks/expansion_contraction.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `failure_cases` | `experiment/data/input_snapshots/subgraph_structure/root/failure_cases` |
| `noise_model_d_rounds_perfect_meas_shots_expansion_contraction_max_shots` | `experiment/data/input_snapshots/subgraph_structure/root/failure_cases/{noise_model}-{d_rounds},perfect meas,shots={expansion_contraction_max_shots}` |
| `bfts_cases_csv` | `experiment/data/input_snapshots/subgraph_structure/root/failure_cases/{noise_model}-{d_rounds},perfect meas,shots={expansion_contraction_max_shots}/bfts_cases.csv` |

1. Read the explicitly named historical failure-case snapshot and select the source seeds.
2. Generate the notebook’s depth-1/depth-2 structures and weight-2/3/4 event variants.
3. Evaluate the configured p sweep and record generated-event outcomes and structural records.
4. Write separate aggregate expansion/contraction and probability summaries from the per-event records.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

Starting from selected failure structures, this stage asks whether nearby constructions preserve the failure and how that behavior depends on neighborhood depth, error-support weight and p. Depth-1/depth-2 and weight-2/3/4 sections are distinct comparisons; do not pool them into a single structural success fraction.

Hold the source seed and construction rule visible when comparing an expanded or contracted event with its source. Changing the support can change both its syndrome and its prior probability. Counts of failing generated events measure behavior under the construction procedure; probability-weighted summaries address a different question and need their corresponding weights and denominators.

The historical seed snapshot is fixed explicitly. This makes the construction study reproducible with respect to its chosen seeds, but does not make those seeds an unbiased population of physical failures. Conclusions about whether a motif persists under modification require the resulting matched records, not the appearance of a larger or smaller plotted subgraph.

## Current results and conclusions

该部分实验未进行，可参考data_packedup自行分析

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/stats` |
| trace | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/trace` |
| config | `experiment/data/graph_analysis/trapping_sets/03_expansion_contraction/main/config` |
| plots | `experiment/research/trapping_sets/stages/03_expansion_contraction/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

Changing the structure can change the represented event and its probability. Compare matched constructions and conditioning; the historical input snapshot is not silently replaced by a newer producer run.

Use saved aggregate records to assess depth/weight effects before asserting a structural mechanism.

Return to [the research-theme workflow](../../README.md).
