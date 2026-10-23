# Trapping structures at fixed error weight

Compare syndrome-induced and error-induced structures across fixed DEM error weights, and relate them to BPOSD/Tesseract corrections and disagreements.

## Directory and entry

```text
04_low_weight_structure/
  README.md
  notebooks/low_weight_structure.ipynb
  config/
  plots/                    # generated when configured figures are saved
```

Run [low_weight_structure.ipynb](notebooks/low_weight_structure.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |
| `weight_4_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_4_data.csv` |

1. Read the shared weight-2/3/4 datasets while preserving bit-string columns.
2. Collect error-induced and syndrome-induced graph features using the existing graph definitions.
3. Compare outcome-group distributions and WL-hash groups separately for each weight.
4. Inspect representative original-error, BPOSD-correction and Tesseract-correction subgraphs.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

Fixed-weight inputs provide a controlled way to compare structures without mixing every physical error weight. The notebook inspects syndrome-induced and error-induced subgraphs, selected WL-hash groups and the differences among the sampled error, BPOSD correction and Tesseract correction. Those objects explain different aspects of the same decoding input.

A syndrome-induced graph describes the observed detectors; an error-induced graph uses the known support and is an offline explanation. Different supports can share a syndrome. Likewise, two different corrections can be logically equivalent, so a structural difference is not itself a logical failure. Interpret each example together with the syndrome-consistency and logical-outcome labels.

Compare motifs within their weight, code and prior conditions before considering cross-weight changes. Selected hash examples are illustrations rather than estimates of motif prevalence, and the fixed-weight support distribution is not the physical-noise distribution. This branch connects the low-weight acquisition to the other trapping-set hypotheses without assuming it consumes stage 03 outputs.

## Current results and conclusions

该部分实验未进行，可参考data_packedup自行分析

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/stats` |
| trace | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/trace` |
| config | `experiment/data/graph_analysis/trapping_sets/04_low_weight_structure/main/config` |
| plots | `experiment/research/trapping_sets/stages/04_low_weight_structure/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

A syndrome-induced graph, an error-induced graph and a correction-induced graph are different objects. Fixed-weight support statistics are not physical-noise LER estimates; weight-2/3/4 populations must remain separately labeled.

Use low-weight acquisition as the producer, and this stage as its structural consumer rather than duplicating data collection.

Return to [the research-theme workflow](../../README.md).
