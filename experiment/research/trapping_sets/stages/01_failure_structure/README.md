# Empirical trapping-set failure structures

Characterize graph structures associated with BPOSD/Tesseract outcome classes and identify empirical motifs worth testing, without treating every graph feature as a causal trapping set.

## Directory and entry

```text
01_failure_structure/
  README.md
  notebooks/failure_structure.ipynb
  config/
  plots/                    # generated when configured figures are saved
```

Run [failure_structure.ipynb](notebooks/failure_structure.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `failure_cases` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases` |
| `ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bfts_cases.csv` |
| `ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bstf_cases.csv` |
| `ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bs_cases.csv` |
| `ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/ca_bf_cases.csv` |

1. Read the four directional outcome populations produced by failure-region sampling.
2. Build detector-induced/Tanner subgraphs, compute structural features and compare distributions, correlations and PCA.
3. Group graphs using WL hashes/isomorphism and inspect feature/GED-near clusters.
4. Form classical/symmetry candidates and test constructed adversarial events with decoders; keep aggregate summaries separate from per-case traces.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

The starting question is whether decoder disagreement groups have recurring local structures. The notebook compares graph features, their correlations and PCA, then examines isomorphism, WL hashes, feature distances and graph-edit proximity. These tools organize candidate motifs at different resolutions; they are not interchangeable tests of a trapping mechanism.

Compare the two disagreement directions with both-success and both-failure controls under the recorded sampling conditions. A motif found only in a selected failure example may reflect error weight, graph size or case selection. The structural analysis therefore motivates explicit candidate and adversarial decoding tests rather than treating a cluster as a proof of causality.

Keep the graph representation with each result: an error-induced subgraph, a syndrome-induced subgraph and the residual between error and correction describe different objects. A shared WL hash alone does not establish graph isomorphism, and graph isomorphism alone does not guarantee identical decoding behavior when priors or embedding differ. Classical parity-graph candidates and stabilizer-symmetry candidates retain separate interpretations.

## Current results and conclusions

该部分实验未进行，可参考data_packedup自行分析

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/stats` |
| trace | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/trace` |
| config | `experiment/data/graph_analysis/trapping_sets/01_failure_structure/main/config` |
| plots | `experiment/research/trapping_sets/stages/01_failure_structure/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

WL hashes and visual similarity are grouping tools, not proofs of isomorphism or decoding failure. Outcome-selected cases are not an unbiased physical-noise distribution.

Candidate enumeration formalizes structures separately; low-weight structure compares controlled support weights.

Return to [the research-theme workflow](../../README.md).
