# Enumerate and test trapping-set candidates

Turn empirical graph observations into explicit candidate structures and test their decoder behavior using recorded error supports.

## Directory and entry

```text
02_candidates/
  README.md
  notebooks/candidates.ipynb
  config/
  plots/                    # generated when configured figures are saved
```

Run [candidates.ipynb](notebooks/candidates.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `failure_cases_with_errors` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors` |
| `ca_bfts_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bfts_cases.csv` |
| `ca_bstf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bstf_cases.csv` |
| `ca_bs_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bs_cases.csv` |
| `ca_bf_cases_csv` | `experiment/data/failure_analysis/failure_regions/root/main/trace/circuit_level/cases/circuit_level-1,perfect meas,shots=1000000/with_errors/ca_bf_cases.csv` |

1. Read the with-errors failure-case inputs from the producer path in the manifest.
2. Construct the graph/PCM representations required by classical and stabilizer-symmetry candidate definitions.
3. Enumerate, label and inspect candidate supports, then evaluate the notebook’s selected decoding examples.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

Empirical motifs become useful hypotheses only when their defining supports can be constructed and tested. This stage examines classical parity-graph candidates and stabilizer-symmetry candidates against the recorded error-support population. Candidate membership is a structural property; failure is an outcome of a particular decoder, prior and input.

Read enumeration together with the selected decoding examples. A candidate that does not fail under those settings is a counterexample to an unconditional trapping claim, not necessarily an enumeration error. Conversely, a failing example does not show that every embedding or every prior with the same motif fails. Retain candidate definition, support, syndrome and decoder setting when comparing with the empirical motifs from stage 01.

The current stage provides candidate analysis rather than a completed population-level validation. The absence of current aggregate results constrains the conclusions that can be drawn from the implemented enumeration.

## Current results and conclusions

该部分实验未进行，可参考data_packedup自行分析

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/stats` |
| trace | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/trace` |
| config | `experiment/data/graph_analysis/trapping_sets/02_candidates/main/config` |
| plots | `experiment/research/trapping_sets/stages/02_candidates/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

The stage name describes candidate analysis; enumeration alone is not evidence that every candidate traps the configured decoder. Keep classical parity-graph candidates distinct from stabilizer/logical effects.

Compare candidate-level evidence with empirical stage-01 motifs and explicit expansion/contraction tests.

Return to [the research-theme workflow](../../README.md).
