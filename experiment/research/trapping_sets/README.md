# Trapping-set structure studies

Organize Tanner-graph structural experiments around explicit failure cases, candidate constructions and fixed-weight comparisons.

## Folder structure

```text
trapping_sets/
  README.md
  stages/
    01_failure_structure/
      README.md
      notebooks/
      config/
    02_candidates/
      README.md
      notebooks/
      config/
    03_expansion_contraction/
      README.md
      notebooks/
      config/
    04_low_weight_structure/
      README.md
      notebooks/
      config/
```

Figures are saved under each stage’s configured local output, while shared data remain under `experiment/data/`. Each stage README identifies its input versions, workflow, actual evidence and limitations.

## Reading and execution order

Failure-region samples feed empirical failure/candidate stages; the expansion study uses its fixed historical snapshot; low-weight acquisition feeds the fixed-weight structural stage. The display order is a reading order, not a claim that every stage consumes the preceding output.

| Stage | Question and entry |
| --- | --- |
| [01_failure_structure](stages/01_failure_structure/README.md) | Characterize graph structures associated with BPOSD/Tesseract outcome classes and identify empirical motifs worth testing, without treating every graph feature as a causal trapping set. |
| [02_candidates](stages/02_candidates/README.md) | Turn empirical graph observations into explicit candidate structures and test their decoder behavior using recorded error supports. |
| [03_expansion_contraction](stages/03_expansion_contraction/README.md) | Measure how seeded local structures expand/contract and how generated-event decoding changes with depth, weight and physical error rate. |
| [04_low_weight_structure](stages/04_low_weight_structure/README.md) | Compare syndrome-induced and error-induced structures across fixed DEM error weights, and relate them to BPOSD/Tesseract corrections and disagreements. |

Use conda `qec` with working directory `Decoders/`. Initialize the selected stage independently, confirm its named input datasets, then run acquisition/diagnostics before aggregate tables and plots. Follow actual data dependencies rather than executing all notebooks in one shared namespace. A matching completed cache may be reused where the stage’s own configuration allows it.

## How the structural hypotheses are tested

The four stages have complementary roles: discover empirical motifs, construct candidate supports, perturb local structures, and compare fixed-weight examples. The scientific progression is from association to a testable decoding hypothesis, rather than from a graph picture directly to a general trapping-set claim. The input manifests determine the actual data dependencies; this progression does not require a shared notebook kernel.

Interpret structural similarity with its graph representation, embedding, priors and logical labels. A recurring motif can suggest a candidate, but its trapping behavior requires an explicit decoder outcome. Generated-event fractions and fixed-weight fractions also differ from physical-noise LER. Analyses without saved aggregates remain open questions.

## What the current results support

该部分实验未进行，可参考data_packedup自行分析

## Configuration and results

Each `config/README.md` documents the local files, values, path contract and code-level constraints. No root-level configuration document is required. `stats/` is for aggregates, `trace/` for per-case/iteration data, and figures are local. Native BP research traces use the existing Cython adapters; historical Python timings and current native timings are not interchangeable. Existing producer IDs and dataset paths have not changed.

Return to [the research index](../README.md).
