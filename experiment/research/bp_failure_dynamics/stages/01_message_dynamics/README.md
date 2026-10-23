# BP message dynamics and the limits of stability

Determine whether BP message instability diagnoses failure of the syndrome constraint, logical decoding failure, or both. Separate syndrome consistency, logical correctness and exact recovery of the sampled physical error.

## Directory and entry

```text
01_message_dynamics/
  README.md
  notebooks/message_dynamics.ipynb
  config/
```

Run [message_dynamics.ipynb](notebooks/message_dynamics.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |

1. Read the shared weight-2/weight-3 cases and reconstruct code matrices, priors and syndromes.
2. Validate the instrumented native BP trace, then build matched success/failure cohorts and trace message/LLR histories.
3. Measure late residuals, sign flips, period-2 behavior, extrinsic conflict and their localized subgraphs.
4. Form the four syndrome/logical outcome groups; compute conditional enrichment, residual AUROC and exact sampled-error-match audits.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

The motivating hypothesis is that persistent, spatially concentrated message motion explains BP failures. The first task is to determine which kind of failure it explains. For correction `c`, syndrome consistency means `H c = s` over GF(2); logical correctness means that the correction has the sampled error's logical observables. Recovering the exact sampled error is a third, stricter event. A small message residual does not imply either of the latter two events.

An edge residual compares the same directed message in consecutive iterations. Prior-normalized residuals account for the input LLR scale; the late window and threshold come from the local configuration. Compare message motion with posterior-LLR motion because incoming messages can cancel: stable beliefs can coexist with moving individual messages. Sign flips, period-2 behavior and the selected unstable subgraph describe different aspects of that motion.

| Analysis | Comparison and control | Scientific interpretation |
| --- | --- | --- |
| Representative trajectories | Message histories and posterior beliefs for the same decoder input | Identify oscillation or cancellation without treating a single case as population evidence |
| Matched outcomes | Success and failure cases with controlled distance, prior and error-weight cohort | Test whether instability is enriched in failure rather than merely in harder inputs |
| Four outcome groups | Cross syndrome consistency with logical correctness | Separate failure to satisfy constraints from selection of a wrong logical sector |
| Conditional prediction | Residual AUROC for inconsistency, then for logical error within consistent cases | Test the stronger claim that instability diagnoses logical failure after constraints are satisfied |
| Exact-error audit | Compare correction bits with sampled error bits separately from logical labels | Avoid counting a logically equivalent correction as a decoding failure |

A fixed-budget diagnostic trace may continue after the first syndrome-consistent iterate; the benchmark decoder may stop there. Keep those endpoints separate. Repeated physical supports can share a syndrome and logical label, so preserve the notebook's deduplication and cohort denominators when interpreting the aggregate tables.

## Current results and conclusions

Saved residual-prediction tables show that late message motion distinguishes syndrome inconsistency more clearly than logical error among syndrome-consistent corrections. The stable-endpoint summary also contains syndrome-consistent but logically wrong cases. Convergence is therefore insufficient evidence of logical correctness.

The exact-error audit separates recovery of the sampled physical support from recovery of its logical sector. These conclusions concern the saved fixed-weight population; they do not establish an online logical-error detector or a physical-noise LER estimate.

Evidence: [residual_prediction_auc.csv](../../../../data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/residual_prediction_auc.csv), [stable_valid_logical_wrong_fraction.csv](../../../../data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/stable_valid_logical_wrong_fraction.csv), [exact_error_match_summary.csv](../../../../data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats/exact_error_match_summary.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/stats` |
| trace | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/trace` |
| config | `experiment/data/bp_dynamics/bp_failure_dynamics/01_message_dynamics/main/config` |
| plots | `experiment/research/bp_failure_dynamics/stages/01_message_dynamics/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

These are diagnostic results on selected fixed-weight cases, not physical-noise LER estimates or a deployable logical-error detector. Exact physical-error recovery is stricter than logical success. Native timings and historical Python timings have different meanings.

Use these distinctions when reading stage 02; its taxonomy must not equate convergence with success.

Return to [the research-theme workflow](../../README.md).
