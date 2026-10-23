# From probe gain to useful decoding feedback

Determine whether short-probe instability reduction predicts repair and whether a frozen feedback policy improves held-out stochastic decoding.

## Directory and entry

```text
02_dynamics_feedback/
  README.md
  notebooks/dynamics_feedback.ipynb
  config/
```

Run [dynamics_feedback.ipynb](notebooks/dynamics_feedback.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `instability_graph` | `experiment/data/input_snapshots/instability_guided_bp/02_dynamics_feedback/instability_graph` |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `trace` | `experiment/data/input_snapshots/instability_guided_bp/02_dynamics_feedback/instability_graph/trace` |
| `graph_intervention_population_cases_csv` | `experiment/data/input_snapshots/instability_guided_bp/02_dynamics_feedback/instability_graph/trace/graph_intervention_population_cases.csv` |

1. Read the explicit historical stage-01 input snapshot and evaluate raw gain versus eventual repair.
2. Ablate metrics, combine overlap/gain, compare screened and uniform candidate pools, and examine trajectory redirection.
3. Evaluate early intervention, fixed-region feedback, stopping, gating and candidate budgets.
4. Select policies on development data, then run the cross-p/distance and final stochastic comparisons.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

The preceding stage leaves two explanations open: instability reduction may signal a useful trajectory change, or the score may reward smoothness and candidate-pool bias. This follow-up separates the metric, the available actions and the decision policy before evaluating a full decoder.

Raw instability changes preserve a common measurement scale across the before/after probe. Percentile scores recomputed within each state can change their reference distribution, so a percentile decrease need not be an absolute decrease. The metric-ablation blocks test this distinction rather than treating all gain definitions as equivalent.

| Blocks | Main controlled comparison | Question |
| --- | --- | --- |
| A–C | Raw gain, alternative metrics, overlap and their joint predictors | Is response informative beyond where the candidate acts? |
| D–E | Hybrid targeting and a uniformly constructed candidate pool | Does an apparent benefit survive candidate-generation bias? |
| F–G | Stabilization versus redirection, including early stable-wrong cases | Does reducing motion actually change the logical outcome? |
| H–K | Fixed-region feedback, stopping, intervention gating and budget | Can repeated decisions avoid moving targets, wasted probes and regressions? |
| L–N | Candidate structure, failure type, p and distance | Which associations generalize beyond the initial cohort? |
| O | Frozen full stochastic decoder comparison | Is the final accuracy/cost tradeoff competitive on the whole population? |

Candidate observations from one decoder input are related; grouped validation prevents views of the same input from appearing on both sides of a prediction comparison. Failure-cohort repair rates and full-population LER answer different questions. In particular, favorable local repair statistics do not override a worse final stochastic LER, and feedback thresholds should not be justified using the final benchmark's truth labels.

## Current results and conclusions

The saved candidate audit gives only weak discrimination between useful and ineffective interventions from trajectory-gain scores. Hybrid targeting improves conditional repair on the selected cohort, but it also causes more regressions than the random control.

The saved final feedback comparison does not outperform BPOSD in the illustrated high-noise regime. Local targeting gains therefore do not establish overall decoder superiority; the candidate-level diagnostic and the complete feedback decoder must be evaluated separately.

Evidence: [raw_gain_vs_repair.csv](../../../../data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/raw_gain_vs_repair.csv), [hybrid_targeting_ablation.csv](../../../../data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/hybrid_targeting_ablation.csv), [final_stochastic_benchmark.csv](../../../../data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/final_stochastic_benchmark.csv), [unbiased_candidate_validation.csv](../../../../data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats/unbiased_candidate_validation.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/stats` |
| trace | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/trace` |
| config | `experiment/data/decoder_intervention/instability_guided_bp/02_dynamics_feedback/main/config` |
| plots | `experiment/research/instability_guided_bp/stages/02_dynamics_feedback/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

This follow-up consumes its named historical input snapshot, not automatically the latest stage-01 output. Candidate-level correlation, selected-failure repair and unbiased LER are different evaluations.

Treat the negative final comparison as a constraint on the feedback claim; a better selection signal is still needed.

Return to [the research-theme workflow](../../README.md).
