# Instability-guided graph interventions

Test whether local message-instability information helps choose redundant checks, damping or scaling that repair BP failures, rather than merely reducing a diagnostic residual.

## Directory and entry

```text
01_graph_intervention/
  README.md
  notebooks/graph_intervention.ipynb
  config/
```

Run [graph_intervention.ipynb](notebooks/graph_intervention.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |

1. Trace shared fixed-weight cases and classify failures; localize transient instability.
2. Compare unstable-variable regions with residual-operator regions.
3. Generate redundant-row candidates and probe each under equal budgets; compare gain targeting, overlap targeting, random and anti-targeted controls.
4. Evaluate local/global damping and scaling, combined interventions, iterative adaptation and cross-distance behavior.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

Even a weak global failure predictor may identify where a useful intervention should act. This stage therefore tests instability as a spatial response signal: does a candidate redundant check that overlaps an unstable region, or reduces its short-probe instability, improve decoding on that same input?

| Experiment block | Comparison | Purpose |
| --- | --- | --- |
| Transient maps and overlap | Compare unstable variables/edges with the residual error structure | Establish spatial association before proposing an action |
| Candidate short probes | Apply candidate rows, measure the probe response, then evaluate the full outcome | Determine whether local response predicts repair rather than merely smoother messages |
| Targeting ablation | Compare overlap/gain targeting with random and weight-matched random selection | Separate information in the target from candidate size and probing effort |
| Local dynamics controls | Compare local/global damping and min-sum scaling with graph changes | Test whether numerical stabilization accounts for the benefit |
| Iterative adaptation | Recompute intervention decisions and inspect failure-type effects | Determine which trajectories remain actionable after a first change |

Comparisons retain the same syndrome, priors and original logical decoding problem. Added redundant checks require the corresponding syndrome entries; a changed graph must not silently become a changed constraint problem. Ground-truth error support can explain overlap offline but is not an input to a practical targeting rule.

Judge a targeting strategy using both repairs of baseline failures and regressions on baseline successes. Matching the candidate/probe budget and row-weight distribution is necessary to distinguish useful targeting from simply trying more or larger changes. The next stage tests whether a short-probe gain remains meaningful under stricter metric and candidate-pool controls.

## Current results and conclusions

The saved ablation shows that targeted interventions repair more selected failures than the weight-matched random control, while also introducing regressions on baseline successes. Overlap targeting repairs slightly more failures than instability-reduction targeting despite the latter achieving larger probe-instability reductions.

Reducing the diagnostic is therefore not equivalent to improving decoding. The result motivates controlled targeting studies but does not by itself establish a whole-population LER advantage.

Evidence: [targeting_ablation.csv](../../../../data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/targeting_ablation.csv), [unified_intervention_comparison.csv](../../../../data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/unified_intervention_comparison.csv), [instability_residual_overlap.csv](../../../../data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/instability_residual_overlap.csv), [transient_instability_summary.csv](../../../../data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats/transient_instability_summary.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/stats` |
| trace | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/trace` |
| config | `experiment/data/decoder_intervention/instability_guided_bp/01_graph_intervention/main/config` |
| plots | `experiment/research/instability_guided_bp/stages/01_graph_intervention/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

Repeated candidate evaluations are not independent physical shots. Repair/regression use different denominators, and some overlap-based controls use information unavailable to an online decoder.

Stage 02 tests predictive value, candidate-pool bias, useful trajectory redirection and held-out feedback policies.

Return to [the research-theme workflow](../../README.md).
