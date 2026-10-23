# Logical-sector ambiguity and susceptibility

Test whether stable logical failures are associated with competition between logical sectors, and whether targeted perturbations reveal useful sector changes.

## Directory and entry

```text
03_logical_degeneracy/
  README.md
  notebooks/instability_degeneracy_experiments.ipynb
  config/
```

Run [instability_degeneracy_experiments.ipynb](notebooks/instability_degeneracy_experiments.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |

1. Reconstruct BP outcome classes from the shared fixed-weight inputs and form matched cohorts.
2. Compute logical-sector masses by exact nullspace enumeration where the configured problem size permits it.
3. Compare top-two sector gaps, entropy and true-sector gaps, including paired tests and AUROC estimates.
4. Perturb instability-selected, low-reliability and random variables; summarize valid sector switches and discovery of the correct sector.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

A stable wrong correction motivates a different question from oscillation: does the posterior distribute substantial probability across competing logical sectors? For each syndrome, the exact calculation sums the configured independent error-prior probabilities within each sector, producing masses `Z_l` and normalized probabilities `p_l = Z_l / sum(Z)`. This is a statement about the notebook's DEM prior model.

The top-two log-mass gap is `log(Z_top1) - log(Z_top2)`. A smaller gap indicates stronger sector competition. Logical entropy `-sum(p_l log(p_l))`, normalized by the log of the number of sectors, measures the distribution across all sectors. The true-sector margin uses the known sampled logical label and is therefore an offline diagnostic. These quantities are different from the BP hard-decision confidence.

| Experiment | Control | What it can establish |
| --- | --- | --- |
| Exact sector masses | Audit coset enumeration and probability normalization | The ambiguity measurements are computed consistently for the chosen finite code |
| Matched stable endpoints | Compare stable-wrong and stable-correct cases under the matching procedure | Association between logical ambiguity and stable failure, with paired uncertainty |
| Prediction comparison | Compare gap/entropy with late residual on the same cases | Whether sector ambiguity explains distinctions that message motion misses |
| Logical susceptibility | Compare instability-targeted, random and reliability-based perturbations | Whether perturbation location changes the chance of discovering another consistent sector |

For susceptibility, keep the fraction of syndrome-consistent outputs, the fraction changing logical sector and the fraction becoming correct separate. A perturbation that breaks the syndrome constraint is not a successful sector switch. Exact enumeration and truth-based matching are diagnostic tools; neither is assumed available to an online decoder, and correlation with ambiguity does not establish that degeneracy causes every stable failure.

## Current results and conclusions

On the saved matched population, stable wrong endpoints exhibit stronger logical-sector competition than stable correct endpoints. Sector-gap and entropy diagnostics distinguish the selected outcome classes more clearly than late message residuals. This supports logical ambiguity as a useful diagnostic association beyond convergence.

The exact-mass audit supports the finite-instance calculation. Neither the matched comparison nor its probability-normalization checks establish causality, scalability to arbitrary codes, or an online policy that can inspect the true sector.

Evidence: [logical_ambiguity_auc.csv](../../../../data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_auc.csv), [logical_ambiguity_paired_comparison.csv](../../../../data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_ambiguity_paired_comparison.csv), [exact_logical_mass_summary.csv](../../../../data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/exact_logical_mass_summary.csv), [logical_susceptibility_group_comparison.csv](../../../../data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats/logical_susceptibility_group_comparison.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/stats` |
| trace | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/trace` |
| config | `experiment/data/bp_dynamics/bp_failure_dynamics/03_logical_degeneracy/main/config` |
| plots | `experiment/research/bp_failure_dynamics/stages/03_logical_degeneracy/plots/main` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

The single-feature AUROC rows are not cross-validated; the table explicitly labels the combined cross-validated models. Perturbation comparisons are condition-specific and do not establish a universally superior targeting rule.

Use ambiguity as a candidate mechanism to test in the gap/PCM studies; do not treat it as an already proven universal cause.

Return to [the research-theme workflow](../../README.md).
