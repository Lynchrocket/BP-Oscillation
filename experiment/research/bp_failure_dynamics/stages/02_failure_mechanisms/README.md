# BP failure taxonomy, stopping and repair

Separate dynamical syndrome failures, premature first-valid failures and persistent logical failures, then test which can be repaired by stopping changes or Adaptive BP.

## Directory and entry

```text
02_failure_mechanisms/
  README.md
  notebooks/bp_failure_mechanisms.ipynb
  config/
```

Run [bp_failure_mechanisms.ipynb](notebooks/bp_failure_mechanisms.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `error_weight` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight` |
| `weight_2_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_2_data.csv` |
| `weight_3_data_csv` | `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_3_data.csv` |

1. Reconstruct the weight-2/weight-3 diagnostic population and validate native trace semantics.
2. Classify complete trajectories by syndrome consistency, first-valid behavior, long-run stability and logical sector.
3. Compare stopping rules and stable endpoints; characterize the residual-operator Tanner structures.
4. Evaluate Adaptive BP by failure type, add the distance-9 weight-4 cohort, and run longer-BP, restart and random-redundant-view controls at comparable edge-update budgets.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

Stage 01 establishes why a single stable/unstable label is insufficient. This stage asks whether a failure is caused by an unresolved constraint, an unfortunate stopping time, or persistence in the wrong logical sector. Classification uses the full diagnostic trajectory; it is an offline explanation, not a runtime policy with access to truth.

| Failure family | Diagnostic distinction | Intervention being tested |
| --- | --- | --- |
| Dynamic syndrome failure | The relevant BP endpoint does not satisfy the syndrome | Additional iteration, restart or a changed graph may resolve the constraint |
| Premature first-consistent failure | The first consistent correction is wrong, but a later consistent correction is correct | Changing the stopping rule may suffice without changing the graph |
| Stable persistent logical failure | A wrong logical sector persists with small late message motion | Stabilization alone is unlikely to address the error |
| Unstable persistent logical failure | Wrong logical output persists while messages remain active | Test trajectory redirection rather than assuming eventual convergence fixes the sector |

The stopping experiment measures when syndrome consistency and logical correctness first occur and whether they persist. The residual-operator analysis studies the difference between correction and sampled error, connecting a logical failure to its Tanner subgraph. Adaptive BP is then evaluated separately by failure type rather than averaged across unrelated mechanisms.

The equal-budget comparison includes longer plain BP, restarted BP and random redundant views. Its purpose is to separate extra computation, access to alternative graph representations and the benefit of adaptive selection. Interpret repairs on the selected failure cohort together with those controls: a large repair fraction alone does not establish a targeting advantage or an improvement in unbiased LER. Cross-distance comparisons also retain their own sampling and failure-type denominators.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

Restart controls now use independent Gaussian V2C initialization noise with a reproducible seed per sample/restart. Priors, PCM, syndrome and the allocated sweep budget are unchanged. `restart_message_sigma=0` reproduces deterministic initialization. The prior restart repair counts are historical.

## Current results and conclusions

The historical equal-budget comparison shows that Adaptive BP and random redundant graph views repair more of the selected failures than simply extending plain BP. The strong random-view control prevents attributing the whole benefit to adaptive targeting. Failure-type and stopping-rule tables support separating unresolved syndrome constraints, premature stopping and persistent logical errors.

The saved restart comparison predates randomized message initialization and cannot establish the effectiveness of the corrected restart control. The revised restart experiment has no regenerated result tables: 该部分实验未进行，可参考data_packedup自行分析

Evidence: [adaptive_equal_budget_comparison.csv](../../../../data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/main/stats/adaptive_equal_budget_comparison.csv), [failure_type_population_summary.csv](../../../../data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/main/stats/failure_type_population_summary.csv), [stopping_rule_benchmark.csv](../../../../data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/main/stats/stopping_rule_benchmark.csv), [cross_distance_summary.csv](../../../../data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/main/stats/cross_distance_summary.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/stats` |
| trace | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/trace` |
| config | `experiment/data/bp_dynamics/bp_failure_dynamics/02_failure_mechanisms/principle_fixes_v1/config` |
| plots | `experiment/research/bp_failure_dynamics/stages/02_failure_mechanisms/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

The table measures repair within an already selected failure cohort. It is not an unbiased LER or proof that a particular trapping structure causes failure. Compare actual/configured edge updates alongside timing.

Stage 03 investigates logical-sector ambiguity behind stable logical failures rather than assuming another stopping rule resolves them.

Return to [the research-theme workflow](../../README.md).
