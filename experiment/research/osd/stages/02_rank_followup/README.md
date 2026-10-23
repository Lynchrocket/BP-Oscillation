# Repair-variable ranks and snapshot selection

Explain which repair variables are misplaced in the reliability ordering and determine whether practical rank/snapshot rules can capture the oracle opportunity.

## Directory and entry

```text
02_rank_followup/
  README.md
  notebooks/rank_followup.ipynb
  config/
```

Run [rank_followup.ipynb](notebooks/rank_followup.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

The manifest selects the following inputs explicitly. A historical snapshot remains a historical snapshot; it is not automatically replaced by a newer producer result.

| Input key | Repository-relative source |
| --- | --- |
| `reliability_dynamics` | `experiment/data/input_snapshots/osd_reliability/02_rank_followup/reliability_dynamics` |
| `event_screen_csv` | `experiment/data/input_snapshots/osd_reliability/02_rank_followup/reliability_dynamics/{_input_run}/event_screen.csv` |
| `unbiased_cases_csv` | `experiment/data/input_snapshots/osd_reliability/02_rank_followup/reliability_dynamics/{_input_run}/unbiased_cases.csv` |

1. Read the named reliability-dynamics input snapshot and reconstruct trace-native gap cohorts.
2. Measure repair Top-K coverage, mis-ranking features and oracle rank/snapshot ceilings.
3. Fit rank/snapshot rules on development data, evaluate pairwise models and BP parameter sensitivity.
4. Run the separate gap-regime benchmark and selective OSD frontier, tracking repairs and regressions.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

This stage turns an instability/reliability association into a more precise repair-ranking question. The diagnostic repair set is the support of `e_BP XOR e_OSD_success`, the difference between a BP correction and a successful OSD correction. It is not necessarily the sampled physical error support, and depends on which successful correction is used.

Ranks increase with reliability, so rank 1 is least reliable. Repair Top-K coverage measures how much of that repair set lies among the K least-reliable variables. These labels explain ranking failures offline; a runtime ranking model cannot inspect the successful correction or true logical label of the evaluated shot.

| Comparison | What is varied | What is held conceptually separate |
| --- | --- | --- |
| Mis-ranking predictors | Magnitude, instability and trajectory features | Association versus useful held-out prediction |
| Oracle rank ceiling | Give the ordering unavailable repair information | Ordering opportunity versus a practical rank rule |
| Oracle snapshot ceiling | Choose a successful time using truth | Snapshot opportunity versus causal snapshot selection |
| Practical models | Runtime snapshot rules, rank corrections and pairwise ranking | Development fitting versus held-out evaluation |
| Final benchmark/frontier | Apply frozen rules and selective OSD | Conditional repair versus regressions, LER and invocation cost |

BP parameter scans test whether a different soft-information generator improves these rankings before changing OSD itself. Preserve shot-level grouping when several snapshots originate from one input. The reported practical negative result remains central: an oracle opportunity is not evidence that the current dynamic-rank model can exploit it. Sparse distance-9 events also cannot establish cross-distance generalization.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

## Current results and conclusions

Saved snapshot-oracle results show that useful soft information can occur before the final BP endpoint. The practical dynamic-rank rule does not realize that opportunity in the illustrated benchmark and can worsen decoding relative to the baseline.

Oracle access to truth is therefore an opportunity bound, not evidence that the current selector works. The larger-distance gap population is too sparse to support generalization. These conclusions describe the historical OSD implementation.

Evidence: [experiment_status.csv](../../../../data/osd/osd_reliability/02_rank_followup/main/stats/experiment_status.csv), [snapshot_oracle_ceiling.csv](../../../../data/osd/osd_reliability/02_rank_followup/main/stats/snapshot_oracle_ceiling.csv), [oracle_rank_ceiling.csv](../../../../data/osd/osd_reliability/02_rank_followup/main/stats/oracle_rank_ceiling.csv), [final_gap_regime_benchmark.csv](../../../../data/osd/osd_reliability/02_rank_followup/main/stats/final_gap_regime_benchmark.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/osd/osd_reliability/02_rank_followup/principle_fixes_v1/stats` |
| trace | `experiment/data/osd/osd_reliability/02_rank_followup/principle_fixes_v1/trace` |
| config | `experiment/data/osd/osd_reliability/02_rank_followup/principle_fixes_v1/config` |
| plots | `experiment/research/osd/stages/02_rank_followup/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

Gap closure may be negative or exceed one depending on the paired LER denominator; it is not a probability. The small distance-9 event count prevents a robust distance-9 conclusion.

This stage motivates better causal selection, not deployment of the current dynamic-rank rule.

Return to [the research-theme workflow](../../README.md).

## JSON execution switches

Execution switches are read from the stage's JSON configuration; environment variables do not override them.

| Key | Default | Values and behavior |
| --- | --- | --- |
| `smoke_mode` | `false` | `false`: use the configured experiment budget. `true`: apply the notebook's existing reduced sampling, iteration and analysis settings and redirect generated outputs to a temporary directory. |
| `reuse_existing_data` | `true` | `true`: reuse available compatible cached results. `false`: recompute the corresponding results. Existing backend and algorithm version checks still apply. |

Set these booleans before running the notebook's configuration cell. Smoke mode does not implicitly disable cache reuse. Temporary smoke outputs are isolated from the configured research outputs.

`config/reliability_rank_followup.json` controls smoke mode for both sections and cache reuse for the main section. `config/followup_config.json` independently controls cache reuse for the follow-up section through its own `reuse_existing_data` field. It has no separate `smoke_mode` field.
