# PCM diversity: candidate generation versus selection

Test whether equivalent PCM changes expose alternative BP fixed points, how this relates to logical near-degeneracy, and whether non-oracle selection can turn generated candidates into lower LER.

## Directory and entry

```text
02_pcm_fixed_points/
  README.md
  notebooks/pcm_fixed_points.ipynb
  config/
```

Run [pcm_fixed_points.ipynb](notebooks/pcm_fixed_points.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

This stage generates its own sampled populations; its path manifest has no external dataset input. Related earlier stages provide scientific motivation, not hidden in-memory variables.

1. Build the unbiased population and a separately enriched gap/control cohort, then compute exact sector margins.
2. Part A compares basis replacement with redundant-row addition; Part B relates margins to basins, PCM sensitivity and oscillation.
3. Part C varies distinct targeted PCM ensembles and view budgets; Part D compares sector diversity, oracle generation and non-oracle selectors.
4. Part E freezes PCM repair and OSD1 handoff policies and evaluates the final unbiased LER.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

This stage tests a proposed chain: logical-sector competition may permit multiple BP basins, and changing the PCM may change which basin is reached. Each link needs its own measurement. A correlation among ambiguity, restart sensitivity and PCM sensitivity does not alone prove that chain.

Two transformations must remain distinct. An invertible row transformation replaces the basis, `H' = T H`, with syndrome `s' = T s`. Redundant addition appends combinations, `H' = [H; R H]`, with `s' = [s; R s]`. Operations are over GF(2). Both preserve the constraint set when performed consistently; their graph size and update cost can differ. Logical observables and error priors retain their original meaning.

| Part | Purpose and comparison | Main interpretation |
| --- | --- | --- |
| A: mechanism cleanup | Separate basis replacement from redundant addition; pair successful/failing views of the same shot | Identify which representation change affects the endpoint without pooling unlike graph operations |
| B: ambiguity and basins | Relate sector margin to restart outcomes, PCM sensitivity and oscillation | Test each association in the proposed mechanism separately |
| C: distinct ensembles | Compare random, instability, ambiguity and conflict targeting as view budget increases | Measure useful diversity rather than repeatedly counting the same generated view |
| D: generation versus selection | Compare existence of a correct candidate with non-oracle selection from that candidate pool | Separate failure to produce a solution from failure to select it |
| E: frozen decoder | Evaluate PCM repair and OSD1 handoff on the unbiased final population | Determine whether diagnostic opportunity survives practical selection and regression costs |

A generation oracle succeeds if a correct candidate exists in the pool. A deployable selector must choose without truth; for the same pool its success cannot exceed that oracle. Growing oracle coverage is therefore insufficient evidence of a better decoder. Paired view analyses share a physical shot and must not be interpreted as independent new shots. Keep distinct-view count, graph/update budget and the baseline-success population visible when interpreting the final comparison.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

A fixed point now requires a stable hard decision and small tail-window changes in LLR, C2V and V2C messages. Period-2 checks include both message directions. Within each shot and PCM role, eligible endpoints are clustered by their full terminal message/LLR state, with equal hard decisions and complete-link tolerance checks. `fixed_point_count` and `basin_entropy` describe the observed endpoint clusters; `restart_sector_count` and `sector_entropy` separately describe logical sectors. Neither finite restart counts nor the finite tail window prove that all mathematical attractors have been found.

## Current results and conclusions

The saved partial studies show that alternative PCM constructions can expose additional candidate corrections. Increasing the view budget improves oracle generation in the illustrated ambiguity-targeted comparison, while the selection tables retain a gap between generating a correct candidate and choosing it without truth.

These historical observations support separating representation, candidate generation and practical selection. They do not validate the revised full-message-state definition of fixed points or basins. The complete final-policy evaluation is unavailable: 该部分实验未进行，可参考data_packedup自行分析

Evidence: [A1_basis_vs_redundant.csv](../../../../data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/main/stats/A1_basis_vs_redundant.csv), [C1_C2_distinct_targeted_scaling.csv](../../../../data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/main/stats/C1_C2_distinct_targeted_scaling.csv), [D5_generation_selection_gap.csv](../../../../data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/main/stats/D5_generation_selection_gap.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/stats` |
| trace | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/trace` |
| config | `experiment/data/decoder_gap/gap_mechanisms/02_pcm_fixed_points/principle_fixes_v1/config` |
| plots | `experiment/research/gap_mechanisms/stages/02_pcm_fixed_points/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

These are preliminary selected-cohort results while the final benchmark is unfinished. They do not establish a final LER advantage or a scalable general-purpose PCM selector.

Wait for frozen-selection and unbiased Part-E outputs before making a decoder-performance claim.

Return to [the research-theme workflow](../../README.md).

## JSON execution switches

Execution switches are read from the stage's JSON configuration; environment variables do not override them.

| Key | Default | Values and behavior |
| --- | --- | --- |
| `smoke_mode` | `false` | `false`: use the configured experiment budget. `true`: apply the notebook's existing reduced sampling, iteration and analysis settings and redirect generated outputs to a temporary directory. |
| `reuse_existing_data` | `true` | `true`: reuse available compatible cached results. `false`: recompute the corresponding results. Existing backend and algorithm version checks still apply. |

Set these booleans before running the notebook's configuration cell. Smoke mode does not implicitly disable cache reuse. Temporary smoke outputs are isolated from the configured research outputs.

Configuration: [next_gap_mechanism_pcm_fixed_point_experiments.json](config/next_gap_mechanism_pcm_fixed_point_experiments.json).
