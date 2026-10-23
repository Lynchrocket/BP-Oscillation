# Attribute mechanisms within the decoder gap

Distinguish premature stopping, PCM sensitivity, restart basins, update-rule sensitivity and logical ambiguity using paired interventions on actual gap cases.

## Directory and entry

```text
01_mechanism_attribution/
  README.md
  notebooks/mechanism_attribution.ipynb
  config/
```

Run [mechanism_attribution.ipynb](notebooks/mechanism_attribution.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

This stage generates its own sampled populations; its path manifest has no external dataset input. Related earlier stages provide scientific motivation, not hidden in-memory variables.

1. Generate an unbiased physical-shot population, then construct the separately labeled enriched mechanism cohort.
2. Audit strict early stop and OSD-order decomposition; compare equivalent PCM views and random restarts.
3. Probe update rules, damping/scaling, exact logical-sector margins and prior perturbations.
4. Combine overlapping mechanism labels, compare targeted rows/subgraphs, and evaluate frozen repair policies in a separate paired benchmark.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

The primary diagnostic cohort consists of shots where standard BPOSD fails and Tesseract succeeds. The goal is to learn which changes repair those same shots. A forced-OSD repair, a graph-view repair and a restart repair can overlap; they are intervention sensitivities rather than mutually exclusive causal labels.

| Mechanism probe | Intervention on the same input | Interpretation boundary |
| --- | --- | --- |
| Early acceptance and search order | Continue BP or apply the tested OSD orders | Distinguish stopping opportunity from search-depth sensitivity |
| PCM dependence | Decode equivalent parity-check representations with transformed syndromes | Test sensitivity to graph representation while preserving the constraint problem |
| Multiple basins | Repeat initialization/scheduling experiments | Test whether alternative trajectories reach different consistent sectors |
| Update dynamics | Change the tested update rule, damping or scaling | Determine whether the endpoint depends on the numerical dynamics |
| Logical competition | Compute sector masses and prior sensitivity for the tractable code | Relate repairs to ambiguity rather than assuming a unique high-confidence sector |
| Targeted repair | Compare ambiguity-informed rows and other repair variants | Test whether a diagnostic association supports useful intervention |

Retain per-shot intervention outcomes before aggregating mechanism frequency or overlap. Candidate selection may use the permitted decoder features; sampled logical truth is reserved for evaluation and oracle comparisons. The final cross-p/cross-distance benchmark has a different role from the enriched diagnostic cohort: it measures repairs and newly introduced failures on a population that includes baseline successes.

Equivalent PCM experiments in this stage motivate the next stage's stricter separation between replacing a basis and adding redundant checks. A result attributed to generic PCM sensitivity here does not by itself identify which of those operations is responsible.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

## Current results and conclusions

The historical enriched cohort exhibits overlapping PCM-sensitivity and alternative-basin labels, while premature stopping explains only part of the selected gaps. This motivates examining graph representation and trajectory alternatives rather than attributing all disagreements to the first-consistent stopping rule.

The separately sampled population tables show that disagreement depends on the noise regime. Enriched-cohort mechanism frequencies cannot be interpreted as physical-noise LER or as an exclusive partition of causes. These observations require regeneration before they can be attributed to the corrected OSD implementation.

Evidence: [mechanism_frequency.csv](../../../../data/decoder_gap/gap_mechanisms/01_mechanism_attribution/main/stats/mechanism_frequency.csv), [population_summary.csv](../../../../data/decoder_gap/gap_mechanisms/01_mechanism_attribution/main/stats/population_summary.csv), [early_stop_summary.csv](../../../../data/decoder_gap/gap_mechanisms/01_mechanism_attribution/main/stats/early_stop_summary.csv), [final_benchmark.csv](../../../../data/decoder_gap/gap_mechanisms/01_mechanism_attribution/main/stats/final_benchmark.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/stats` |
| trace | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/trace` |
| config | `experiment/data/decoder_gap/gap_mechanisms/01_mechanism_attribution/principle_fixes_v1/config` |
| plots | `experiment/research/gap_mechanisms/stages/01_mechanism_attribution/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

These interventions show sensitivity, not a unique causal explanation for every case. Compare current native results under their recorded backend/order; do not assume bitwise identity with historical Python trajectories.

Stage 02 separates generating a useful candidate through PCM diversity from selecting that candidate without access to ground truth.

Return to [the research-theme workflow](../../README.md).

## JSON execution switches

Execution switches are read from the stage's JSON configuration; environment variables do not override them.

| Key | Default | Values and behavior |
| --- | --- | --- |
| `smoke_mode` | `false` | `false`: use the configured experiment budget. `true`: apply the notebook's existing reduced sampling, iteration and analysis settings and redirect generated outputs to a temporary directory. |
| `reuse_existing_data` | `true` | `true`: reuse available compatible cached results. `false`: recompute the corresponding results. Existing backend and algorithm version checks still apply. |

Set these booleans before running the notebook's configuration cell. Smoke mode does not implicitly disable cache reuse. Temporary smoke outputs are isolated from the configured research outputs.

Configuration: [bposd_tesseract_gap_mechanism_attribution_experiments.json](config/bposd_tesseract_gap_mechanism_attribution_experiments.json).
