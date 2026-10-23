# Locate the BPOSD–Tesseract gap

Establish paired decoder outcomes before optimizing stopping, OSD order or handoff. Determine whether the sampled regime actually contains the target BPOSD-fail/Tesseract-success events.

## Directory and entry

```text
00_gap_exploration/
  README.md
  notebooks/gap_exploration.ipynb
  config/
```

Run [gap_exploration.ipynb](notebooks/gap_exploration.ipynb) with conda `qec` and working directory `Decoders/`. Run its initialization and definition cells before the scientific sections. Every stage owns its notebook namespace; a preceding notebook need not have been run in the same kernel.

## Inputs and execution flow

This stage generates its own sampled populations; its path manifest has no external dataset input. Related earlier stages provide scientific motivation, not hidden in-memory variables.

1. Generate paired physical shots for the configured distances/p values and classify shared successes, shared failures and disagreements.
2. Trace BP and audit OSD order, endpoint/snapshot reliability and oracle opportunities.
3. Evaluate stopping rules, false-stop features, selective audits and BP soft-information variants.
4. Freeze proposed policies and compare paired LER, OSD invocation and cost on the held-out population.

The notebook headings follow this workflow; the configuration reference is [config/README.md](config/README.md). It documents the active JSON files, actual defaults and initialization constraints.

## Scientific rationale and comparison design

This stage tests the proposed explanation that standard BPOSD sometimes accepts a wrong syndrome-consistent BP result before OSD can use its soft information. That explanation predicts a specific paired event: standard decoding fails, while a forced OSD call and the reference decoder succeed on the same shot. It is a hypothesis to screen for, not an assumed property of every error-rate regime.

| Experiment block | Controlled comparison | Question answered |
| --- | --- | --- |
| Paired gap and OSD order | Identical physical shots, PCM, syndrome, priors and logical truth | Is there a measurable gap, and which orders repair the observed disagreement? |
| BP snapshots and ordering | Change the handoff snapshot or reliability input while retaining the decoding problem | Is useful OSD information already present earlier in the trajectory? |
| BP-only stopping | Compare implementable stopping features with offline trajectory/oracle labels | Can a causal stopping rule avoid a false stop? |
| Selective OSD and OSD-0 audit | Track invocation rate, repairs and regressions, including the zero-order reference | Is the gain due to selection, the OSD search or inconsistent baseline semantics? |
| Soft-generator and final benchmark | Freeze the tested rule before paired accuracy/cost evaluation | Does improved diagnostic behavior translate into a practical decoder? |

An oracle that chooses a successful snapshot using truth is an upper bound on an available opportunity, not a runtime method. Likewise, a selected disagreement cohort cannot estimate physical-noise LER. If a sampled regime contains no target gap events, its screening result is informative but downstream conditional estimates cannot substantiate a mechanism for that missing cohort.

When reporting gap closure, use `(LER_base - LER_method) / (LER_base - LER_reference)` only with its paired estimates and a nonzero denominator. It can be negative or greater than one and should not be read as a success probability.

## Implementation revision and historical results

The corrected implementation writes to the `principle_fixes_v1` run selected by the local manifest. Existing `main` results and input snapshots are preserved. Any saved observations discussed below describe the earlier implementation and do not validate this revision. No result tables have been regenerated for the corrected implementation.

Corrected-implementation evaluation: 该部分实验未进行，可参考data_packedup自行分析

OSD now follows the benchmark `osd_cs` definition: signed BP column ordering, original channel-prior cost `sum(log(1/p_i))` over selected faults, all single non-pivot flips and pairs among the first configured-order non-pivot positions. Order zero tests only the initial solution. Magnitude/learned/oracle rank scores remain explicitly identified experimental inputs; they do not replace channel priors in candidate scoring. This matches the benchmark search, not exact posterior logical-sector decoding.

## Current results and conclusions

The saved low-noise paired population is dominated by shared decoder successes and contains no target BPOSD-failure/Tesseract-success class. The classifier status reports insufficient events. This regime therefore cannot establish that a learned gap classifier or selective repair rule closes the target gap.

A missing target class does not prove that the mechanism is absent at other noise settings, and no observed failures does not imply a zero underlying error probability. These observations describe the historical OSD implementation.

Evidence: [gap_class_summary.csv](../../../../data/decoder_gap/gap_mechanisms/00_gap_exploration/main/stats/gap_class_summary.csv), [gap_decoder_ler.csv](../../../../data/decoder_gap/gap_mechanisms/00_gap_exploration/main/stats/gap_decoder_ler.csv), [experiment_status.csv](../../../../data/decoder_gap/gap_mechanisms/00_gap_exploration/main/stats/experiment_status.csv), [final_benchmark.csv](../../../../data/decoder_gap/gap_mechanisms/00_gap_exploration/main/stats/final_benchmark.csv).

## Outputs and how to interpret them

| Output | Current repository-relative location |
| --- | --- |
| stats | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/stats` |
| trace | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/trace` |
| config | `experiment/data/decoder_gap/gap_mechanisms/00_gap_exploration/principle_fixes_v1/config` |
| plots | `experiment/research/gap_mechanisms/stages/00_gap_exploration/plots/principle_fixes_v1` |

`stats/` contains summary/aggregate tables; `trace/` contains per-case/iteration records and caches. A reserved path is not evidence that a file has been generated. Keep distance, p, cohort selection and denominator with every metric.

## Limits and next step

Oracle opportunities and held-out policy performance must be separated. Gap-closure quantities with an absent/zero denominator are not interpretable evidence of improvement. Inspect the status file before interpreting fitted policies.

Stage 01 uses a higher-p/enriched gap cohort to study mechanisms; its selected-cohort fractions must not be interpreted as physical event frequencies.

Return to [the research-theme workflow](../../README.md).
