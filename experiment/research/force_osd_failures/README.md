# Disagreements when OSD is forced

Determine whether BPOSD–Tesseract disagreements remain when the BPOSD path is forced to execute OSD, and inspect the retained physical-error populations rather than attributing all differences to BP early stopping.

## Structure and entry

```text
force_osd_failures/
  README.md
  notebooks/force_osd_failure_analysis.ipynb
  config/paths.json
```

Run [force_osd_failure_analysis.ipynb](notebooks/force_osd_failure_analysis.ipynb) from `Decoders/` with conda `qec`. [config/README.md](config/README.md) documents the actual inputs/outputs; scientific settings remain in code cells.

## Procedure

1. Read the named historical one-round/perfect-measurement benchmark snapshot and select the notebook's sampling p/code conditions.
2. Call `collect_samples_with_errors(..., build_decoders='force_osd')` separately for the four outcome groups. The helper constructs the forced-OSD BPOSD and Tesseract pair.
3. Save non-empty batches with one CSV header and inspect detector/error weights and decoder predictions. Keep noise settings and selection conditions with each row.
4. Plot retained error-weight distributions and inspect the directional disagreement groups. This is a comparison of sampled cases, not an equal-budget LER benchmark.

`bfts` means BPOSD fail/Tesseract success; `bstf` the reverse; `bs` both success; `bf` both fail. These outcomes are evaluated against true logical observables by the collector, not by assuming syndrome consistency implies logical correctness.

## Current results and conclusion

The saved case files retain disagreements in both directions even when OSD is forced. BP early stopping alone therefore does not explain all observed BPOSD–Tesseract differences in this acquisition.

These are selected and capped physical-error cases. Their relative sizes do not establish an overall decoder ranking, an unbiased LER, or a specific cause of OSD failure. Supporting case files are under `experiment/data/failure_analysis/force_osd_failures/root/main/trace/force_OSD/`.

## Outputs and limits

The current result files are per-case **trace** CSVs. This run has no aggregate CSV in the reserved `stats/` directory, and its completion does not imply that every declared output path has been populated. Figures are local under the manifest's `plots` path.

The collector uses stochastic sampling without a fixed sampler seed in this path, so a new run need not reproduce the same rows/counts. A mechanism claim needs controlled interventions or the paired studies in [gap_mechanisms](../gap_mechanisms/README.md), rather than inference from these distributions alone.
