# Fixed-weight decoding failures and scaling sensitivity

Acquire reusable DEM error-support cases, inspect why BPOSD and Tesseract disagree at fixed weight, and compare BP/BPOSD/Tesseract across code distance, p and minimum-sum scaling. Structural trapping-set comparisons are consumers in [trapping_sets](../trapping_sets/README.md), not additional acquisition here.

## Structure and entry

```text
low_weight_failures/
  README.md
  notebooks/low_weight_failures.ipynb
  config/
    paths.json
    runtime.json
    prior_scaling.json
    prior_scaling_paths.json
  plots/
    prior_scaling/
```

Use [low_weight_failures.ipynb](notebooks/low_weight_failures.ipynb), conda `qec`, and working directory `Decoders/`. The [configuration reference](config/README.md) describes the exact knobs, allowed prior-scaling values and paths. Run each part's initialization before its own calculation; the prior-scaling part constructs fresh code/matrix/decoder contexts rather than inheriting a preceding part's distance.

## Experiment flow

1. **Fixed-weight acquisition:** read the named historical failure-case snapshot to select conditions; enumerate weight-2/3/4 DEM error supports and save errors, syndromes, corrections, observables and diagnostic features.
2. **Case diagnostics:** compare LLR weights and diameters and inspect representative BP beliefs, error explanations and cross-layer hyperedges. A pictured case is illustrative, not a population conclusion.
3. **Weight probability/discrepancy:** sample fixed supports at increasing weights and summarize paired BPOSD/Tesseract outcomes. Keep its fixed seed and support-sampling distribution with the result.
4. **Independent prior/scaling comparison:** use exhaustive weight-2 supports for both configured distances and a fixed-seed weight-4 support sample for the larger distance. Reuse the same support IDs across p/alpha within each cohort; sweep the two named p grids and scaling settings, then compute discrete failure/recovery boundaries.
5. **Threshold-suppression examples:** run the final explicitly configured examples separately from the above scans. Do not present a finite-case boundary as an asymptotic code threshold.

BP-only uses the existing Cython trace with parallel minimum-sum, sequential accumulation and first-syndrome-consistent stopping. BPOSD is a separate existing C++ decoder measurement including OSD fallback; Tesseract is evaluated independently of alpha. For every new per-support measurement, `success` requires both `syndrome_consistent` (`H @ correction == syndrome` over GF(2)) and `logical_match`. The batch decoder's final convergence flag is not a substitute for per-case success.

## Current results and conclusions

The saved parameter summaries show that BP-only, BPOSD and Tesseract have different residual failure behavior on the same fixed-support cohorts. OSD fallback substantially reduces BP-only failures in the illustrated comparison, while the decoders can still disagree on particular supports. This does not establish a universal decoder ranking or a physical-noise LER.

Failure and recovery can vary non-monotonically with scaling. The boundary tables therefore record transitions and censoring rather than assuming one guaranteed threshold. Comparing cohorts of different distances and error weights does not isolate a pure distance effect.

Evidence: [parameter_summary.csv](../../data/failure_analysis/low_weight_failures/prior_scaling/main/stats/parameter_summary.csv), [critical_alpha_distribution.csv](../../data/failure_analysis/low_weight_failures/prior_scaling/main/stats/critical_alpha_distribution.csv) and [fixed_weight_discrepancy.csv](../../data/failure_analysis/low_weight_failures/root/main/stats/fixed_weight_discrepancy.csv). Alpha behavior can be non-monotone: the boundary output records first failure, recovery transitions and censoring without interpolating a guaranteed threshold.

## Datasets, consumers and limits

- Shared acquisition: `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/weight_{2,3,4}_data.csv`. BP dynamics, quotient and trapping-set studies explicitly read these inputs.
- Prior/scaling: `experiment/data/failure_analysis/low_weight_failures/prior_scaling/main/{stats,trace,config}/`. Aggregate parameter/boundary tables are separate from per-support corrections, stable support IDs and checkpoints.
- Prior/scaling figures: `plots/prior_scaling/`; other diagnostic figures retain their explicitly configured local paths.

Historical prior-scaling filenames that omit distance are not aliases for the new d=5/d=9 cohorts; the new analysis does not consume them. Support-conditional fractions, selected physical failure examples and physical-noise LER are different populations. Cache reuse follows `reuse_existing_data` and the available parameter-point files; configuration, source and binary hashes are not checked.
