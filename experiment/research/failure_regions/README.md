# Failure-region sampling and case visualization

Locate parameter regimes where BPOSD and Tesseract disagree, retain physical samples with detector/error/observable data, and inspect representative decoding trajectories. This is both a data producer and a case-analysis workspace; it is not one sequential notebook.

## Structure and entries

```text
failure_regions/
  README.md
  notebooks/
    ca_bposd_tesseract_failure_region.ipynb
    failure_case_analysis.ipynb
    belief_visualization.ipynb
  config/paths.json
  plots/
```

| Notebook | Work performed | Execution order |
| --- | --- | --- |
| [ca_bposd_tesseract_failure_region.ipynb](notebooks/ca_bposd_tesseract_failure_region.ipynb) | Read/collect benchmark curves; select gap regimes; acquire outcome-labeled cases under code-capacity, phenomenological and circuit-level models | Producer for the structural studies |
| [failure_case_analysis.ipynb](notebooks/failure_case_analysis.ipynb) | Run a separate one-round/perfect-measurement case acquisition; count outcome groups, inspect error-weight distributions and selected belief traces | Independent analysis branch with its own main-case paths |
| [belief_visualization.ipynb](notebooks/belief_visualization.ipynb) | Trace selected beliefs and display static, Plotly and animated graph views across the four outcome groups | Needs the named failure cases and trapping-structure outputs |

The `CA_*` names in this legacy entry are source/API names for the fixed triangular codes. Use conda `qec` and kernel working directory `Decoders/`.

## Data flow and experiment procedure

1. Initialize the selected notebook and read its [explicit path manifest](config/paths.json). Benchmark references are historical input snapshots.
2. For acquisition, select the notebook's noise/round/measurement settings and benchmark-defined p regime. Benchmark and case-acquisition cells have their own sampling budgets and selection conditions.
3. Generate cases/counts with the configured BPOSD and Tesseract. Retain the actual code, p/error components, syndrome, error support, true logical observables and decoder predictions where the collector supplies them.
4. Use `trapping_sets` for structural feature extraction. The graph-visualization notebook reads those feature files explicitly rather than reconstructing a conclusion from filenames or images.
5. Analyze each acquisition branch separately, then inspect representative belief/graph examples. A historical example missing from the new sample is skipped with an explanation.

Outcome abbreviations are `bfts` = BPOSD fails/Tesseract succeeds, `bstf` = BPOSD succeeds/Tesseract fails, `bs` = both succeed, and `bf` = both fail. Saved examples can be selected, deduplicated or capped; counts of retained cases are not automatically physical LER estimates.

## Current results and conclusion

The saved main-case acquisition contains both BPOSD-failure/Tesseract-success and the reverse disagreement direction, together with shared successes and failures. These retained populations provide contrasting examples; their selection and capping prevent interpreting their relative sizes as unbiased LER or a universal decoder ranking.

The multi-noise acquisition is incomplete. The downstream structural/visualization pipeline has no completed current evaluation: 该部分实验未进行，可参考data_packedup自行分析

## Inputs, outputs and dependencies

The two acquisition branches retain different output directories beneath `experiment/data/failure_analysis/failure_regions/root/main/`. In particular, the multi-noise producer uses `trace/<noise_model>/cases/`, while the independent main-case analysis uses `trace/main_cases/trace/<case-setting>/`. They are not interchangeable aliases. Benchmark/aggregate tables stay under `stats/`; per-case records and animations are separate.

Local figures and animations stay under `plots/`. The manifest's inputs include benchmark snapshots and the graph-feature paths produced by [trapping_sets](../trapping_sets/README.md). Detailed keys, placeholders and exact current paths are documented in [config/README.md](config/README.md).

The dependency chain for visualization is acquisition → structural extraction → visualization. Replaying a few representative traces cannot establish a population-level failure mechanism. This README records observed output availability; the unfinished producer and visualization have no completed conclusions.
