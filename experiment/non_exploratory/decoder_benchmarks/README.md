# Decoder logical-error-rate benchmarks

[decoder_benchmark.ipynb](notebooks/decoder_benchmark.ipynb) constructs noisy triangular-code circuits, collects or resumes Sinter decoder statistics, and plots logical error rate (LER) against physical error probability. It compares decoder outcomes and sampling costs under explicitly separated circuit/noise/schedule settings. Later sections reuse saved results to inspect the BPOSD–Tesseract gap and produce figures.

## Experiments in the notebook

| Study | Comparison |
| --- | --- |
| Noise model | One-/three-round circuit-level circuits, one-/two-round phenomenological circuits, and code-capacity noise |
| Circuit noise components | Perfect measurement; additionally disabled idle noise; additionally disabled preparation/reset noise |
| CNOT schedule | Default builder schedule versus `CircuitBuilder(..., schedule=1)` |
| BP schedule | Default parallel minimum-sum updates versus `schedule='serial'`, including selected CNOT/noise variants |
| Minimum-sum scaling | `alpha=0.8` and `1.2` for one-round perfect-measurement circuits; `alpha=0.8` for three rounds |
| OSD invocation | Standard BPOSD versus `force_osd=True`, with Tesseract as a reference; includes a distance-9 case |
| Adaptive decoding | AdaptiveBPOSD, AdaptiveBP, MBBP-LD and ExpandedAdaptiveBP against BPOSD/Tesseract where the corresponding section constructs them |
| Quotient splitting | `QuotientSplit(max_iter=50, internal_iter=5, damping=0.2)` against BPOSD/Tesseract |

The initial code set is the fixed 5-bit/distance-3 and 9-bit/distance-5 triangular codes. Dedicated later sections add `TriangularCode(3)` (distance 9). Historical `CA_*` labels remain in task metadata and saved CSVs so that existing curves and resume records remain compatible. The actual `ErrorRate(...)`, builder and decoder constructor calls define each case; headings and historical filenames are descriptive labels, not configuration parsers.

The usual p grid is `np.geomspace(1e-8, 0.1, 12)`. Most earlier sections use `max_shots=1_000_000_000` and `max_errors=2000`; later Adaptive BP/MBBP-LD/expanded and distance-9 sections use their explicitly assigned budgets, including `50_000_000` shots. These are per-task stopping limits, not a guarantee of that many completed samples.

BPOSD normally uses `max_iter=60`, `bp_method='minimum_sum'`, `osd_order=10` and `osd_method='osd_cs'`. Tesseract normally uses `beam=8`, `pq_limit=200000` and the existing `det_penalty=0.0` argument. Adaptive BP experiments use their own constructor settings; the dedicated distance-9 case sets `adaptive_max_views=8`. Every `ErrorRate` component omitted by a section inherits p; explicitly zero components remain zero. Different case settings must not be pooled as one controlled comparison.

## How to run

Use conda `qec` and set the kernel working directory to the `Decoders/` repository root.

1. Run the main initialization, imports, code/style definitions and plotting helpers.
2. Select a benchmark subsection and run its settings cell before its task-construction/collection and plotting cells. Some subsections change the code set or decoder settings; rerun the relevant setup when switching between them.
3. Collection calls `sinter.collect` with `multiprocessing.cpu_count()` workers and `save_resume_filepath=save_file`. Reusing the same filename resumes that experiment; use a separate output path when deliberately changing scientific settings.
4. To inspect existing results without sampling, run the initialization/settings and the cells that read the corresponding CSV, then their plotting cells. The `Initial Benchmark` section reads the main baseline CSV instead of collecting it again. Its input must already exist. The final `Benchmark Figures` section reads the saved three-round circuit-level CSV.

The three sections have independent initialization cells but share explicit local path/runtime manifests. Their scope-specific variables keep the baseline analysis and final figure setup separate from the main benchmark loop.

## Inputs, outputs and interpretation

[config/paths.json](config/paths.json) names the saved benchmark CSVs and collection filename templates. [config/README.md](config/README.md) documents the schema, every current path key, template arguments and logging control. No configuration reference outside this experiment is required.

Results live under `experiment/data/benchmark/decoder_benchmark/root/main/stats/benchmark/`; figures stay under this experiment's `plots/` using the manifest's explicit destinations. The shared producer name remains singular for compatibility with existing datasets. `trace/` and effective `config/` output roots are reserved in the manifest; their presence does not imply that each section writes a trace or configuration snapshot.

A Sinter record contains task metadata, shots, logical errors, discards, elapsed seconds and related counters. Task metadata identifies code, decoder and p; main collection tasks also record the noise model. LER is errors divided by shots in the notebook's summary helper. That helper additionally exposes seconds and milliseconds per shot. Curve grouping follows `(code, decoder)` within the selected experiment; keep its noise model, round count and other conditions with the result.

Plots use logarithmic axes. The main curve plots set `highlight_max_likelihood_factor=10`; these likelihood regions should not be relabeled as 95% confidence intervals. The initial analysis pivots BPOSD and Tesseract LERs by code/p and computes `gap_bposd_minus_tess = BPOSD - TESSERACT`; a positive gap means the saved BPOSD LER is higher. It does not establish a failure mechanism. Timing columns are measurements of that execution and should not be treated as a controlled backend-speed comparison across unrelated runs.

## Qualitative observations from saved data

The saved paired benchmark tables contain both positive and negative BPOSD-minus-Tesseract error-rate differences across code and noise settings. No single ranking follows across every condition. Some low-noise settings contain no observed failures; equal empirical rates there do not establish equivalent decoders or a zero underlying failure probability.

Noise components, round count, BP/CNOT schedule and decoder settings define separate comparisons. Differences between files cannot be assigned to one mechanism unless those conditions are controlled. The benchmark provides performance observations; the research studies test possible explanations. These saved benchmarks were retained from the backup and were not regenerated with the recent implementation changes.
