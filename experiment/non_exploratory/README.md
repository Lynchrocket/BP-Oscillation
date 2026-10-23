# Decoder benchmarks and implementation checks

This directory contains four independent notebooks. Benchmarks measure decoder logical error rates over sampled circuits; checks inspect specific decoding cases, circuit graphs or agreement with saved noise-model reference curves.

| Notebook | Purpose | Main inputs and outputs |
| --- | --- | --- |
| [Decoder benchmark](decoder_benchmarks/notebooks/decoder_benchmark.ipynb) | Compare BPOSD, Tesseract and BP variants across noise models, schedules and decoder settings | Collect/resume Sinter CSVs; plot LER curves and the BPOSD–Tesseract gap |
| [Adaptive BP cases](check/notebooks/adaptive.ipynb) | Inspect two fixed weight-2 syndromes, then evaluate Adaptive BP on saved weight-2 cases | Read saved error/correction records; display corrections, LLRs, convergence and failed cases |
| [Matching graph](check/notebooks/matching_graph.ipynb) | Inspect triangular-code circuits, detector matching graphs, Tanner graphs and parity-check matrices | Construct circuits locally; display graphs/matrices and save selected SVG diagrams |
| [Noise-model sanity](check/notebooks/check_noise_model.ipynb) | Compare current BPOSD/Tesseract simulations with two saved reference curves | Read reference tables; collect/resume Sinter statistics; plot fitted LER intervals |

Use conda `qec` and start the notebook kernel from the `Decoders/` repository root. Run the selected notebook's initialization before its analysis. The notebook files remain separate; their detailed procedures, settings and output meanings are described in [decoder_benchmarks/README.md](decoder_benchmarks/README.md) and [check/README.md](check/README.md).

Editable path manifests stay beside the notebooks under each group's `config/`. Their fields are documented in the corresponding configuration README. Shared datasets remain under `experiment/data/<purpose>/<producer>/<stage>/<run>/`; local figures remain under the owning group's `plots/`. A saved CSV is an input/result artifact, not evidence that every current notebook section has been executed with the same settings.
