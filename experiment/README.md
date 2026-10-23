# Experiment index

The current directory contains research themes, decoder benchmarks, implementation checks, shared helpers and datasets. Open the relevant local README before selecting a notebook.

## Directory structure

```text
experiment/
  research/
    bp_failure_dynamics/
    failure_regions/
    force_osd_failures/
    gap_mechanisms/
    instability_guided_bp/
    low_weight_failures/
    osd/
    quotient/
    trapping_sets/
  non_exploratory/
    decoder_benchmarks/
    check/
  misc/
    figures.ipynb
  utils/
  data/
  data_packedup/
```

## Research entries

| Theme | Scope and entry layout |
| --- | --- |
| [BP failure dynamics](research/bp_failure_dynamics/README.md) | Message diagnostics, failure mechanisms and logical ambiguity in separate stages |
| [Failure regions](research/failure_regions/README.md) | Case acquisition, failure analysis and belief visualization in separate notebooks |
| [Forced-OSD failures](research/force_osd_failures/README.md) | Independent forced-OSD case-analysis notebook |
| [Gap mechanisms](research/gap_mechanisms/README.md) | Gap exploration, mechanism attribution and PCM fixed-point stages |
| [Instability-guided BP](research/instability_guided_bp/README.md) | Graph-intervention and dynamics-feedback stages |
| [Low-weight failures](research/low_weight_failures/README.md) | Fixed-support acquisition and prior/scaling studies in one notebook |
| [OSD](research/osd/README.md) | Reliability/rank and handoff/gating branches with separate stage notebooks |
| [Quotient](research/quotient/README.md) | Independent QuotientBP and QuotientSplit notebooks |
| [Trapping sets](research/trapping_sets/README.md) | Failure structures, candidates, expansion/contraction and low-weight structures |

See [research/README.md](research/README.md) for scientific dependencies and reading order. Staged themes keep `README.md`, `notebooks/` and `config/` inside each `stages/<stage>/` directory. Standalone themes keep these directly under the theme.

## Benchmarks, checks and schematic

- [Decoder benchmarks](non_exploratory/decoder_benchmarks/README.md): collection/resume, saved-curve analysis and benchmark figures.
- [Implementation checks](non_exploratory/check/README.md): Adaptive BP examples, matching graphs and noise-model reference comparison. The notebooks remain independent.
- [Schematic notebook](misc/figures.ipynb): separate illustrative figures.

## Execution and data sharing

Use conda `qec`, set the notebook kernel working directory to `Decoders/`, and run the selected notebook's initialization first. Stage-local configuration READMEs describe the JSON files and path manifests actually used by each notebook. The quotient notebooks and grouped checks use separate named manifests in their shared `config/` directories.

Generated datasets are described in [data/README.md](data/README.md). Shared output paths preserve purpose, producer, stage and run, with separate `stats/`, `trace/` and effective `config/` locations. [data_packedup/](data_packedup/) retains historical datasets. An output path declared in configuration may not exist until a result is written.

Figures and animations use local destinations named by each manifest; `plots/` directories may be created when needed. Quotient has no separate plots directory. Reusable program code lives in [utils/](utils/README.md). Do not run concurrent writers into the same stage/run.

## Execution state and evidence

Experiment execution is stopped. Research READMEs summarize qualitative observations from saved data and identify unperformed or incomplete evaluations. Historical results do not establish the behavior of later algorithm revisions. Numeric results remain in data files; editable settings remain in local configuration or notebook code.
