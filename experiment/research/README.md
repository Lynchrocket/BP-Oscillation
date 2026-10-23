# Research experiments: navigation and execution flow

Start with the research question, then open the theme README and the relevant notebook/stage. This directory contains case acquisition, diagnostic/mechanism studies and experimental decoder policies; their outputs are not interchangeable benchmarks.

## Where to go

| Theme | Question / role | Current layout |
| --- | --- | --- |
| [low_weight_failures](low_weight_failures/README.md) | Acquire fixed-weight cases; compare prior/scaling sensitivity | One notebook with independent scientific parts |
| [failure_regions](failure_regions/README.md) | Acquire physical failure/disagreement cases and visualize selected trajectories | Three notebooks sharing explicit paths |
| [force_osd_failures](force_osd_failures/README.md) | Test whether disagreements remain with forced OSD | One notebook |
| [bp_failure_dynamics](bp_failure_dynamics/README.md) | Separate instability, syndrome consistency and logical correctness; classify/repair failures | Three stage notebooks |
| [instability_guided_bp](instability_guided_bp/README.md) | Test graph interventions and whether probe gain predicts useful feedback | Two stage notebooks |
| [osd](osd/README.md) | Reliability/rank branch and causal handoff/gating branch | Four stage notebooks in two research branches |
| [gap_mechanisms](gap_mechanisms/README.md) | Locate paired gaps, attribute mechanisms, separate PCM generation from selection | Three stage notebooks |
| [trapping_sets](trapping_sets/README.md) | Empirical motifs, candidate constructions and fixed-weight structures | Four stage notebooks |
| [quotient](quotient/README.md) | Explore quotient trajectories, splitting, damping and freezing | Two independent notebooks; no stage reorganization |

## Directory convention

For staged themes:

```text
<theme>/
  README.md                 # research flow and cross-stage interpretation
  stages/<stage>/
    README.md               # purpose, steps, observed results, limits
    notebooks/              # independent entry and initialization
    config/                 # local field/path reference and JSON files
    plots/                  # generated local figures, where configured
```

The standalone themes keep `notebooks/`, `config/` and README directly under the theme. Quotient has no plots subfolder. A configured output directory may be absent until its first result is written.

## Data dependencies and practical run order

1. Prepare the input producer needed by the chosen study: low-weight acquisition supplies shared fixed-support tables; failure-region acquisition supplies physical disagreement/failure cases.
2. BP failure dynamics progresses from diagnostics to taxonomy/repair and then logical-sector ambiguity. Instability-guided BP proceeds from interventions to feedback validation, but its follow-up reads a named historical snapshot.
3. The OSD reliability branch proceeds to rank/snapshot follow-up, while online handoff proceeds to intervention gating. The branch-local `01`/`02` prefixes do not define one linear four-stage pipeline.
4. Gap exploration, attribution and PCM fixed-point analysis form a scientific progression; each generates its own populations. Trapping-set stages use the input producer/snapshot recorded in each manifest. Failure-region belief visualization waits for the structural feature files it explicitly reads.

This is a dependency/reading guide, not an instruction to run all notebooks in one kernel. Use conda `qec`, launch from `Decoders/`, and run each chosen notebook's initialization before its scientific sections. Preserve input versions, seeds and the selected noise/cohort settings when comparing results.

## Reading the conclusions

Every stage README separates its research question from saved evidence and states the scope and limitations of its qualitative conclusions. Experiments are stopped. Sections without available result data explicitly identify the missing evaluation; historical tables do not validate later implementation changes.

`experiment/data/<purpose>/<producer>/<stage>/<run>/stats/` holds aggregates; `trace/` holds detailed records/caches, and `config/` holds generated settings. Figures remain local. Fixed-weight fractions, enriched-failure repair rates, oracle ceilings and unbiased physical-noise LER answer different questions. Syndrome consistency alone is not logical success; example plots and historical outputs are not proof of a universal mechanism.

Each local `config/README.md` describes current files, defaults, path keys and initializer constraints. No root-level JSON document is needed. Historical source names in imports or shared producer paths remain where required for compatibility; they are not additional experiments.
