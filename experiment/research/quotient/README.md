# Quotient BP and quotient-splitting explorations

Investigate how quotient representations, damping and selective freezing change decoding trajectories on saved low-weight cases. QuotientBP and QuotientSplit remain two independent notebooks in one unchanged folder layout.

## Structure and entries

```text
quotient/
  README.md
  notebooks/
    quotient.ipynb
    quotient_split.ipynb
  config/
    quotient_bp_paths.json
    quotient_split_paths.json
```

| Entry | What it does |
| --- | --- |
| [quotient.ipynb](notebooks/quotient.ipynb) | Build quotient models, inspect hard-decision trajectories and damping, compare with the C++ decoder path, then try reliability-frozen quotient decoding |
| [quotient_split.ipynb](notebooks/quotient_split.ipynb) | Inspect low-weight quotient-splitting trajectories/damping and compare the existing C++ QuotientSplit implementation |

Use conda `qec` and working directory `Decoders/`. Initialize only the notebook you intend to run; its manifest is independent of the other notebook. [config/README.md](config/README.md) describes paths and the notebook-owned scientific settings.

## Experiment flow

1. Load the shared weight-2/weight-3 tables from the low-weight producer, preserving bit strings. The examples use the fixed code and one-round setup declared in the notebook.
2. Reconstruct the circuit/DEM, PCM and priors at the selected p; identify the saved error, syndrome and true logical sector.
3. Inspect quotient and original hard-decision trajectories and the notebook's damping variants. The custom trajectory helpers here are intentional experimental implementations, not a claim that every update runs through native Cython BP.
4. Run the notebook's C++ backend comparison. In the QuotientBP frozen section, the explicit example uses `h_max_iter=50`, `d_max_iter=60`, zero h/d damping and `max_unfrozen_variables=20`, then compares the correction's logical observables with truth.

## Current conclusion and evidence limit

The notebooks provide exploratory trajectories and case diagnostics, but no saved aggregate comparison establishes an overall LER or runtime improvement. An individual trajectory or a frozen-case loop does not replace a paired physical-noise benchmark.

Aggregate performance evaluation: 该部分实验未进行，可参考data_packedup自行分析

## Inputs and outputs

Both notebooks read `experiment/data/failure_analysis/low_weight_failures/root/main/trace/error_weight/`. Their output roots remain separate:

- `experiment/data/quotient_decoding/quotient_bp/root/main/`
- `experiment/data/quotient_decoding/quotient_split/root/main/`

`stats`, `trace` and `config` in the manifests are reserved output locations, not proof that files exist. The current notebooks do not save figures and there is no plots subfolder. The required but unused `outputs.plots` manifest field points to the quotient folder itself.
