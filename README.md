# Decoders

Quantum error correction decoders, noise models and decoding experiments.

## Repository structure

```text
Decoders/
  src/
    codes/                  # code definitions and construction helpers
      __init__.py
      ca_code.py            # triangular-code definitions
      codes.py              # CSS and other code families
      BB_code.py
    noises.py
    noises_ca.py
    decoders/
      cpp_backend/
      cython/
      naive/
      lib/
        outer/              # upstream submodules, including clean ldpc
        ldpc_custom/        # project-owned LDPC extensions and build preparation
  experiment/
    research/               # research themes and stage notebooks
    non_exploratory/        # decoder benchmarks and implementation checks
    misc/                   # schematic notebook
    utils/                  # reusable experiment helpers
    data/                   # shared experiment datasets
    data_packedup/           # retained historical datasets
  tests/
  report/                   # report artifacts
  environment.yml
  pyproject.toml
  setup.py
```

## Entry points

- [Experiment index and data organization](experiment/README.md)
- [Research questions, stages and workflows](experiment/research/README.md)
- [Decoder benchmarks and implementation checks](experiment/non_exploratory/README.md)
- [Shared experiment helpers](experiment/utils/README.md)
- [Shared datasets](experiment/data/README.md)
- [Code definitions and import paths](src/codes/README.md)
- [External decoders and custom LDPC installation](src/decoders/lib/README.md)

Use conda `qec` for Python and pip operations. Notebook kernels use the `Decoders/` repository root as their working directory. Each experiment or stage documents its configuration beside its notebooks under `config/`; there is no root-level configuration reference.

## Environment and installation

From the repository root, initialize the submodules and create the environment:

```bash
git submodule update --init --recursive
conda env create -f environment.yml
conda activate qec
```

`environment.yml` pins the scientific, notebook and build-tool dependencies. It does not build local decoder repositories. Follow [the decoder dependency instructions](src/decoders/lib/README.md) to install BPOTF, Tesseract and the custom `ldpc_modified` package before installing this project:

```bash
conda activate qec
python -m pip install --no-user --no-build-isolation -e ".[bench,dev]"
```

`bench` adds notebook analysis dependencies; `dev` adds development tools. For BB/LightStim code construction, install the local LightStim repository as described in the dependency instructions, then select `.[bench,dev,bb]`. The `bb` extra does not install GPU or server extras.

The project metadata keeps NumPy compatible with beliefmatching and accepts the development Stim/Tesseract versions used by this checkout. Local distributions must be installed before resolving project requirements; an installed package name alone does not reproduce the custom sources. The main package includes `codes`, the project decoder packages and the `noises`/`noises_ca` modules. External source repositories and `ldpc_custom` build inputs remain separate from the main wheel.
