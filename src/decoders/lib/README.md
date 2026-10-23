# Decoder library dependencies

`outer/` contains external repositories referenced through the root `.gitmodules`. The clean `outer/ldpc` checkout is pinned by the main repository. Local LDPC additions and modifications are ordinary project files under [ldpc_custom/](ldpc_custom/README.md).

Initialize dependencies with `git submodule update --init --recursive`. Follow the custom LDPC README to assemble and install `ldpc_modified` in conda `qec`; downloading the upstream repository alone does not install the custom decoders.

The custom build is generated under `build/ldpc_modified/`, outside the tracked source tree. It combines upstream code with the project-owned additions rather than committing a renamed copy of every upstream file. Python wrappers in `src/decoders/cpp_backend/` continue to import `ldpc_modified`. The independent project implementation in `src/decoders/cython/` is unchanged.

## Installation order

Create the `qec` environment using the root [environment.yml](../../../environment.yml) before building local dependencies. Run commands from `Decoders/` after `conda activate qec`.

1. Initialize submodules with `git submodule update --init --recursive`.
2. Install BPOTF using `python -m pip install --no-user --no-build-isolation ./src/decoders/lib/outer/BPOTF`. The environment supplies its CMake, NumPy/SciPy and stub-generation prerequisites.
3. Build and install the Tesseract Python wheel following [its source-build instructions](outer/tesseract-decoder/README.md#hacking-on-the-python-module-locally). Use the active Python version and an appropriate package version for this checkout. Its Bazel/Bazelisk toolchain is separate from the Conda scientific environment.
4. Prepare and install `ldpc_modified` following [ldpc_custom/README.md](ldpc_custom/README.md).
5. Install the main project with `python -m pip install --no-user --no-build-isolation -e ".[bench,dev]"`.

The main decoder package eagerly imports its C++ wrappers, so BPOTF, Tesseract and `ldpc_modified` are required before a normal `decoders` import. None of them is built automatically by the root package setup. Their distribution requirements in `pyproject.toml` record dependencies; they do not replace these source-build steps.

For the optional BB/LightStim helpers, install `./src/decoders/lib/outer/LightStim` in `qec` with pip before selecting the main project's `bb` extra. BivariateBicycleCodes is a source submodule used by the BB helper, not an additional PyPI distribution to install. Do not install all of LightStim's optional GPU/server/decoder extras unless those workflows are needed.
