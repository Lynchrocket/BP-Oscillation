# Custom LDPC extensions

This directory contains project-owned extensions to a pinned upstream LDPC revision. The upstream source is downloaded as the clean [outer/ldpc](../outer/ldpc/) submodule. The installed Python package remains `ldpc_modified`, independently of the standard `ldpc` package.

```text
ldpc_custom/
  README.md
  JSON.md
  config.json
  prepare.py
  setup.py
  pyproject.toml
  src_cpp/                       # additional decoder implementations
  src_python/ldpc_modified/       # additional bindings and modified upstream bindings
  python_test/                   # retained custom regression checks
```

The local C++ sources cover Adaptive BP/OSD, expanded/structural/region variants, MBBP-LD, quotient decoding and recursive/correlation refiners. The modified upstream BPOSD binding retains `force_osd`. The package initializer and build files retain all existing extension registrations. Project Cython decoders in `src/decoders/cython/` are separate and unchanged; this dependency still uses its own Cython bindings to expose C++ implementations.

## Prepare and install

From the `Decoders/` repository root:

```bash
git submodule update --init --recursive
conda activate qec
python src/decoders/lib/ldpc_custom/prepare.py
python -m pip install --no-user --no-build-isolation --no-deps ./build/ldpc_modified
```

Use the project `qec` environment with its build/runtime dependencies already available. The preparation script checks [config.json](config.json), the exact upstream commit and a clean submodule. It exports tracked upstream files, renames Python/Cython identifiers and C++ namespaces, then overlays local sources. Header include filenames are preserved. The upstream license, README and bundled headers come from the fixed revision.

All assembled sources, generated type stubs, Cython C++, binaries and build metadata belong in the ignored `build/ldpc_modified/` directory. Modify the sources here, then prepare and install again. Do not edit generated copies. Preparation replaces its previously marked build directory; do not run it concurrently with compilation or while using an editable installation from that directory. Install again after every preparation. No source-control history or previously generated binaries are needed to reproduce the assembly.

## Configuration and validation

[JSON.md](JSON.md) documents both build-control fields. There are no experiment settings in this configuration. Decoder arguments, stopping behavior and existing `ldpc_modified` imports remain unchanged.

Run the retained custom checks after installation:

```bash
conda activate qec
PYTHONPATH=src:. python -m pytest src/decoders/lib/ldpc_custom/python_test -q
```

These are dependency checks, not experiment runs. Existing test failures must be compared against the pre-migration baseline before attributing them to this source-layout change. Preserve existing copyright and license notices when extending or redistributing this code.
