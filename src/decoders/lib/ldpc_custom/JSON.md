# LDPC build-control JSON

The local [config.json](config.json) is read by [prepare.py](prepare.py). Both keys are required and additional keys are rejected. Paths are fixed relative to this directory: upstream `../outer/ldpc`, custom sources here, and output `<repository>/build/ldpc_modified`.

| Key | Default | Allowed values and meaning | Dependencies / conflicts |
| --- | --- | --- | --- |
| `upstream_commit` | `d3429964cd4ffe1abfc041c6ec8b8425cb174f40` | A full lowercase hexadecimal Git commit ID. Selects the exact upstream revision exported into the build. | Must equal the checked-out upstream HEAD and the submodule must be clean. Changing it requires reviewing local-source compatibility and updating the parent repository's gitlink. No fallback or automatic checkout is performed. |
| `print_level` | `2` | Integer `1`: no preparation output; `2`: progress bars and final output path; `3`: progress bars, individual file operations and final path. Booleans are rejected. | Controls preparation messages only. Git errors still propagate; it does not change build inputs, algorithms or experiment data. |

The generated build-owner marker is implementation metadata, not an additional user configuration file. Its owner path prevents preparation from replacing an unrelated directory. An unmarked or invalid build directory must be inspected and moved aside explicitly. Preparation and compilation are sequential operations; concurrent writers to the same build directory are unsupported.
