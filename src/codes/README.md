# Code definitions

```text
codes/
  README.md
  __init__.py
  ca_code.py
  codes.py
  BB_code.py
```

[ca_code.py](ca_code.py) defines the fixed triangular codes and the recursive triangular-code family. [codes.py](codes.py) defines the CSS wrapper and the Iceberg, Steane, surface, toric and bivariate-bicycle families. [BB_code.py](BB_code.py) contains the existing LightStim-backed BB construction helpers and requires that dependency when imported.

Use the `src/` import root, as configured for the project and notebook initialization:

```python
from codes.ca_code import CA_5_2_3, CA_9_2_5, TriangularCode
from codes.codes import CSS, SteaneCode
```

Generic code classes are also exported by the package, so `from codes import CSS` remains available. The old top-level `ca_code` module is not retained. Both `ca_code` and `src.ca_code` imports must use `codes.ca_code` instead.

Moving a definition changes its module path. Existing serialized objects referring to old modules may require migration, and source/module-based experiment fingerprints can reject historical caches. Do not edit saved fingerprints to make an old cache appear compatible; retained results remain historical observations. This directory change does not rerun experiments or modify code matrices.
