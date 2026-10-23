# Tesseract decoding visualization

This directory contains the log converter, browser viewer and an example visualization dataset:

```text
tess_viz/
  README.md
  to_json.py
  index.html
  logfile.json
```

## Convert a decoder log

Build and invoke the decoder using the instructions in [the Tesseract dependency](../../../src/decoders/lib/outer/tesseract-decoder/README.md). Enable `--verbose` to include coordinates, and retain its `Error`, `Detector`, `activated_errors` and `activated_detectors` records in a text log. This directory is a viewer, not a separate decoder build tree.

From the `Decoders/` repository root, activate conda `qec` and convert a saved log:

```bash
conda activate qec
python experiment/utils/tess_viz/to_json.py decoder_log.txt -o experiment/utils/tess_viz/logfile.json
```

`decoder_log.txt` is the input log supplied by the user. The explicit output path writes the JSON beside the viewer. Without `-o`, the converter uses its default output filename in the current working directory.

## View the result

Open [index.html](index.html) in a browser. The viewer attempts to load [logfile.json](logfile.json) beside it; its file picker can load another JSON produced by [to_json.py](to_json.py).
