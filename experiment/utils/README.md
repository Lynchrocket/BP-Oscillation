# Shared experiment utilities

Program code and purpose-classified datasets are shared; editable configuration and figures belong to each experiment stage. Helpers are ordinary functions/classes receiving explicit configuration, cache and callback arguments. They do not read notebook namespaces or retain experiment state in module globals.

## Modules

| Module | Purpose | Research users |
| --- | --- | --- |
| [__init__.py](__init__.py) | Package initialization. | `bp_failure_dynamics`, `failure_regions`, `force_osd_failures`, `gap_mechanisms`, `instability_guided_bp`, `low_weight_failures`, `osd`, `quotient`, `trapping_sets` |
| [bb_model.py](bb_model.py) | Build vendored BB noise models with explicit parameters. | — |
| [belief_tracking.py](belief_tracking.py) | Belief trajectories and decoding visualizations. | `failure_regions`, `force_osd_failures`, `low_weight_failures` |
| [benchmark_plotting.py](benchmark_plotting.py) | Benchmark tables, gap intervals and plotting helpers. | `failure_regions`, `force_osd_failures` |
| [bp_message_dynamics.py](bp_message_dynamics.py) | Native BP message trajectories and dynamics metrics. | `bp_failure_dynamics`, `instability_guided_bp`, `quotient`, `trapping_sets` |
| [bp_trace.py](bp_trace.py) | Shared native BP contexts, trace execution and backend cache markers. | `bp_failure_dynamics`, `gap_mechanisms`, `instability_guided_bp`, `osd`, `quotient`, `trapping_sets` |
| [data_paths.py](data_paths.py) | Dataset catalog resolution and stage input/output paths. | `bp_failure_dynamics`, `failure_regions`, `force_osd_failures`, `gap_mechanisms`, `instability_guided_bp`, `low_weight_failures`, `osd`, `quotient`, `trapping_sets` |
| [decoder_setup.py](decoder_setup.py) | Construct BP, BPOSD and force-OSD decoders. | `failure_regions`, `force_osd_failures`, `low_weight_failures`, `quotient`, `trapping_sets` |
| [dem_graph.py](dem_graph.py) | Parse detector error models, construct Tanner graphs and extract coordinates or active detector IDs. | `failure_regions`, `trapping_sets` |
| [gap_mechanisms.py](gap_mechanisms.py) | Gap attribution, OSD candidates and sector-cost analysis. | `gap_mechanisms`, `osd` |
| [graph_intervention.py](graph_intervention.py) | Redundant-check interventions and conditional outcome metrics. | `instability_guided_bp` |
| [graph_struct.py](graph_struct.py) | Tanner graph construction and structural features. | `failure_regions`, `force_osd_failures`, `low_weight_failures`, `quotient`, `trapping_sets` |
| [graph_vis.py](graph_vis.py) | Graph visualizations and fault overlays. | `failure_regions`, `force_osd_failures`, `low_weight_failures`, `quotient`, `trapping_sets` |
| [misc.py](misc.py) | Bit conversion, numerical helpers, table I/O, aggregate paths and algorithm cache version checks. | `bp_failure_dynamics`, `failure_regions`, `force_osd_failures`, `gap_mechanisms`, `instability_guided_bp`, `low_weight_failures`, `osd`, `quotient`, `trapping_sets` |
| [osd_reliability.py](osd_reliability.py) | Reliability ordering, endpoint analysis and OSD comparisons. | `osd` |
| [prior_scaling.py](prior_scaling.py) | Prior scaling and low-weight failure analysis. | `low_weight_failures` |
| [quotient_split.py](quotient_split.py) | Quotient-split graph and decoding helpers. | `quotient` |
| [sample.py](sample.py) | Sample collection and decoder outcome counting. | `failure_regions`, `force_osd_failures`, `low_weight_failures` |
| [tess_viz/to_json.py](tess_viz/to_json.py) | Export graph visualization data. | — |
| [trajectory_quotient_bp.py](trajectory_quotient_bp.py) | Trajectory-based quotient BP calculations. | `quotient` |
| [trapping_sets.py](trapping_sets.py) | Trapping-set candidates, fault conversion, structural collection and aggregate outputs. | `trapping_sets` |

Module header comments list the first-level research folders that use the module, directly or through another helper. Modules without research consumers have no research folder list. Type annotations use explicit imports or forward references rather than postponed annotation evaluation.

Variant suffixes preserve pre-existing differences in parsing and data schemas; they are not interchangeable. `sample` imports decoder constructors from `decoder_setup`. Common I/O helpers live in `misc`, and trapping-set candidate and collection helpers live together in `trapping_sets`.

## Ordinary notebook execution

Notebooks contain their scientific cell sources directly. Shared helpers are ordinary Python functions and classes with explicit configuration, cache, decoder and callback arguments. Notebook adapters expose ordinary argument signatures; no helper reads a notebook namespace. Merged sections use distinct variable names.

## Shared data and sweeps

[data_paths.py](data_paths.py) supplies explicit catalog resolution and classified output paths. Figures remain local. Parameter sweeps, benchmark task construction and plotting loops are written directly in their owning notebooks. These loops use notebook variables directly; they do not pass notebook dependencies into wrappers or return dictionaries of local variables. Exhaustive enumeration and random sampling retain separate paths. Message-specific metric factories preserve the richer directional schema. The Plotly graph helper has one implementation in [graph_vis.py](graph_vis.py), including the fault overlays requested by trapping-set analyses.

The grouped OSD studies read explicit input/output paths from their stage-local `config/paths.json`. Output ownership is no longer inferred from directory names.

- [trapping_sets.py](trapping_sets.py): shared weight-specific structural collection, separate outcome/structure summaries, local Plotly output and expansion/contraction record summaries. Existing graph feature definitions are reused.

[bb_model.py](bb_model.py) exposes `build_bb_model(...)` with explicit preset, noise, cycle count, output filename and print level. It preserves the vendored BB model calculation without runtime source rewriting.

Import functions and classes under their original names when unambiguous. Use conventional module abbreviations such as `np`, `pd` and `plt` without experiment prefixes. Aliases remain only when an existing binding, a local scope or a different implementation needs to be distinguished. Import names shared across sections refer to the same module/object; scientific configuration and mutable state keep their section-specific names.

Notebook kernels run from the `Decoders/` repository root. Their initialization uses `Path.cwd()` directly and adds the three project import locations in one statement; it does not search ancestors or validate the launch directory.

## Native BP research trajectories

`BP.from_pcm(H, priors)` in `src/decoders/cython/bp.py` prepares a native
parity-check graph. `BP.trace(syndrome, iterations, ...)` performs the message
updates in `bp_core.pyx`. The BP message-dynamics, graph-intervention, OSD and gap
helpers only adapt parameters and trace schemas. Their contexts inherit
`bp_trace.BPContext`, whose cached decoder assumes that `H` and `priors` remain
unchanged. Create a new context for a changed graph or channel.

A full trace contains `c2v`, `v2c`, `llr`, `hard`, `valid`,
`runtime_seconds` and `first_valid`. Index zero is the initial state; edges use
`np.nonzero(H)` row-major order. Every call owns its message/history arrays, so
later calls do not overwrite earlier traces. `record_messages=False` omits the
edge histories (their shape is `(0, num_edges)`) while retaining node histories.
The existing experiment trace classes and CSV columns are unchanged.

| API option | Default | Meaning |
| --- | --- | --- |
| `iterations` | decoder `max_iter` | Positive number of complete sweeps. |
| `update_rule` | `min_sum` | `min_sum` or `sum_product`. Research sum-product retains `1e-15` clipping; normal product-sum decode retains `1e-12`. |
| `scaling`, `offset` | `1.0`, `0.0` | Positive scalar/per-edge min-sum scaling and non-negative offset; neither modifies sum-product. |
| `damping` | `0.0` | Scalar/per-edge C2V damping in `[0, 1)`; V2C uses the updated extrinsic LLR directly. |
| `schedule` | `parallel` | `parallel`, check-ordered `serial`, `random_serial_fixed`, or `random_serial`; the latter reshuffles a fresh canonical check order each sweep. |
| `rng` | generator seeded with zero | NumPy generator used for check order and initial perturbations. |
| `prior_llr` | channel LLR | Optional finite length-`n` prior override. |
| `initial_message_sigma` | `0.0` | Non-negative Gaussian V2C initialization noise. |
| `stop_early` | `False` | Stop after a completed sweep satisfies the syndrome; otherwise continue through the iteration budget. |
| `record_messages` | `True` | Record both edge directions as well as node states. |
| `include_initial_valid` | `True` | Allow `first_valid=0`; the handoff adapters retain their original positive-iteration convention. |
| `accumulation`, `serial_accumulation` | `sequential` | `sequential`, `sum_then_prior`, or `numpy`; preserve each experiment's floating-point reduction order. |
| `combine_scaling_damping` | `False` | Preserve the OSD/handoff coefficient order `(1-damping)*scaling*raw`; requires zero offset. |

Normal `decode_to_correction`, BP/BPOSD stopping rules, and the two-item legacy
kernel return values retain their previous behavior. Quotient research code and
its partition-function BP have not been changed.

Timing now measures native per-shot initialization, history allocation and
sweeps. It excludes graph preparation and Python argument preparation. Compare
implementations with the same iteration count, stopping rule and recorded
histories. Historical Python timings are not interchangeable with native timings.

`bp_trace.check_bp_cache` records `config/bp_backend.json` beneath the existing
manifest's configuration output directory. Its `bp_backend` value is
`cython_research_bp_v1`; the same generated metadata key is added to the effective
configuration. This is implementation metadata,
not a new user configuration option. Existing non-empty stats/trace directories
without a matching marker are rejected before reuse or overwrite, including when
`reuse_existing_data=False`. Select fresh output locations to run the native
experiment. Existing results are preserved; the migration does not rerun them.
Smoke paths must be activated before checking the backend and remain temporary.

## Benchmark-compatible OSD and result versions

`ordered_statistics_postprocess(ordering_values, syndrome, H, order, priors, osd_method="osd_cs")` in the native Cython OSD module requires channel priors separately from ordering values. It matches the benchmark signed-value column ordering and `sum(log(1/p_i))` candidate cost. Order zero tests the initial solution; a positive order tests all single non-pivot flips and pairs within the first `min(order,k)` non-pivot positions. `osd_candidate_count(k,order)` reports that budget. The compiled `OSD_VERSION` is `cpp_cs_v1`. Unsupported methods and invalid inputs raise errors; there is no unweighted fallback. Deliberately transformed reliability scores remain experimental ordering variants.

`misc.check_algorithm_cache(paths, config, versions)` records implementation-owned `algorithm_versions` in the effective config and `config/algorithm_versions.json`. Nonempty results with missing or mismatched versions are rejected before reuse. Existing datasets are retained; revised stages use their manifest's separate run paths.

The message-dynamics trace adapter accepts `rng` and `initial_message_sigma` as keyword arguments. Their default behavior is unchanged (`sigma=0`); randomized restart controls supply reproducible independent seeds without changing the channel prior.

## Configuration snapshots and cache reuse

Prior-scaling and the OSD/gap stage initializers record plain configuration snapshots without computing or comparing SHA256 fingerprints. Existing snapshots are retained when `reuse_existing_data=true`; a missing snapshot or `reuse_existing_data=false` writes the current `config` object. Rank follow-up additionally records the full `base_config` object. Legacy snapshots may retain unused fingerprint fields; they do not need migration. Prior-scaling no longer hashes its implementation, native binary or imported source classes.

PCM E1 completion markers now contain `complete: true`. Readers accept legacy markers carrying additional fingerprint fields while retaining trace-file and progress-completeness checks. Backend and algorithm-version markers remain independently enforced. SHA256 uses that determine reproducible dataset splits, seeds or audit identifiers are unchanged. No new user configuration key controls this behavior.
