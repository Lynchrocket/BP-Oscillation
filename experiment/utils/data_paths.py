# bp_failure_dynamics, failure_regions, force_osd_failures, gap_mechanisms, instability_guided_bp, low_weight_failures, osd, quotient, trapping_sets
"""Purpose-classified datasets with experiment-local configuration and plots."""
import json
from pathlib import Path
from string import Formatter
from tempfile import TemporaryDirectory
PURPOSES = {
    'decoder_benchmark': 'benchmark', 
    'failure_regions': 'failure_analysis', 
    'low_weight_failures': 'failure_analysis', 
    'force_osd_failures': 'failure_analysis', 
    'bp_failure_dynamics': 'bp_dynamics', 
    'subgraph_structure': 'graph_analysis', 
    'trapping_sets': 'graph_analysis', 
    'ca_matching_graphs': 'graph_analysis', 
    'instability_guided_bp': 'decoder_intervention', 
    'adaptive_bp_repair': 'decoder_intervention', 
    'osd_reliability': 'osd', 
    'osd_handoff': 'osd', 
    'gap_mechanisms': 'decoder_gap', 
    'quotient_bp': 'quotient_decoding', 
    'quotient_split': 'quotient_decoding', 
    'bb_quotient_split': 'quotient_decoding', 
    'noise_model_sanity_check': 'noise_validation',
}
KINDS = {'stats', 'trace', 'config'}

def contained_path(root: Path, relative: str) -> Path:
    """Reject absolute paths, parent traversal and escaping symlinks."""
    value = Path(relative)
    assert relative and (not value.is_absolute()) and ('..' not in value.parts), relative
    target = (root / value).resolve()
    assert target.is_relative_to(root.resolve()), target
    return target

class DataCatalog:
    """Resolve an explicitly selected dataset version; never choose a fallback."""

    def __init__(self, data_root: Path):
        self.root = Path(data_root).resolve()
        content = json.loads((self.root / 'catalog.json').read_text())
        assert type(content['schema_version']) is int and content['schema_version'] == 1
        self.datasets = content['datasets']
        assert isinstance(self.datasets, dict)
        for key, item in self.datasets.items():
            assert isinstance(key, str) and key
            assert item['kind'] in KINDS, item
            assert item['producer'] in PURPOSES, item
            assert all((isinstance(item[field], str) and item[field] for field in ('stage', 'run', 'path'))), item
            contained_path(self.root, item['path'])

    def resolve(self, dataset_id: str, *, kind: str) -> Path:
        assert kind in KINDS and dataset_id in self.datasets, dataset_id
        item = self.datasets[dataset_id]
        assert item['kind'] == kind, (dataset_id, kind, item['kind'])
        target = contained_path(self.root, item['path'])
        assert target.exists(), (dataset_id, target)
        return target

class ExperimentPaths:
    """Read explicit repository-relative input/output paths from local JSON."""

    def __init__(self, project_root: Path, manifest_path: Path):
        self.project_root = Path(project_root).resolve()
        self.manifest_path = Path(manifest_path).resolve()
        self.manifest = json.loads(self.manifest_path.read_text())
        assert type(self.manifest['schema_version']) is int and self.manifest['schema_version'] == 1
        for section in ('inputs', 'outputs', 'templates'):
            assert isinstance(self.manifest[section], dict), section
            for key, value in self.manifest[section].items():
                assert isinstance(key, str) and key and isinstance(value, str), (key, value)
                contained_path(self.project_root, value)
        assert {'run', 'stats', 'trace', 'config', 'plots'} <= self.manifest['outputs'].keys()
        self._temporary = None

    def output(self, name: str, /, **parameters) -> Path:
        section = 'templates' if parameters else 'outputs'
        assert name in self.manifest[section], name
        value = self.manifest[section][name]
        required = {field for _, field, _, _ in Formatter().parse(value) if field is not None}
        assert required == parameters.keys(), (name, required, parameters.keys())
        target = contained_path(self.project_root, value.format(**parameters))
        if self._temporary is not None:
            run = contained_path(self.project_root, self.manifest['outputs']['run'])
            plots = contained_path(self.project_root, self.manifest['outputs']['plots'])
            if target.is_relative_to(run):
                target = Path(self._temporary.name) / 'data' / target.relative_to(run)
            elif target.is_relative_to(plots):
                target = Path(self._temporary.name) / 'plots' / target.relative_to(plots)
            else:
                assert False, ('Output is outside the declared run/plots roots', name, target)
        return target

    @property
    def base(self) -> Path:
        """Explicit parent of the configured run, for named historical runs."""
        return self.output('run').parent

    def run(self, run: str='main', *, smoke_mode: bool=False) -> Path:
        assert type(smoke_mode) is bool, smoke_mode
        if smoke_mode and self._temporary is None:
            self._temporary = TemporaryDirectory(prefix='decoders-smoke-')
        assert run and len(Path(run).parts) == 1 and (run not in {'.', '..'}), run
        target = self.output('run')
        return target if run == 'main' else contained_path(target.parent, run)

    def stats(self, run: str='main') -> Path:
        return self.output('stats') if run == 'main' else self.run(run) / 'stats'

    def trace(self, run: str='main') -> Path:
        return self.output('trace') if run == 'main' else self.run(run) / 'trace'

    def plots(self, run: str='main') -> Path:
        assert run and len(Path(run).parts) == 1 and (run not in {'.', '..'}), run
        target = self.output('plots')
        return target if run == 'main' else contained_path(target.parent, run)

    def input(self, name: str, /, **parameters) -> Path:
        """Resolve a named input path and assert it exists."""
        assert name in self.manifest['inputs'], name
        value = self.manifest['inputs'][name]
        required = {field for _, field, _, _ in Formatter().parse(value) if field is not None}
        assert required == parameters.keys(), (name, required, parameters.keys())
        target = contained_path(self.project_root, value.format(**parameters))
        assert target.exists(), (name, target)
        return target

    def close(self):
        """Delete temporary smoke outputs; main outputs are never removed."""
        if self._temporary is not None:
            self._temporary.cleanup()
            self._temporary = None
