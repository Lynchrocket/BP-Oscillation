# bp_failure_dynamics, gap_mechanisms, instability_guided_bp, osd, quotient, trapping_sets
"""Native BP graph caching and provenance for research notebooks."""
from functools import cached_property
import json

from decoders.cython.bp import BP, BP_TRACE_VERSION


class BPContext:
    """Cache a native graph for a context whose H and priors stay fixed."""

    @cached_property
    def bp_decoder(self):
        return BP.from_pcm(self.H, self.priors)


def check_bp_cache(data_paths, config):
    """Reject results timed with a different BP implementation before reuse."""
    marker = data_paths.output('config') / 'bp_backend.json'
    has_results = any(
        folder.exists() and any(folder.iterdir())
        for folder in (data_paths.output('stats'), data_paths.output('trace'))
    )
    if has_results:
        assert marker.is_file(), (
            'Historical Python BP results have no native backend marker. '
            'Use a fresh output location; historical results are not overwritten.'
        )
        assert json.loads(marker.read_text())['bp_backend'] == BP_TRACE_VERSION, 'BP backend cache mismatch'
    config['bp_backend'] = BP_TRACE_VERSION
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({'bp_backend': BP_TRACE_VERSION}, indent=2) + '\n')
