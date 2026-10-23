# low_weight_failures
"""Fixed-support comparisons of distance, physical noise and min-sum scaling."""
from dataclasses import dataclass
from itertools import combinations
from math import comb
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from beliefmatching import detector_error_model_to_check_matrices
from tqdm.auto import tqdm

from codes.ca_code import CA_9_2_5, TriangularCode
from src.noises_ca import CircuitBuilder, ErrorRate
from src.decoders.cpp_backend.bp_osd import BPOSD
from src.decoders.cpp_backend.tesseract import Tesseract
from src.decoders.cython.bp import BP, BP_TRACE_VERSION


@dataclass
class Context:
    code: object
    dem: object
    H: np.ndarray
    L: np.ndarray
    priors: np.ndarray
    bp: BP


def validate_config(config):
    assert type(config['schema_version']) is int and config['schema_version'] == 1
    assert type(config['print_level']) is int and config['print_level'] in [1, 2, 3]
    assert type(config['reuse_existing_data']) is bool
    assert type(config['random_seed']) is int and config['random_seed'] >= 0
    assert type(config['max_iter']) is int and config['max_iter'] > 0
    assert type(config['d_rounds']) is int and config['d_rounds'] > 0
    assert 0 <= config['p_meas'] < 0.5
    assert type(config['osd_order']) is int and config['osd_order'] >= 0
    assert config['osd_method'] == 'osd_cs'
    assert type(config['tesseract_beam']) is int and config['tesseract_beam'] > 0
    assert type(config['tesseract_pq_limit']) is int and config['tesseract_pq_limit'] > 0
    grid = config['alpha_grid']
    assert 0 < grid['start'] < grid['stop']
    assert type(grid['count']) is int and grid['count'] >= 2
    alphas = alpha_values(config)
    assert np.isfinite(alphas).all() and np.all(alphas > 0)
    assert config['low_p_alphas'] and np.all(np.isfinite(config['low_p_alphas']))
    assert np.all(np.asarray(config['low_p_alphas']) > 0)
    assert len(set(config['low_p_alphas'])) == len(config['low_p_alphas'])
    for grid in config['p_grids'].values():
        assert 0 < grid['start'] < grid['stop'] < 0.5
        assert type(grid['count']) is int and grid['count'] >= 2
    assert set(config['p_grids']) == {'high_p', 'low_p'}
    groups = set()
    for case in config['cases']:
        key = (case['distance'], case['weight'])
        assert type(case['distance']) is int and case['distance'] in [5, 9]
        assert key not in groups
        assert type(case['weight']) is int and case['weight'] > 0
        assert case['sampling'] in ['enumerate', 'sample']
        if case['sampling'] == 'enumerate':
            assert case['sample_count'] is None
        else:
            assert type(case['sample_count']) is int and case['sample_count'] > 0
        groups.add(key)
    assert groups


def alpha_values(config):
    grid = config['alpha_grid']
    return np.unique(np.r_[np.linspace(grid['start'], grid['stop'], grid['count']), grid['extra']])


def build_context(distance, p, config):
    code = CA_9_2_5() if distance == 5 else TriangularCode(3)
    assert code.d == distance
    circuit = CircuitBuilder(code).get_circuit_level_circuit(ErrorRate(p=float(p), p_meas=config['p_meas']), d_rounds=config['d_rounds'])
    dem = circuit.detector_error_model()
    matrices = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
    H = matrices.check_matrix.toarray().astype(np.uint8)
    L = matrices.observables_matrix.toarray().astype(np.uint8)
    priors = np.asarray(matrices.priors, dtype=float)
    return Context(code, dem, H, L, priors, BP.from_pcm(H, priors))


def make_supports(num_errors, case, config):
    weight = case['weight']
    assert weight <= num_errors
    total = comb(num_errors, weight)
    count = total if case['sampling'] == 'enumerate' else case['sample_count']
    assert count <= total
    visible = config['print_level'] > 1
    if case['sampling'] == 'enumerate':
        supports = list(tqdm(combinations(range(num_errors), weight), total=total, desc=f"Supports d={case['distance']} w={weight}", disable=not visible))
    else:
        rng = np.random.default_rng(np.random.SeedSequence([config['random_seed'], case['distance'], weight]))
        seen, supports = set(), []
        with tqdm(total=count, desc=f"Supports d={case['distance']} w={weight}", disable=not visible) as progress:
            while len(supports) < count:
                support = tuple(sorted(rng.choice(num_errors, weight, replace=False).tolist()))
                if support not in seen:
                    seen.add(support)
                    supports.append(support)
                    progress.update()
    errors = np.zeros((count, num_errors), dtype=np.uint8)
    errors[np.arange(count)[:, None], np.asarray(supports)] = 1
    return errors, total


def bit_text(bits):
    return ''.join(map(str, np.asarray(bits, dtype=np.uint8).tolist()))


def outcome(context, correction, syndrome, truth):
    consistent = bool(np.array_equal(context.H @ correction % 2, syndrome))
    logical_match = bool(np.array_equal(context.L @ correction % 2, truth))
    return consistent, logical_match, consistent and logical_match


def tesseract_outcomes(context, syndromes, truths, p, config):
    compiled = Tesseract(context.code, float(p), beam=config['tesseract_beam'], pq_limit=config['tesseract_pq_limit'], seed=config['random_seed']).compile_decoder_for_dem(dem=context.dem)
    results = []
    for syndrome, truth in tqdm(
        zip(syndromes, truths), total=len(syndromes), desc='Tesseract reference',
        leave=False, disable=config['print_level'] < 3,
    ):
        correction = np.zeros(context.H.shape[1], dtype=np.uint8)
        correction[list(compiled.decoder.decode_to_errors(syndrome))] = 1
        results.append(outcome(context, correction, syndrome, truth))
    return np.asarray(results, dtype=bool)


def evaluate_point(context, errors, syndromes, truths, reference, p, alpha, config):
    decoder = BPOSD(
        context.code, float(p), max_iter=config['max_iter'],
        ms_scaling_factor=float(alpha), osd_order=config['osd_order'],
        osd_method=config['osd_method'],
    ).compile_decoder_for_dem(dem=context.dem).decoder
    rows = []
    for sample_id in tqdm(range(len(errors)), desc='Fixed supports', leave=False, disable=config['print_level'] < 3):
        syndrome, truth = syndromes[sample_id], truths[sample_id]
        trace = context.bp.trace(
            syndrome, config['max_iter'], scaling=float(alpha), stop_early=True,
            record_messages=False, accumulation='sequential',
        )
        endpoint = trace.first_valid if trace.first_valid is not None else len(trace.hard) - 1
        bp_correction = trace.hard[endpoint]
        bposd_correction = decoder.decode_to_correction(syndrome)
        row = {
            'sample_id': sample_id,
            'p': float(p),
            'alpha': float(alpha),
            'error': bit_text(errors[sample_id]),
            'syndrome': bit_text(syndrome),
            'truth': bit_text(truth),
            'bp_correction': bit_text(bp_correction),
            'bposd_correction': bit_text(bposd_correction),
            'bp_iterations': endpoint,
        }
        outcomes = {
            'bp': outcome(context, bp_correction, syndrome, truth),
            'bposd': outcome(context, bposd_correction, syndrome, truth),
            'tesseract': reference[sample_id],
        }
        for method, values in outcomes.items():
            for field, value in zip(['syndrome_consistent', 'logical_match', 'success'], values):
                row[f'{method}_{field}'] = bool(value)
        rows.append(row)
    return pd.DataFrame(rows)


def save_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.csv.tmp')
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def summarize_point(frame):
    result = {'samples': len(frame)}
    for method in ['bp', 'bposd', 'tesseract']:
        failures = int((~frame[f'{method}_success']).sum())
        result[f'{method}_failures'] = failures
        result[f'{method}_failure_rate'] = failures / len(frame)
        result[f'{method}_syndrome_inconsistent'] = int((~frame[f'{method}_syndrome_consistent']).sum())
        result[f'{method}_logical_mismatches'] = int((~frame[f'{method}_logical_match']).sum())
    gap = ~frame.bposd_success & frame.tesseract_success
    result['bposd_fail_tesseract_success'] = int(gap.sum())
    result['bposd_fail_tesseract_success_rate'] = float(gap.mean())
    result['bposd_success_tesseract_fail'] = int((frame.bposd_success & ~frame.tesseract_success).sum())
    result['both_success'] = int((frame.bposd_success & frame.tesseract_success).sum())
    result['both_fail'] = int((~frame.bposd_success & ~frame.tesseract_success).sum())
    result['mean_bp_iterations'] = float(frame.bp_iterations.mean())
    return result


def alpha_boundaries(frames, alphas):
    rows = []
    for method in ['bp', 'bposd']:
        failures = np.stack([~frame[f'{method}_success'].to_numpy() for frame in frames], axis=1)
        for sample_id, states in enumerate(failures):
            failed = np.flatnonzero(states)
            first = int(failed[0]) if len(failed) else len(alphas)
            rows.append({
                'sample_id': sample_id,
                'decoder': method,
                'alpha_first_fail': float(alphas[first]) if first < len(alphas) else np.nan,
                'alpha_robust': float(alphas[first - 1]) if first > 0 else np.nan,
                'left_censored': first == 0,
                'right_censored': first == len(alphas),
                'failure_transitions': int(np.count_nonzero(states[1:] != states[:-1])),
            })
    return pd.DataFrame(rows)


def prepare_run(paths, config):
    validate_config(config)
    effective = dict(config, backends={
        'bp': BP_TRACE_VERSION,
        'bp_accumulation': 'sequential',
        'bposd': 'cpp_backend.BPOSD',
        'tesseract': 'cpp_backend.Tesseract',
    })
    marker = paths.output('config') / 'effective_config.json'
    if not marker.is_file() or not config['reuse_existing_data']:
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(json.dumps({'config': effective}, indent=2) + '\n')


def run_sweep(paths, config):
    prepare_run(paths, config)
    summaries, boundary_summaries = [], []
    for case in config['cases']:
        case_keys = {'distance': case['distance'], 'weight': case['weight']}
        reference_context = build_context(case['distance'], 0.01, config)
        errors, total_supports = make_supports(reference_context.H.shape[1], case, config)
        support_frame = pd.DataFrame({
            'sample_id': np.arange(len(errors)),
            'error': [bit_text(error) for error in errors],
            'distance': case['distance'],
            'weight': case['weight'],
            'sampling': case['sampling'],
            'num_errors': reference_context.H.shape[1],
        })
        save_csv(support_frame, paths.output('supports_csv', **case_keys))
        metadata = dict(case_keys, sampling=case['sampling'], total_supports=total_supports)
        references = {}
        point_count = sum(
            grid['count'] * (len(alpha_values(config)) if name == 'high_p' else len(config['low_p_alphas']))
            for name, grid in config['p_grids'].items()
        )
        with tqdm(total=point_count, desc=f"d={case['distance']} w={case['weight']}", disable=config['print_level'] == 1) as progress:
            for grid_name, grid in config['p_grids'].items():
                alphas = alpha_values(config) if grid_name == 'high_p' else np.sort(config['low_p_alphas'])
                for p_index, p in enumerate(np.geomspace(grid['start'], grid['stop'], grid['count'])):
                    context = build_context(case['distance'], p, config)
                    assert np.array_equal(context.H, reference_context.H) and np.array_equal(context.L, reference_context.L), 'DEM support ordering changed across p'
                    syndromes = errors @ context.H.T % 2
                    truths = errors @ context.L.T % 2
                    point_frames = []
                    for alpha_index, alpha in enumerate(alphas):
                        parameters = dict(case_keys, grid=grid_name, p_index=p_index, alpha_index=alpha_index)
                        target = paths.output('point_cases_csv', **parameters)
                        if config['reuse_existing_data'] and target.is_file():
                            frame = pd.read_csv(target, dtype={key: str for key in ['error', 'syndrome', 'truth', 'bp_correction', 'bposd_correction']})
                            assert len(frame) == len(errors)
                        else:
                            if float(p) not in references:
                                references[float(p)] = tesseract_outcomes(context, syndromes, truths, p, config)
                            frame = evaluate_point(context, errors, syndromes, truths, references[float(p)], p, alpha, config)
                            frame = frame.assign(**case_keys, p_grid=grid_name)
                            save_csv(frame, target)
                        summaries.append(dict(metadata, p_grid=grid_name, p=float(p), alpha=float(alpha), **summarize_point(frame)))
                        point_frames.append(frame)
                        progress.set_postfix(p=f'{p:g}', alpha=f'{alpha:g}', refresh=False)
                        progress.update()
                    boundaries = alpha_boundaries(point_frames, alphas)
                    boundaries = boundaries.assign(**case_keys, p_grid=grid_name, p=float(p))
                    save_csv(boundaries, paths.output('boundary_cases_csv', **case_keys, grid=grid_name, p_index=p_index))
                    for (method, first, transitions), group in boundaries.groupby(['decoder', 'alpha_first_fail', 'failure_transitions'], dropna=False):
                        boundary_summaries.append(dict(
                            metadata,
                            p_grid=grid_name,
                            p=float(p),
                            decoder=method,
                            alpha_first_fail=first,
                            failure_transitions=int(transitions),
                            samples=len(group),
                            fraction=len(group) / len(errors),
                        ))
                    if config['print_level'] >= 3:
                        print(dict(case_keys, p_grid=grid_name, p=float(p), samples=len(errors)))
    summary = pd.DataFrame(summaries)
    boundary_summary = pd.DataFrame(boundary_summaries)
    save_csv(summary, paths.output('summary_csv'))
    save_csv(boundary_summary, paths.output('boundary_summary_csv'))
    if config['print_level'] >= 2:
        print({'prior_scaling_parameter_points': len(summary), 'summary': str(paths.output('summary_csv'))})
    return summary, boundary_summary


def plot_comparisons(summary, paths):
    paths.plots().mkdir(parents=True, exist_ok=True)
    for (distance, weight, grid), frame in summary.groupby(['distance', 'weight', 'p_grid']):
        fig, axes = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
        for ax, method in zip(axes, ['bp', 'bposd', 'tesseract']):
            matrix = frame.pivot(index='alpha', columns='p', values=f'{method}_failure_rate').sort_index()
            x = np.log10(matrix.columns.to_numpy(dtype=float))
            x_mid = (x[:-1] + x[1:]) / 2
            x_edges = 10 ** np.r_[2 * x[0] - x_mid[0], x_mid, 2 * x[-1] - x_mid[-1]]
            y = matrix.index.to_numpy(dtype=float)
            if len(y) == 1:
                y_edges = np.array([0.5 * y[0], 1.5 * y[0]])
            else:
                y_mid = (y[:-1] + y[1:]) / 2
                y_edges = np.r_[2 * y[0] - y_mid[0], y_mid, 2 * y[-1] - y_mid[-1]]
            mesh = ax.pcolormesh(x_edges, y_edges, matrix.to_numpy(), shading='flat', vmin=0, vmax=1)
            ax.set(
                xscale='log',
                xlabel='Physical error rate p',
                ylabel='Min-sum scaling alpha',
                title=method.upper(),
            )
            fig.colorbar(mesh, ax=ax, label='Conditional failure fraction')
        fig.suptitle(f'd={distance}, weight={weight}, {grid}: fixed DEM error supports')
        fig.savefig(paths.output('phase_plot', distance=int(distance), weight=int(weight), grid=grid), dpi=160)
        plt.show()


def example_case(summary, paths, distance, weight):
    selected = summary[(summary.distance == distance) & (summary.weight == weight) & (summary.bposd_fail_tesseract_success > 0)]
    for _, point in selected.sort_values(['p_grid', 'p', 'alpha']).iterrows():
        grids = summary[(summary.distance == distance) & (summary.weight == weight) & (summary.p_grid == point.p_grid)]
        p_index = int(np.flatnonzero(np.sort(grids.p.unique()) == point.p)[0])
        alpha_index = int(np.flatnonzero(np.sort(grids.alpha.unique()) == point.alpha)[0])
        frame = pd.read_csv(paths.output('point_cases_csv', distance=distance, weight=weight, grid=point.p_grid, p_index=p_index, alpha_index=alpha_index), dtype={'error': str, 'syndrome': str, 'truth': str})
        examples = frame[~frame.bposd_success & frame.tesseract_success]
        if len(examples):
            return examples.iloc[0]
    return None


def plot_example(example, config, paths):
    context = build_context(int(example.distance), float(example.p), config)
    syndrome = np.fromiter(map(int, example.syndrome), dtype=np.uint8)
    trace = context.bp.trace(syndrome, config['max_iter'], scaling=float(example.alpha), record_messages=False)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(trace.llr, alpha=0.5)
    ax.set(xlabel='BP iteration', ylabel='Posterior LLR', title=f'd={int(example.distance)}, weight={int(example.weight)}, p={example.p:g}, alpha={example.alpha:g}')
    fig.savefig(paths.output('example_plot', distance=int(example.distance), weight=int(example.weight)), dpi=160, bbox_inches='tight')
    plt.show()
