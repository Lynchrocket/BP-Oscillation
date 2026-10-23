# failure_regions, force_osd_failures
"""Reusable experiment functions with explicit arguments and dependencies."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def stats_to_df(stats_list):
    rows = []
    for s in stats_list:
        shots = int(s.shots)
        errors = int(s.errors)
        rows.append({'code': s.json_metadata['code'], 'decoder': s.json_metadata['decoder'], 'p': float(s.json_metadata['p']), 'shots': shots, 'errors': errors, 'logical_error_rate': errors / shots if shots > 0 else np.nan, 'seconds': float(s.seconds), 'ms_per_shot': float(s.seconds) / shots * 1000.0 if shots > 0 else np.nan})
    return pd.DataFrame(rows)

def plot_gap_interval(summary: pd.DataFrame, x_scale=('log', 'log'), y_scale=('log', 'log')):
    plot_codes = sorted(summary['code'].unique())
    fig, axes = plt.subplots(len(plot_codes), 2, figsize=(12, 4 * len(plot_codes)), squeeze=False)
    for r, code_name in enumerate(plot_codes):
        d = summary[summary['code'] == code_name].sort_values('p').copy()
        d['gap_bposd_minus_tess'] = d['BPOSD'] - d['TESSERACT']
        ax = axes[r, 0]
        ax.plot(d['p'], d['BPOSD'], marker='o', label='BPOSD')
        ax.plot(d['p'], d['TESSERACT'], marker='s', label='Tesseract')
        ax.set_xscale(x_scale[0])
        ax.set_yscale(y_scale[0])
        ax.set_title(f'{code_name}: LER vs p')
        ax.set_xlabel('physical error rate p')
        ax.set_ylabel('logical error rate')
        ax.grid(True, which='both', ls='--', alpha=0.4)
        ax.legend()
        axg = axes[r, 1]
        axg.plot(d['p'], d['gap_bposd_minus_tess'], marker='d', color='tab:red')
        axg.axhline(0.0, color='black', lw=1)
        axg.set_xscale(x_scale[1])
        axg.set_yscale(y_scale[1])
        axg.set_title(f'{code_name}: gap = BPOSD - Tesseract')
        axg.set_xlabel('physical error rate p')
        axg.set_ylabel('gap')
        axg.grid(True, which='both', ls='--', alpha=0.4)
        pos = d[d['gap_bposd_minus_tess'] > 0]
        if len(pos) > 0:
            axg.scatter(pos['p'], pos['gap_bposd_minus_tess'], color='tab:green', zorder=3, label='gap>0')
            axg.legend()
    plt.tight_layout()
    plt.show()

def group_func(stat):
    """
    One curve for each (code, decoder) pair.
    """
    return (stat.json_metadata['code'], stat.json_metadata['decoder'])

def plot_args_func(index, curve_id, curve_stats, *, code_style, decoder_style):
    code, decoder = curve_id
    c = code_style[code]
    d = decoder_style[decoder]
    return {'color': c['color'], 'marker': d['marker'], 'linestyle': d['linestyle'], 'linewidth': 1.4, 'markersize': 4.5, 'label': f"{d['label']}, {c['label']}"}

def extract_gap_intervals(df_code: pd.DataFrame, gap_col='gap_bposd_minus_tess', min_gap=0.0):
    d = df_code.sort_values('p').reset_index(drop=True).copy()
    d['is_gap'] = d[gap_col] > min_gap
    intervals = []
    start_idx = None
    for i, v in enumerate(d['is_gap'].tolist()):
        if v and start_idx is None:
            start_idx = i
        if not v and start_idx is not None:
            seg = d.iloc[start_idx:i]
            intervals.append({'p_start': float(seg['p'].iloc[0]), 'p_end': float(seg['p'].iloc[-1]), 'n_points': int(len(seg)), 'max_gap': float(seg[gap_col].max()), 'mean_gap': float(seg[gap_col].mean())})
            start_idx = None
    if start_idx is not None:
        seg = d.iloc[start_idx:]
        intervals.append({'p_start': float(seg['p'].iloc[0]), 'p_end': float(seg['p'].iloc[-1]), 'n_points': int(len(seg)), 'max_gap': float(seg[gap_col].max()), 'mean_gap': float(seg[gap_col].mean())})
    return pd.DataFrame(intervals)

