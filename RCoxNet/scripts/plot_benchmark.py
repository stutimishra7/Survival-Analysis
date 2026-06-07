"""
benchmark_defaults_figure.py — Publication-ready grouped bar C-index plot comparing
RCoxNet against 5 baseline models using published default hyperparameters.

Input:
  Results/dl_defaults/all_defaults_combined.csv
      Baseline models (Cox-EN, Cox-nnet, SurvivalNet, DeepSurv, DeepHit)
      run with published default HPs — 20 repeats × 4 cancers.

  Results/Pipeline_AllPPI_<CANCER>/final_model_cindex.csv
      RCoxNet AllPPI pipeline final test C-index — 20 repeats × 4 cancers.

Output:
  Results/dl_defaults/benchmark_defaults_figure.pdf
  Results/dl_defaults/benchmark_defaults_figure.png
  Results/dl_defaults/benchmark_defaults_figure.svg
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm
import logging

logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# Font setup (matches benchmark_figure.py)
available_fonts = {f.name for f in fm.fontManager.ttflist}
mpl.rcParams['font.sans-serif'] = ['Arial' if 'Arial' in available_fonts else 'Liberation Sans']
mpl.rcParams['font.family']       = 'sans-serif'
mpl.rcParams['pdf.fonttype']      = 42
mpl.rcParams['ps.fonttype']       = 42
mpl.rcParams['font.size']         = 7
mpl.rcParams['axes.titlesize']    = 8
mpl.rcParams['axes.labelsize']    = 7
mpl.rcParams['xtick.labelsize']   = 7
mpl.rcParams['ytick.labelsize']   = 7
mpl.rcParams['legend.fontsize']   = 6.5
mpl.rcParams['axes.linewidth']    = 0.8
mpl.rcParams['xtick.major.width'] = 0.8
mpl.rcParams['ytick.major.width'] = 0.8
mpl.rcParams['xtick.major.size']  = 3
mpl.rcParams['ytick.major.size']  = 3

# Paths
ROOT      = os.path.dirname(os.path.abspath(__file__))
DEFAULTS  = os.path.join(ROOT, 'Results', 'dl_defaults')
ALLPPI    = os.path.join(ROOT, 'Results')
OUT_DIR   = DEFAULTS

CANCERS   = ['BRCA', 'GBM', 'LUNG', 'OV']

# map CSV model names to display labels
BASELINE_MAP = {
    'Cox-EN_defaults':      'Cox-EN',
    'Cox-nnet_defaults':    'Cox-nnet',
    'SurvivalNet_defaults': 'SurvivalNet',
    'DeepSurv_defaults':    'DeepSurv',
    'DeepHit_defaults':     'DeepHit',
}

# Display order: RCoxNet first (highlighted), then baselines
MODEL_ORDER = ['RCoxNet', 'Cox-EN', 'Cox-nnet', 'SurvivalNet', 'DeepSurv', 'DeepHit']

# Palette — matte, well-separated, user-anchored:
#   RCoxNet:  matte green  #48A14D
#   DeepHit:  matte burnt orange  #C04000
MODEL_CFG = {
    'RCoxNet':     dict(color='#48A14D', label='RCoxNet'),      # matte forest green
    'Cox-EN':      dict(color='#5B8DB8', label='Cox-EN'),        # matte steel blue
    'Cox-nnet':    dict(color='#9B6BB5', label='Cox-nnet'),      # matte violet
    'SurvivalNet': dict(color='#D4A017', label='SurvivalNet'),  # matte gold
    'DeepSurv':    dict(color='#4E9E8E', label='DeepSurv'),     # matte teal
    'DeepHit':     dict(color='#C04000', label='DeepHit'),       # matte burnt orange
}

# Load baseline results
def load_baselines():
    path = os.path.join(DEFAULTS, 'all_defaults_combined.csv')
    if not os.path.exists(path):
        raise FileNotFoundError(f'Baseline results not found: {path}')
    df = pd.read_csv(path)
    df['model'] = df['model'].map(BASELINE_MAP)
    df = df.dropna(subset=['model', 'test_c_index'])
    return df[['cancer', 'model', 'repeat', 'test_c_index']]

# Load RCoxNet AllPPI results
def load_rcoxnet():
    rows = []
    for cancer in CANCERS:
        path = os.path.join(ALLPPI, f'Pipeline_AllPPI_{cancer}', 'final_model_cindex.csv')
        if not os.path.exists(path):
            print(f'[warn] RCoxNet result not found: {path}')
            continue
        df = pd.read_csv(path)
        for rep, ci in enumerate(df['test_ci'].values, start=1):
            rows.append({'cancer': cancer, 'model': 'RCoxNet',
                         'repeat': rep, 'test_c_index': float(ci)})
    if not rows:
        raise FileNotFoundError('No RCoxNet Pipeline_AllPPI results found.')
    return pd.DataFrame(rows)

# Combine & summarise
def build_summary(data):
    summary = (data.groupby(['cancer', 'model'])['test_c_index']
               .agg(mean='mean', std='std', n='count')
               .reset_index())
    return summary

# Plot
def make_figure(summary, out_dir):
    present = [m for m in MODEL_ORDER if m in summary['model'].unique()]

    n_cancers = len(CANCERS)
    n_models  = len(present)
    group_w   = 0.80
    bar_w     = group_w / n_models
    offsets   = np.linspace(-(group_w - bar_w) / 2,
                             (group_w - bar_w) / 2,
                             n_models)

    fig, ax = plt.subplots(figsize=(4.5, 3.0), dpi=300)

    for mi, model_key in enumerate(present):
        color  = MODEL_CFG[model_key]['color']
        is_rcox = model_key == 'RCoxNet'

        for ci, cancer in enumerate(CANCERS):
            row = summary[(summary['cancer'] == cancer) &
                          (summary['model']  == model_key)]
            if row.empty:
                continue
            mean = row['mean'].values[0]
            std  = row['std'].values[0]
            x    = ci + offsets[mi]

            ax.bar(x, mean,
                   width      = bar_w * 0.88,
                   color      = color,
                   edgecolor  = 'none',
                   linewidth  = 0,
                   zorder     = 3        if is_rcox else 2)

            ax.errorbar(x, mean, yerr=std,
                        fmt='none',
                        color='#333333',
                        capsize=2.0, capthick=0.9,
                        linewidth=0.9,
                        zorder=4)

    # Y axis
    all_vals = summary['mean'] + summary['std']
    ymax = max(0.90, round(float(all_vals.max()) + 0.06, 1))
    ax.set_ylim(0.40, ymax)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.10))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.tick_params(axis='y', which='minor', length=2, width=0.6)

    # X axis
    ax.set_xticks(range(n_cancers))
    ax.set_xticklabels(CANCERS, fontweight='bold', fontsize=7)
    ax.set_xlim(-0.55, n_cancers - 0.45)

    # Labels
    ax.set_xlabel('Cancer type', labelpad=4, fontsize=7)
    ax.set_ylabel('C-index (mean ± SD, 20 repeats)', labelpad=4, fontsize=7)

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(axis='both', which='major',
                   width=0.8, length=3, pad=2, colors='black')
    ax.grid(False)

    # Legend
    handles = [
        mpatches.Patch(
            facecolor  = MODEL_CFG[m]['color'],
            edgecolor  = 'none',
            linewidth  = 0,
            label      = MODEL_CFG[m]['label'])
        for m in present
    ]
    ax.legend(handles=handles,
              loc='upper center', bbox_to_anchor=(0.5, 1.02),
              ncol=n_models, frameon=False,
              handlelength=1.0, handleheight=0.75,
              columnspacing=0.6, handletextpad=0.35,
              borderpad=0, fontsize=6.5)

    plt.tight_layout(pad=0.6)

    os.makedirs(out_dir, exist_ok=True)
    for ext in ('pdf', 'png', 'svg'):
        out = os.path.join(out_dir, f'benchmark_defaults_figure.{ext}')
        kw  = dict(bbox_inches='tight', transparent=(ext != 'png'))
        if ext == 'png':
            kw['dpi'] = 300
        else:
            kw['format'] = ext
        fig.savefig(out, **kw)
        print(f'Saved → {out}')

    plt.close(fig)
    return summary

# Summary table
def print_summary(summary):
    print('\nC-index summary (mean ± std, 20 repeats)')
    present = [m for m in MODEL_ORDER if m in summary['model'].unique()]
    pivot   = summary.copy()
    pivot['value'] = pivot.apply(
        lambda r: f"{r['mean']:.4f} ± {r['std']:.4f}", axis=1)
    tbl = (pivot.pivot(index='model', columns='cancer', values='value')
                .reindex(present)[CANCERS])
    print(tbl.to_string())

    print('\nRCoxNet improvement over best baseline')
    baselines = [m for m in present if m != 'RCoxNet']
    for cancer in CANCERS:
        rcox = summary[(summary['cancer']==cancer) &
                       (summary['model']=='RCoxNet')]['mean'].values
        if not len(rcox):
            continue
        best_bl = summary[(summary['cancer']==cancer) &
                          (summary['model'].isin(baselines))]['mean'].max()
        best_model = summary[(summary['cancer']==cancer) &
                             (summary['model'].isin(baselines))].loc[
                        summary[(summary['cancer']==cancer) &
                                (summary['model'].isin(baselines))]['mean'].idxmax(),
                        'model']
        print(f'  {cancer}: RCoxNet={rcox[0]:.4f}  '
              f'best_baseline={best_bl:.4f} ({best_model})  '
              f'Δ=+{rcox[0]-best_bl:.4f}')

# Main
if __name__ == '__main__':
    print('Loading baseline results ...')
    baselines = load_baselines()
    print(f'  {len(baselines)} rows  |  '
          f'models: {sorted(baselines["model"].unique())}')

    print('Loading RCoxNet AllPPI results ...')
    rcoxnet = load_rcoxnet()
    print(f'  {len(rcoxnet)} rows  |  cancers: {sorted(rcoxnet["cancer"].unique())}')

    data    = pd.concat([baselines, rcoxnet], ignore_index=True)
    summary = build_summary(data)

    print('Generating figure ...')
    make_figure(summary, OUT_DIR)
    print_summary(summary)
