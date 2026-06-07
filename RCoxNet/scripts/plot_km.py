"""
Plot Kaplan-Meier survival curves from pre-computed risk scores.

Risk scores must be generated first by running pipeline_rcoxnet.ipynb
(Steps 1-8), which saves:
    results/km_plots/km_scores_{CANCER}_notebook.csv

Usage
-----
    python scripts/plot_km.py --cancer BRCA
    python scripts/plot_km.py --cancer ALL
    python scripts/plot_km.py --cancer BRCA --out results/km_plots/km_BRCA.pdf
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

CANCERS = ['BRCA', 'GBM', 'LUNG', 'OV']

HIGH_COLOR = '#C0392B'
LOW_COLOR  = '#1565C0'


def plot_km(cancer: str, scores_csv: str, ax=None, save_path: str = None):
    df = pd.read_csv(scores_csv)

    scores       = df['avg_risk_score'].values
    s_min, s_max = scores.min(), scores.max()
    scaled       = 2 * (scores - s_min) / (s_max - s_min) - 1

    thresh = float(np.median(scaled))
    high   = scaled >= thresh
    low    = ~high

    t_hi = df.loc[high, 'OS_MONTHS'].values
    t_lo = df.loc[low,  'OS_MONTHS'].values
    e_hi = df.loc[high, 'OS_STATUS'].values
    e_lo = df.loc[low,  'OS_STATUS'].values

    p        = logrank_test(t_hi, t_lo,
                            event_observed_A=e_hi,
                            event_observed_B=e_lo).p_value
    pval_str = f'{p:.5f}' if p >= 1e-5 else f'{p:.2e}'

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(4.5, 4))

    kmf = KaplanMeierFitter()
    kmf.fit(t_hi, e_hi)
    kmf.plot_survival_function(ax=ax, ci_show=False,
                               color=HIGH_COLOR, linewidth=1.4,
                               show_censors=True,
                               censor_styles={'marker': '+', 'ms': 4, 'mew': 0.8})
    kmf.fit(t_lo, e_lo)
    kmf.plot_survival_function(ax=ax, ci_show=False,
                               color=LOW_COLOR, linewidth=1.4,
                               show_censors=True,
                               censor_styles={'marker': '+', 'ms': 4, 'mew': 0.8})

    ax.get_legend().remove()
    handles = [
        Line2D([0], [0], color=HIGH_COLOR, linewidth=1.4,
               label=f'PI≥{thresh:.3f}  (n={high.sum()})'),
        Line2D([0], [0], color=LOW_COLOR,  linewidth=1.4,
               label=f'PI<{thresh:.3f}   (n={low.sum()})'),
    ]
    ax.legend(handles=handles, loc='upper right', frameon=False, fontsize=8)
    ax.text(0.97, 0.55, f'p-value: {pval_str}',
            transform=ax.transAxes, fontsize=8,
            va='top', ha='right', color='#222222')

    ax.set_title(cancer, fontsize=11, fontweight='bold')
    ax.set_xlabel('Overall Survival in Months', fontsize=9)
    ax.set_ylabel('Survival Probability', fontsize=9)
    ax.set_xlim(left=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    print(f'{cancer}  log-rank p = {p:.4e}')

    if standalone:
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, bbox_inches='tight', dpi=300)
            print(f'Saved → {save_path}')
        else:
            plt.show()
        plt.close(fig)

    return p


def plot_all(root: str, save_path: str = None):
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    for ax, cancer in zip(axes, CANCERS):
        csv = os.path.join(root, 'results', 'km_plots',
                           f'km_scores_{cancer}_notebook.csv')
        if not os.path.exists(csv):
            print(f'Missing: {csv}  (run Steps 1-8 in pipeline_rcoxnet.ipynb first)')
            continue
        plot_km(cancer, csv, ax=ax)

    axes[0].set_ylabel('Survival Probability', fontsize=9)
    for ax in axes[1:]:
        ax.set_ylabel('')

    plt.tight_layout()
    out = save_path or os.path.join(root, 'results', 'km_plots', 'Figure2_KM.pdf')
    fig.savefig(out, bbox_inches='tight', dpi=300)
    print(f'Saved → {out}')
    plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Plot KM curves from pre-computed risk scores')
    parser.add_argument('--cancer', default='ALL', choices=CANCERS + ['ALL'])
    parser.add_argument('--out', default=None, help='output file path (.pdf or .png)')
    parser.add_argument('--root', default=None, help='repository root directory')
    args = parser.parse_args()

    root = args.root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if args.cancer == 'ALL':
        plot_all(root, save_path=args.out)
    else:
        csv = os.path.join(root, 'results', 'km_plots',
                           f'km_scores_{args.cancer}_notebook.csv')
        if not os.path.exists(csv):
            print(f'Missing: {csv}\nRun Steps 1-8 in pipeline_rcoxnet.ipynb first.')
            return
        out = args.out or os.path.join(
            root, 'results', 'km_plots', f'km_{args.cancer}.pdf')
        plot_km(args.cancer, csv, save_path=out)


if __name__ == '__main__':
    main()
