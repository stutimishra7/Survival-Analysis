"""
Ablation study + KM survival plots for all 4 cancers (BRCA, GBM, LUNG, OV).

For each cancer:
  1. Loads data, builds PPI, computes RWR with best (alpha, K)
  2. Trains CoxMLP / ClinicalOnly / RWROnly (20 repeats) — exact match Cell 32
  3. Plots C-index bar + box + per-repeat comparison
  4. Runs 20-repeat inference-time ablation (zero-out RWR / clinical features)
  5. Plots ablation bar charts (absolute C-index + delta)
  6. Runs pooled KM survival analysis (High Risk vs Low Risk)
  7. Plots publication-ready KM curves

Output:
  Results/<cancer>_cindex_bar.png      — C-index comparison bar chart
  Results/<cancer>_cindex_box.png      — C-index boxplot
  Results/ablation_<cancer>.png        — ablation bar chart
  Results/km_pooled_<cancer>.png       — KM survival curve
  Results/cindex_summary.csv           — CoxMLP/ClinOnly/RWROnly C-index
  Results/ablation_summary.csv         — ablation numbers for all cancers
  Results/km_summary.csv               — KM p-values and median survival
"""

import copy, os, sys, warnings
import numpy as np
import pandas as pd
import torch
import torch.optim as optim
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

sys.path.insert(0, '/home/stutik/RCoxNet_Pipeline_Final')

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from rcoxnet.config import (
    CANCER_CONFIGS, CLIN_FEATURES,
    N_REPEATS, TEST_FRAC, VAL_FRAC,
    EPOCHS, CHECKPOINT, GRAD_CLIP,
)
from rcoxnet.data_loader   import load_raw_data
from rcoxnet.mutations     import parse_mutations
from rcoxnet.clinical      import extract_clinical
from rcoxnet.network       import build_transition_matrix
from rcoxnet.rwr           import compute_rwr_scores
from rcoxnet.features      import make_feature_df, select_top_k
from rcoxnet.preprocessing import (
    censoring_time_balanced_split, normalise_genomic,
    normalise_clinical, prep as _prep,
)
from rcoxnet.model import CoxPhRWRNet, cox_loss, c_index
from rcoxnet.train import train_and_eval

ROOT    = '/home/stutik/RCoxNet_Pipeline_Final'
RESULTS = os.path.join(ROOT, 'Results')
device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
dtype   = torch.float32

CANCERS = ['BRCA', 'GBM', 'LUNG', 'OV']
EPOCHS  = 2000

MASK_CONFIGS = {
    'All features (full)': {},
    'Zero RWR':            {'zero_genomic': True},
    'Zero Clinical':       {'zero_clinical': True},
    'Zero AGE':            {'zero_features': ['AGE']},
    'Zero MSI':            {'zero_features': ['msi_score']},
    'Zero TMB':            {'zero_features': ['tmb_score']},
    'Zero AGE+MSI':        {'zero_features': ['AGE', 'msi_score']},
    'Zero AGE+TMB':        {'zero_features': ['AGE', 'tmb_score']},
    'Zero MSI+TMB':        {'zero_features': ['msi_score', 'tmb_score']},
}

# helpers

def prep(Xn, Cn, yt, ye):
    return _prep(device, dtype, Xn, Cn, yt, ye)

def train_one(mk, x_tr, cl_tr, yt_tr, ye_tr,
              x_va, cl_va, yt_va, ye_va,
              x_te, cl_te, yt_te, ye_te,
              lr, l2, seed, return_net=False):
    # exact mirror of notebook train_model (Cell 33)
    torch.manual_seed(seed); np.random.seed(seed)
    net = CoxPhRWRNet(**mk).to(device)
    opt = optim.Adam(net.parameters(), lr=lr, weight_decay=l2)
    sch = optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS, eta_min=1e-6)
    best_val, best_st = 0.0, copy.deepcopy(net.state_dict())
    for ep in range(EPOCHS):
        net.train(); opt.zero_grad()
        loss = cox_loss(net(x_tr, cl_tr), yt_tr, ye_tr)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(net.parameters(), GRAD_CLIP)
        opt.step(); sch.step()
        if ep % CHECKPOINT == 0:
            net.eval()
            with torch.no_grad():
                vl = cox_loss(net(x_va, cl_va), yt_va, ye_va).item()
                vc = c_index(net(x_va, cl_va), yt_va, ye_va).item()
                tl = cox_loss(net(x_tr, cl_tr), yt_tr, ye_tr).item()
            if vc > best_val:
                best_val, best_st = vc, copy.deepcopy(net.state_dict())
    net.load_state_dict(best_st); net.eval()
    with torch.no_grad():
        tc  = c_index(net(x_te, cl_te), yt_te, ye_te).item()
        raw = net(x_te, cl_te).cpu().numpy().flatten()
    if return_net:
        return tc, best_val, raw, net
    return tc, best_val, raw

def eval_with_mask(net, x_te, cl_te, yt_te, ye_te,
                   zero_genomic=False, zero_clinical=False, zero_features=None):
    clin_list = list(CLIN_FEATURES)
    x_in  = torch.zeros_like(x_te)  if zero_genomic  else x_te.clone()
    cl_in = torch.zeros_like(cl_te) if zero_clinical else cl_te.clone()
    if not zero_clinical and zero_features:
        for feat in zero_features:
            if feat in clin_list:
                cl_in[:, clin_list.index(feat)] = 0.0
    net.eval()
    with torch.no_grad():
        return c_index(net(x_in, cl_in), yt_te, ye_te).item()

import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
import matplotlib.ticker as ticker
from matplotlib.colors import to_rgba
from scipy.stats import wilcoxon

# publication rcParams (matches generate_figure4.py exactly)
_available_fonts = {f.name for f in fm.fontManager.ttflist}
_FONT = 'Arial' if 'Arial' in _available_fonts else 'Liberation Sans'
matplotlib.rcParams.update({
    'pdf.fonttype'      : 42,  'ps.fonttype'       : 42,
    'font.family'       : 'sans-serif',
    'font.sans-serif'   : [_FONT],
    'font.size'         : 7,   'axes.titlesize'    : 8,
    'axes.labelsize'    : 7,   'xtick.labelsize'   : 7,
    'ytick.labelsize'   : 7,   'legend.fontsize'   : 6.5,
    'axes.linewidth'    : 0.8, 'xtick.major.width' : 0.8,
    'ytick.major.width' : 0.8, 'xtick.major.size'  : 3,
    'ytick.major.size'  : 3,   'figure.facecolor'  : 'white',
    'axes.facecolor'    : 'white',
})

CINDEX_DIR = os.path.join(RESULTS, 'cindex_plots')
os.makedirs(CINDEX_DIR, exist_ok=True)

MODELS  = ['ClinicalOnly', 'RWROnly', 'CoxMLP']   # CoxMLP rightmost (proposed)
LABELS  = ['Clinical-only', 'RWR-only', 'RCoxNet']
BASELINES_CI = ['ClinicalOnly', 'RWROnly']

MODEL_CFG_CI = {
    'ClinicalOnly': dict(color='#7A9CB8', label='Clinical-only'),
    'RWROnly':      dict(color='#C9A84C', label='RWR-only'),
    'CoxMLP':       dict(color='#6A9E6E', label='RCoxNet'),
}
PANEL_LETTERS = {'BRCA': 'a', 'GBM': 'b', 'LUNG': 'c', 'OV': 'd'}
CANCER_NAMES  = {
    'BRCA': 'Breast (BRCA)', 'GBM': 'Glioblastoma (GBM)',
    'LUNG': 'Lung (LUNG)',   'OV':  'Ovarian (OV)',
}

def _wilcoxon_p(a, b):
    try:
        _, p = wilcoxon(np.array(a), np.array(b), alternative='greater')
        return p
    except Exception:
        return 1.0

def _sig_label(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

def _add_brackets(ax, all_pairs, top_y, dy=0.028, tip=0.010):
    for x1, x2, lbl in all_pairs:
        is_sig = lbl != 'ns'
        color  = '#1a1a1a' if is_sig else '#999999'
        lw     = 0.9       if is_sig else 0.7
        ls     = '-'       if is_sig else '--'
        h = top_y + tip
        ax.plot([x1, x1, x2, x2],
                [h, h + dy*0.45, h + dy*0.45, h],
                lw=lw, color=color, linestyle=ls,
                clip_on=False, zorder=6, solid_capstyle='round')
        ax.text((x1+x2)/2, h + dy*0.45 + 0.005, lbl,
                ha='center', va='bottom',
                fontsize=5.5 if lbl == 'ns' else 6.5,
                fontweight='bold' if is_sig else 'normal',
                fontstyle='italic' if lbl == 'ns' else 'normal',
                color=color, clip_on=False)
        top_y += dy

def plot_cindex(cancer, all_rows):
    ci    = {m: [r['test_c'] for r in all_rows[m]] for m in MODELS}
    means = {m: np.mean(ci[m]) for m in MODELS}
    stds  = {m: np.std(ci[m])  for m in MODELS}
    return means, stds, ci

def make_cindex_figure(all_ci_data, show_dots=True):
    """Combined 1×4 figure matching generate_figure4.py exactly."""
    all_vals = np.concatenate([v for ci in all_ci_data.values()
                                for v in ci.values()])
    n_brackets = len(BASELINES_CI)
    ymin = max(0.35, all_vals.min() - 0.06)
    ymax = min(1.02, all_vals.max() + 0.06 + n_brackets * 0.038)

    fig, axes = plt.subplots(1, 4, figsize=(7.09, 4.4), dpi=300, sharey=True)
    fig.subplots_adjust(wspace=0.06, left=0.10, right=0.985,
                        top=0.76, bottom=0.22)

    for col_i, (ax, cancer) in enumerate(zip(axes, CANCERS)):
        ci     = all_ci_data[cancer]
        vals   = [ci[m] for m in MODELS]
        colors = [MODEL_CFG_CI[m]['color'] for m in MODELS]

        # reference line
        ax.axhline(0.5, color='#cccccc', linewidth=0.6,
                   linestyle=':', zorder=1)

        # boxplot
        bp = ax.boxplot(
            vals,
            patch_artist  = True,
            widths        = 0.52,
            medianprops   = dict(color='#111111', linewidth=1.4, zorder=5),
            whiskerprops  = dict(color='#606060', linewidth=0.75),
            capprops      = dict(color='#606060', linewidth=0.75),
            flierprops    = dict(marker='o', markersize=1.8,
                                 markerfacecolor='#aaaaaa',
                                 markeredgewidth=0, alpha=0.6),
            showfliers    = True,
            zorder        = 2,
        )
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.85)
            patch.set_linewidth(0.5)
            patch.set_edgecolor('#555555')

        # CoxMLP (RCoxNet) — darker green border
        rcox_pos = MODELS.index('CoxMLP')
        bp['boxes'][rcox_pos].set_edgecolor('#2e5e32')
        bp['boxes'][rcox_pos].set_linewidth(0.9)

        # jitter dots
        if show_dots:
            np.random.seed(42)
            for i, (d, color) in enumerate(zip(vals, colors)):
                jx = np.random.uniform(-0.15, 0.15, size=len(d))
                ax.scatter([i+1+j for j in jx], d,
                           color=color, s=4.5, alpha=0.40,
                           edgecolors='none', zorder=3)

        # significance brackets
        rcox_vals = ci['CoxMLP']
        top_data  = max(np.max(v) for v in vals)
        all_pairs = []
        for bl in BASELINES_CI:
            bi   = MODELS.index(bl)
            p    = _wilcoxon_p(rcox_vals, ci[bl])
            lbl  = _sig_label(p)
            span = abs(rcox_pos - bi)
            all_pairs.append((bi+1, rcox_pos+1, lbl, span))
        all_pairs.sort(key=lambda x: x[3])
        all_pairs = [(x1, x2, lbl) for x1, x2, lbl, _ in all_pairs]
        _add_brackets(ax, all_pairs, top_y=top_data, dy=0.028, tip=0.010)

        # axes styling
        ax.set_title(CANCER_NAMES[cancer], fontsize=7.5,
                     fontweight='bold', pad=5)
        ax.set_xticks(range(1, len(MODELS)+1))
        ax.set_xticklabels([MODEL_CFG_CI[m]['label'] for m in MODELS],
                           rotation=40, ha='right',
                           fontsize=6.5, rotation_mode='anchor')
        ax.set_ylim(ymin, ymax)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.10))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
        ax.tick_params(axis='y', which='minor', length=2, width=0.5)
        ax.tick_params(axis='both', which='major',
                       width=0.8, length=3, pad=2, colors='black')

        if col_i == 0:
            ax.set_ylabel('C-index', fontsize=7, labelpad=4)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.8)
        ax.spines['bottom'].set_linewidth(0.8)
        if col_i > 0:
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='y', which='both', left=False)
        ax.grid(False)

    # legend
    handles = [mpatches.Patch(facecolor=MODEL_CFG_CI[m]['color'],
                               edgecolor='none', linewidth=0,
                               label=MODEL_CFG_CI[m]['label'])
               for m in MODELS]
    fig.legend(handles=handles, loc='lower center',
               bbox_to_anchor=(0.54, -0.02), ncol=len(MODELS),
               frameon=False, handlelength=0.9, handleheight=0.7,
               columnspacing=0.8, handletextpad=0.4,
               borderpad=0, fontsize=6.5)

    suffix = 'with_dots' if show_dots else 'no_dots'
    for ext in ('pdf', 'png', 'svg'):
        out = os.path.join(CINDEX_DIR, f'cindex_all_cancers_1x4_{suffix}.{ext}')
        kw  = dict(bbox_inches='tight', transparent=(ext != 'png'))
        if ext == 'png': kw['dpi'] = 300
        else: kw['format'] = ext
        fig.savefig(out, **kw)
        print(f'  Saved → {out}')
    plt.close(fig)

def km_stats(t_high, t_low, e_high, e_low):
    lr   = logrank_test(t_high, t_low,
                        event_observed_A=e_high, event_observed_B=e_low)
    kmh  = KaplanMeierFitter().fit(t_high, e_high)
    kml  = KaplanMeierFitter().fit(t_low,  e_low)
    return lr.p_value, kmh.median_survival_time_, kml.median_survival_time_, kmh, kml

KM_DIR = os.path.join(RESULTS, 'km_plots')
os.makedirs(KM_DIR, exist_ok=True)

KM_COLORS = {'High Risk': '#C0392B', 'Low Risk': '#2471A3'}   # red / blue
# {cancer: {model_key: (kmh, kml, pval, n_high, n_low)}}
all_km_data = {}

KM_MODEL_CONFIGS = {
    'CoxMLP':       dict(x_zero=False, cl_zero=False, label='RWR+Clinical'),
    'ClinicalOnly': dict(x_zero=True,  cl_zero=False, label='Clinical-only'),
    'RWROnly':      dict(x_zero=False, cl_zero=True,  label='RWR-only'),
}

def plot_km(ax, kmh, kml, pval, cancer, n_high, n_low,
            letter=None, show_ylabel=True):
    ps = f'p = {pval:.2e}'

    kmh.plot_survival_function(
        ax=ax, ci_show=True, ci_alpha=0.15,
        color=KM_COLORS['High Risk'], linewidth=0.8)
    kml.plot_survival_function(
        ax=ax, ci_show=True, ci_alpha=0.15,
        color=KM_COLORS['Low Risk'], linewidth=0.8)
    legend = ax.get_legend()
    if legend:
        legend.remove()

    ax.set_title(cancer, fontsize=7, pad=3)
    ax.set_xlabel('OS Months', labelpad=2)
    if show_ylabel:
        ax.set_ylabel('Survival Probability', labelpad=3)
    else:
        ax.set_ylabel('')

    ax.set_ylim(0, 1.05)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # p-value top-right, n counts bottom-right — no legend overlap
    ax.text(0.97, 0.97, ps,
            transform=ax.transAxes, fontsize=5.5, va='top', ha='right',
            color='#222222')
    ax.text(0.97, 0.03,
            f'High n={n_high}  Low n={n_low}',
            transform=ax.transAxes, fontsize=5, va='bottom', ha='right',
            color='#555555')

    if letter:
        ax.text(0.04, 0.97, letter, transform=ax.transAxes,
                fontsize=8, fontweight='bold', va='top', ha='left',
                color='#222222')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.tick_params(axis='both', width=0.4, length=2.0, pad=1.5)
    ax.yaxis.grid(False); ax.xaxis.grid(False)
    ax.set_facecolor('white')

def make_km_figure(km_data, model_key, suffix=''):
    """Combined 1×4 KM figure for one model type."""
    mm    = 1 / 25.4
    FIG_W = 174 * mm
    FIG_H =  72 * mm   # was 62 — extra 10mm prevents legend/x-axis overlap

    fig, axes = plt.subplots(1, 4, figsize=(FIG_W, FIG_H), dpi=300, sharey=True)
    # bottom=0.22 gives ~16mm below axes for x-ticks + label + legend
    fig.subplots_adjust(left=0.065, right=0.995, top=0.88, bottom=0.22, wspace=0.08)

    for col, cancer in enumerate(CANCERS):
        kmh, kml, pval, n_high, n_low = km_data[cancer][model_key]
        plot_km(axes[col], kmh, kml, pval, cancer, n_high, n_low,
                letter=PANEL_LETTERS[cancer], show_ylabel=(col == 0))

    tag = 'km_all_cancers_1x4' if suffix == 'coxmlp' else f'km_{suffix}_1x4'

    # shared legend inside figure boundary — avoids bbox_inches='tight' collision
    handles = [
        mpatches.Patch(facecolor=to_rgba(KM_COLORS['High Risk'], 0.45),
                       edgecolor=KM_COLORS['High Risk'], linewidth=0.6,
                       label='High Risk'),
        mpatches.Patch(facecolor=to_rgba(KM_COLORS['Low Risk'],  0.45),
                       edgecolor=KM_COLORS['Low Risk'],  linewidth=0.6,
                       label='Low Risk'),
    ]
    fig.legend(handles=handles, loc='lower center',
               bbox_to_anchor=(0.525, 0.01), ncol=2, frameon=False,
               handlelength=0.9, handleheight=0.7,
               handletextpad=0.35, columnspacing=0.8, fontsize=6.5)
    out_pdf = os.path.join(KM_DIR, f'{tag}.pdf')
    fig.savefig(out_pdf, format='pdf', bbox_inches='tight', dpi=300)
    fig.savefig(out_pdf.replace('.pdf', '.png'), format='png',
                bbox_inches='tight', dpi=300)
    plt.close()
    print(f'  KM 1×4 saved → {out_pdf}')

def plot_ablation(cancer, abl_results):
    labels     = list(abl_results.keys())
    means_abl  = [np.mean(abl_results[k]) for k in labels]
    stds_abl   = [np.std(abl_results[k])  for k in labels]
    ref        = means_abl[0]
    drops      = [m - ref for m in means_abl]
    colors_abl = ['#2196F3' if k == 'All features (full)'
                  else ('#E53935' if d < -0.01 else '#43A047')
                  for k, d in zip(labels, drops)]

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    ax = axes[0]
    y  = np.arange(len(labels))
    ax.barh(y, means_abl, xerr=stds_abl, color=colors_abl,
            alpha=0.82, capsize=4, edgecolor='white')
    ax.axvline(ref, color='navy', linestyle='--', linewidth=1.5,
               label=f'Full model ({ref:.4f})')
    ax.axvline(0.5, color='gray', linestyle=':', linewidth=1, alpha=0.7,
               label='Random (0.5)')
    for i, (m, s) in enumerate(zip(means_abl, stds_abl)):
        ax.text(m + s + 0.002, i, f'{m:.4f}', va='center', fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Test C-index', fontsize=11)
    ax.set_title(f'{cancer} — Inference-Time Ablation\nAbsolute C-index', fontsize=11)
    ax.legend(fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    ax.barh(y, drops, xerr=stds_abl, color=colors_abl,
            alpha=0.82, capsize=4, edgecolor='white')
    ax.axvline(0, color='black', linewidth=1.2)
    for i, (d, s) in enumerate(zip(drops, stds_abl)):
        ax.text(d + (s + 0.002 if d >= 0 else -s - 0.018), i,
                f'{d:+.4f}', va='center', fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Δ C-index vs full model', fontsize=11)
    ax.set_title(f'{cancer} — Inference-Time Ablation\nDrop from full model', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    fpath = os.path.join(RESULTS, f'ablation_{cancer.lower()}.png')
    plt.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Ablation plot saved → {fpath}')
    return means_abl, stds_abl, drops

# main loop

ablation_summary  = []
km_summary        = []
cindex_summary    = []
clinical_km_summary = []
all_ci_data       = {}   # {cancer: {model: [c-index values]}}
all_clinical_km   = {}   # {cancer: {feature: (kmh, kml, pval, median_val)}}

CLIN_FEAT_LABELS = {
    'AGE':       'Age',
    'msi_score': 'MSI Score',
    'tmb_score': 'TMB Score',
}

for CANCER in CANCERS:
    print(f'\n{"="*60}\n  {CANCER}\n{"="*60}')

    cfg = CANCER_CONFIGS[CANCER]

    # load data
    mut_raw, clin_p_raw, clin_s_raw, ppi_path, cosmic_genes = load_raw_data(ROOT, cfg)
    seed_dict, mut_count_dict = parse_mutations(mut_raw)
    clinical, seed_by_patient, mut_count_by_patient = extract_clinical(
        clin_p_raw, clin_s_raw, seed_dict, mut_count_dict, cfg)

    nodes, node_to_idx, W, _ = build_transition_matrix(ppi_path)

    # best hyperparameters
    GS_CSV = os.path.join(RESULTS, f'Pipeline_AllPPI_{CANCER}/grid_search_results.csv')
    AK_CSV = os.path.join(RESULTS, f'Pipeline_AllPPI_{CANCER}/alpha_k_grid_results.csv')

    gs_best      = pd.read_csv(GS_CSV).sort_values('val_c', ascending=False).iloc[0]
    ak_best      = pd.read_csv(AK_CSV).sort_values('val_mean', ascending=False).iloc[0]
    HIDDEN1      = int(gs_best['hidden1'])
    HIDDEN2      = int(gs_best['hidden2'])
    OUTPUT_NODES = int(gs_best['output_nodes'])
    LR           = float(gs_best['LR'])
    L2           = float(gs_best['L2'])
    best_alpha   = float(ak_best['alpha'])
    best_K       = int(ak_best['K'])

    print(f'  Arch={HIDDEN1},{HIDDEN2},{OUTPUT_NODES}  LR={LR:.0e}  '
          f'L2={L2:.4f}  alpha={best_alpha}  K={best_K}')

    # RWR + feature selection
    all_idx = np.arange(len(nodes))
    rwr_all = compute_rwr_scores(clinical, W, nodes, node_to_idx, seed_by_patient,
                                 feat_idx=all_idx, alpha=best_alpha)
    df_all  = make_feature_df(rwr_all, nodes, clinical)

    X_all    = df_all[nodes].values.astype('float32')
    clin_arr = df_all[CLIN_FEATURES].values.astype('float32')
    yt       = df_all['OS_MONTHS'].values.astype('float32')
    ye       = df_all['OS_STATUS'].values.astype('float32')

    X, _, _, _ = select_top_k(X_all, nodes, yt, ye, k=best_K)
    n_genes    = best_K
    n_patients = X.shape[0]
    print(f'  Patients={n_patients}  Genes={n_genes}  Events={ye.mean():.1%}')

    # clinical feature KM (AGE, MSI, TMB) — matches Cell 48 exactly
    all_clinical_km[CANCER] = {}
    for feat in CLIN_FEAT_LABELS:
        tmp = df_all[['OS_MONTHS', 'OS_STATUS', feat]].dropna()
        median_val = float(tmp[feat].median())
        high = tmp[feat] >= median_val
        low  = ~high
        lr   = logrank_test(tmp.loc[high, 'OS_MONTHS'], tmp.loc[low, 'OS_MONTHS'],
                            event_observed_A=tmp.loc[high, 'OS_STATUS'],
                            event_observed_B=tmp.loc[low, 'OS_STATUS'])
        kmh = KaplanMeierFitter().fit(tmp.loc[high, 'OS_MONTHS'],
                                       tmp.loc[high, 'OS_STATUS'])
        kml = KaplanMeierFitter().fit(tmp.loc[low,  'OS_MONTHS'],
                                       tmp.loc[low,  'OS_STATUS'])
        all_clinical_km[CANCER][feat] = (kmh, kml, lr.p_value,
                                          median_val, high.sum(), low.sum())
        # save per-patient data so replot_km_1x4.py can regenerate without retraining
        clin_csv_path = os.path.join(KM_DIR, f'km_clinical_scores_{CANCER}_{feat}.csv')
        tmp.assign(risk_group=np.where(high, 'High', 'Low')).to_csv(clin_csv_path, index=False)
        print(f'  {feat:12s}: median={median_val:.2f}  p={lr.p_value:.4f}')
        clinical_km_summary.append({
            'cancer':        CANCER,
            'feature':       feat,
            'feature_label': CLIN_FEAT_LABELS[feat],
            'n_total':       len(tmp),
            'median_val':    round(median_val, 4),
            'high_n':        int(high.sum()),
            'low_n':         int(low.sum()),
            'pval':          round(lr.p_value, 6),
            'significance':  '***' if lr.p_value<0.001 else ('**' if lr.p_value<0.01
                              else ('*' if lr.p_value<0.05 else 'ns')),
            'high_med_surv': f'{kmh.median_survival_time_:.1f}m' if np.isfinite(kmh.median_survival_time_) else 'NR',
            'low_med_surv':  f'{kml.median_survival_time_:.1f}m' if np.isfinite(kml.median_survival_time_) else 'NR',
        })

    mk = dict(input_nodes=n_genes, hidden_nodes1=HIDDEN1,
              hidden_nodes2=HIDDEN2, output_nodes=OUTPUT_NODES,
              n_clin=len(CLIN_FEATURES))

    # unified loop — follows notebook exactly
    # Normalise on real data, then zero tensors per variant during training
    # (training-time masking, not inference-time).

    VARIANT_CONFIGS = {
        'CoxMLP':       dict(x_zero=False, cl_zero=False, l2=0.01),
        'ClinicalOnly': dict(x_zero=True,  cl_zero=False, l2=L2),
        'RWROnly':      dict(x_zero=False, cl_zero=True,  l2=0.2),
    }

    all_rows       = {}
    km_score_sum   = {v: np.zeros(n_patients, dtype=np.float64) for v in VARIANT_CONFIGS}
    km_score_count = {v: np.zeros(n_patients, dtype=np.int32)   for v in VARIANT_CONFIGS}
    abl_results    = {k: [] for k in MASK_CONFIGS}

    for variant, vcfg in VARIANT_CONFIGS.items():
        print(f'\n  Training {variant} (l2={vcfg["l2"]}) ...')
        test_cis_v, val_cis_v = [], []

        for rep in range(N_REPEATS):
            rng = np.random.RandomState(rep * 7 + 13)
            tr, va, te = censoring_time_balanced_split(yt, ye, TEST_FRAC, VAL_FRAC, rng)

            # normalise on real data first
            Xtr, Xva, Xte = normalise_genomic(X[tr], X[va], X[te])
            ctr, cva, cte = normalise_clinical(clin_arr[tr], clin_arr[va], clin_arr[te])

            x_tr, cl_tr, yt_tr, ye_tr = prep(Xtr, ctr, yt[tr], ye[tr])
            x_va, cl_va, yt_va, ye_va = prep(Xva, cva, yt[va], ye[va])
            x_te, cl_te, yt_te, ye_te = prep(Xte, cte, yt[te], ye[te])

            # zero tensors during training — matches notebook exactly
            if vcfg['x_zero']:
                x_tr = torch.zeros_like(x_tr)
                x_va = torch.zeros_like(x_va)
                x_te = torch.zeros_like(x_te)
            if vcfg['cl_zero']:
                cl_tr = torch.zeros_like(cl_tr)
                cl_va = torch.zeros_like(cl_va)
                cl_te = torch.zeros_like(cl_te)

            return_net = (variant == 'CoxMLP')
            if return_net:
                tc, val_c, raw, net = train_one(mk, x_tr, cl_tr, yt_tr, ye_tr,
                                                x_va, cl_va, yt_va, ye_va,
                                                x_te, cl_te, yt_te, ye_te,
                                                LR, vcfg['l2'], seed=rep,
                                                return_net=True)
                # inference-time ablation on CoxMLP only
                for mask_label, kwargs in MASK_CONFIGS.items():
                    abl_results[mask_label].append(
                        eval_with_mask(net, x_te, cl_te, yt_te, ye_te, **kwargs))
            else:
                tc, val_c, raw = train_one(mk, x_tr, cl_tr, yt_tr, ye_tr,
                                           x_va, cl_va, yt_va, ye_va,
                                           x_te, cl_te, yt_te, ye_te,
                                           LR, vcfg['l2'], seed=rep,
                                           return_net=False)

            test_cis_v.append(tc)
            val_cis_v.append(val_c)

            sort_idx   = np.argsort(yt[te])[::-1]
            orig_order = np.argsort(sort_idx)
            km_score_sum[variant][te]   += raw[orig_order]
            km_score_count[variant][te] += 1

            if (rep + 1) % 5 == 0:
                print(f'    rep {rep+1:2d}/{N_REPEATS}  C={tc:.4f}')

        all_rows[variant] = [{'rep': r+1, 'val_c': round(v,4), 'test_c': round(t,4)}
                             for r,(t,v) in enumerate(zip(test_cis_v, val_cis_v))]

    mean_c = float(np.mean([r['test_c'] for r in all_rows['CoxMLP']]))

    print(f'\n  {"="*50}  {CANCER} C-INDEX SUMMARY')
    for m, rows in all_rows.items():
        cs = [r['test_c'] for r in rows]
        print(f'  {m:14s}  {np.mean(cs):.4f} ± {np.std(cs):.4f}')

    for m, rows in all_rows.items():
        pd.DataFrame(rows).assign(model=m, cancer=CANCER).to_csv(
            os.path.join(RESULTS, f'{CANCER.lower()}_allppi_{m.lower()}.csv'), index=False)

    means_ci, stds_ci, ci_data = plot_cindex(CANCER, all_rows)
    all_ci_data[CANCER] = ci_data
    print(f'\n  Mean C-index (CoxMLP): {mean_c:.4f}')

    # ablation plot
    means_abl, stds_abl, drops = plot_ablation(CANCER, abl_results)

    print(f'\n  {"Mask condition":<25}  {"Mean":>7}  {"Std":>7}  {"Drop":>7}')
    print(f'  {"-"*55}')
    for k, m, s, d in zip(abl_results.keys(), means_abl, stds_abl, drops):
        print(f'  {k:<25}  {m:7.4f}  {s:7.4f}  {d:+7.4f}')
        ablation_summary.append({'cancer': CANCER, 'mask': k,
                                  'mean_c': round(m, 4), 'std_c': round(s, 4),
                                  'delta': round(d, 4)})

    # save raw per-rep inference ablation values
    abl_rows = []
    for mask_label, vals in abl_results.items():
        for rep_i, v in enumerate(vals):
            abl_rows.append({'rep': rep_i + 1, 'mask': mask_label,
                             'test_c': round(v, 4), 'cancer': CANCER})
    abl_raw_path = os.path.join(RESULTS, f'{CANCER.lower()}_inference_ablation.csv')
    pd.DataFrame(abl_rows).to_csv(abl_raw_path, index=False)
    print(f'  Inference ablation raw data → {abl_raw_path}')

    # pooled KM — all 3 model variants
    def fmt(v): return f'{v:.1f}m' if np.isfinite(v) else 'NR'
    all_km_data[CANCER] = {}
    mm = 1 / 25.4

    for mk_, cfg_ in KM_MODEL_CONFIGS.items():
        never = km_score_count[mk_] == 0
        if never.any():
            print(f'  {mk_}: {never.sum()} patients never tested — excluded')
        mask_      = ~never
        avg_scores = km_score_sum[mk_][mask_] / km_score_count[mk_][mask_]
        t_all, e_all = yt[mask_], ye[mask_]
        thresh     = float(np.median(avg_scores))
        high       = avg_scores >= thresh
        low        = ~high

        pval, med_h, med_l, kmh, kml = km_stats(
            t_all[high], t_all[low], e_all[high], e_all[low])

        n_high, n_low = int(high.sum()), int(low.sum())
        print(f'  {mk_:14s}: p={pval:.2e}  High={fmt(med_h)} (n={n_high})  Low={fmt(med_l)} (n={n_low})')

        all_km_data[CANCER][mk_] = (kmh, kml, pval, n_high, n_low)

        # save per-patient data so replot_km_1x4.py can regenerate without retraining
        model_csv_path = os.path.join(KM_DIR, f'km_scores_{CANCER}_{mk_}.csv')
        pd.DataFrame({
            'OS_MONTHS':      t_all,
            'OS_STATUS':      e_all,
            'avg_risk_score': avg_scores,
            'risk_group':     np.where(high, 'High', 'Low'),
        }).to_csv(model_csv_path, index=False)

        # individual KM plot
        fig, ax = plt.subplots(figsize=(44*mm, 62*mm), dpi=300)
        plot_km(ax, kmh, kml, pval, CANCER, n_high, n_low,
                letter=None, show_ylabel=True)
        plt.tight_layout(pad=0.3)
        fpath_pdf = os.path.join(KM_DIR, f'km_{mk_.lower()}_{CANCER.lower()}.pdf')
        fig.savefig(fpath_pdf, format='pdf', bbox_inches='tight', dpi=300)
        fig.savefig(fpath_pdf.replace('.pdf', '.png'), format='png',
                    bbox_inches='tight', dpi=300)
        plt.close()

        km_summary.append({
            'cancer':        CANCER,
            'model':         KM_MODEL_CONFIGS[mk_]['label'],
            'mean_cindex':   round(mean_c, 4),
            'n_patients':    int(mask_.sum()),
            'pval':          round(pval, 6),
            'high_n':        n_high,
            'low_n':         n_low,
            'high_med_surv': fmt(med_h),
            'low_med_surv':  fmt(med_l),
        })
    print(f'  KM plots saved → {KM_DIR}')

    for m in ['CoxMLP', 'ClinicalOnly', 'RWROnly']:
        cs = [r['test_c'] for r in all_rows[m]]
        cindex_summary.append({
            'cancer': CANCER, 'model': m,
            'mean_c': round(np.mean(cs), 4),
            'std_c':  round(np.std(cs),  4),
        })

# combined 1×4 KM figures — one per model variant
print('\nGenerating combined KM figures ...')
for mk_, cfg_ in KM_MODEL_CONFIGS.items():
    make_km_figure(all_km_data, mk_, suffix=mk_.lower())

# clinical feature KM figures — 3 features × 4 cancers
print('\nGenerating clinical KM figures (AGE, MSI, TMB) ...')
mm    = 1 / 25.4
FIG_W = 174 * mm
FIG_H =  72 * mm   # was 62 — extra 10mm prevents legend/x-axis overlap

for feat, feat_label in CLIN_FEAT_LABELS.items():
    fig, axes = plt.subplots(1, 4, figsize=(FIG_W, FIG_H), dpi=300, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.995, top=0.88, bottom=0.22, wspace=0.08)

    for col, cancer in enumerate(CANCERS):
        ax = axes[col]
        kmh, kml, pval, median_val, n_high, n_low = all_clinical_km[cancer][feat]
        ps = f'p = {pval:.2e}'

        kmh.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.15,
                                    color=KM_COLORS['High Risk'], linewidth=0.8,
                                    label=f'≥{median_val:.2f} (n={n_high})')
        kml.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.15,
                                    color=KM_COLORS['Low Risk'],  linewidth=0.8,
                                    label=f'<{median_val:.2f}  (n={n_low})')

        ax.set_title(cancer, fontsize=7, pad=3)
        ax.set_xlabel('OS Months', labelpad=2)
        if col == 0:
            ax.set_ylabel('Survival Probability', labelpad=3)
        else:
            ax.set_ylabel('')

        ax.set_ylim(0, 1.05)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.text(0.04, 0.97, PANEL_LETTERS[cancer], transform=ax.transAxes,
                fontsize=8, fontweight='bold', va='top', ha='left', color='#222222')

        ax.text(0.97, 0.97, ps, transform=ax.transAxes,
                fontsize=5.5, va='top', ha='right', color='#222222')
        ax.legend(loc='lower left', frameon=False, fontsize=5.5,
                  handlelength=0.8, handletextpad=0.3,
                  borderpad=0.1, labelspacing=0.15)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.4)
        ax.spines['bottom'].set_linewidth(0.4)
        if col > 0:
            ax.tick_params(axis='y', which='both', left=False)
        ax.tick_params(axis='both', width=0.4, length=2.0, pad=1.5)
        ax.set_facecolor('white')

    out_pdf = os.path.join(KM_DIR, f'km_clinical_{feat.lower()}_1x4.pdf')
    fig.savefig(out_pdf, format='pdf', bbox_inches='tight', dpi=300)
    fig.savefig(out_pdf.replace('.pdf', '.png'), format='png',
                bbox_inches='tight', dpi=300)
    plt.close()
    print(f'  Saved → {out_pdf}')

# combined 1×4 boxplot — matches generate_figure4.py exactly
print('\nGenerating combined C-index figures ...')
make_cindex_figure(all_ci_data, show_dots=True)
make_cindex_figure(all_ci_data, show_dots=False)

# save summaries
abl_csv  = os.path.join(RESULTS, 'ablation_summary.csv')
km_csv   = os.path.join(RESULTS, 'km_summary.csv')
ci_csv   = os.path.join(RESULTS, 'cindex_summary.csv')
clin_csv = os.path.join(RESULTS, 'clinical_km_summary.csv')

pd.DataFrame(ablation_summary).to_csv(abl_csv,  index=False)
pd.DataFrame(km_summary).to_csv(km_csv,          index=False)
pd.DataFrame(cindex_summary).to_csv(ci_csv,      index=False)
pd.DataFrame(clinical_km_summary).to_csv(clin_csv, index=False)

print(f'\n\nC-index summary      → {ci_csv}')
print(f'Ablation summary     → {abl_csv}')
print(f'KM summary (models)  → {km_csv}')
print(f'KM summary (clinical)→ {clin_csv}')
print('\nC-index Summary:')
print(pd.DataFrame(cindex_summary).to_string(index=False))
print('\nKM Model Summary:')
print(pd.DataFrame(km_summary).to_string(index=False))
print('\nClinical KM Summary:')
print(pd.DataFrame(clinical_km_summary).to_string(index=False))
