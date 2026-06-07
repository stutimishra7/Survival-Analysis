"""
bench_common.py — shared data preparation for all benchmark models.

Feature space: RWR scores recomputed with each cancer's best (alpha, K),
then top-K genes selected by log-rank p-value — exactly the same feature
space used by RCoxNet AllPPI. Only the model architecture differs.

Splits: 70 / 10 / 20  stratified by event-rate, seed = 42 + rep  (matches RCoxNet).
Normalisation: log10 + z-score fitted on train only, applied to val/test.
"""

import os, sys
import numpy as np
import pandas as pd

BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT      = os.path.dirname(BENCH_DIR)
sys.path.insert(0, ROOT)

from rcoxnet.config import CANCER_CONFIGS, CLIN_FEATURES, META_COLS, TEST_FRAC, VAL_FRAC
from rcoxnet.data_loader import load_raw_data
from rcoxnet.mutations import parse_mutations
from rcoxnet.clinical import extract_clinical
from rcoxnet.network import build_transition_matrix
from rcoxnet.rwr import compute_rwr_scores

# Best (alpha, K) per cancer — from RCoxNet AllPPI grid search results.
# Source: Results/Pipeline_AllPPI_{CANCER}/alpha_k_grid_results.csv, top row.
CANCERS = {
    'BRCA': {'best_alpha': 0.8, 'top_k': 300},
    'GBM':  {'best_alpha': 0.5, 'top_k': 800},
    'OV':   {'best_alpha': 0.5, 'top_k': 300},
    'LUNG': {'best_alpha': 0.8, 'top_k': 600},
}

_DATA_CACHE = {}

# Log-rank p-values (no lifelines dependency)

def _logrank_pvalues(gene_arr, yt, ye):
    from scipy.stats import chi2
    pvals = []
    ye_b  = ye.astype(bool)
    for j in range(gene_arr.shape[1]):
        col  = gene_arr[:, j]
        mask = col > np.median(col)
        g1   = np.where( mask)[0]
        g2   = np.where(~mask)[0]
        if len(g1) < 5 or len(g2) < 5:
            pvals.append(1.0); continue
        times  = np.unique(yt[ye_b])
        O1 = E1 = O2 = E2 = 0.0
        for t in times:
            at1 = (yt[g1] >= t).sum()
            at2 = (yt[g2] >= t).sum()
            n   = at1 + at2
            if n < 2: continue
            d1  = ((yt[g1] == t) & ye_b[g1]).sum()
            d2  = ((yt[g2] == t) & ye_b[g2]).sum()
            d   = d1 + d2
            O1 += d1;  O2 += d2
            E1 += d * at1 / n
            E2 += d * at2 / n
        denom = E1 + E2
        if denom < 1e-12:
            pvals.append(1.0); continue
        stat = (O1-E1)**2/E1 + (O2-E2)**2/E2 if (E1 > 0 and E2 > 0) else 0.0
        pvals.append(float(1 - chi2.cdf(stat, df=1)))
    return np.array(pvals)

def _select_top_k(gene_arr, gene_names, yt, ye, k):
    pvals      = _logrank_pvalues(gene_arr, yt, ye)
    sorted_idx = np.argsort(pvals)
    top_idx    = sorted_idx[:k] if k < len(gene_names) else sorted_idx
    return gene_arr[:, top_idx], [gene_names[i] for i in top_idx]

# RWR feature matrix (cached per cancer)

def _build_rwr_df(cancer):
    if cancer in _DATA_CACHE:
        return _DATA_CACHE[cancer]

    cfg   = CANCER_CONFIGS[cancer]
    alpha = CANCERS[cancer]['best_alpha']
    K     = CANCERS[cancer]['top_k']

    print(f"  [{cancer}] Loading raw data & building PPI...", flush=True)
    mut_raw, clin_p_raw, clin_s_raw, ppi_path, _ = load_raw_data(ROOT, cfg)
    seed_dict, mut_count_dict = parse_mutations(mut_raw)
    clinical, seed_by_patient, _ = extract_clinical(
        clin_p_raw, clin_s_raw, seed_dict, mut_count_dict, cfg)

    nodes, node_to_idx, W, _ = build_transition_matrix(ppi_path)
    all_idx = np.arange(len(nodes))

    print(f"  [{cancer}] Computing RWR (alpha={alpha}, {len(nodes)} genes)...", flush=True)
    rwr_scores = compute_rwr_scores(
        clinical, W, nodes, node_to_idx, seed_by_patient,
        feat_idx=all_idx, alpha=alpha)

    rwr_df = pd.DataFrame(rwr_scores, columns=nodes)
    rwr_df.insert(0, 'SAMPLE_ID', clinical['SAMPLE_ID'].values)
    merged = rwr_df.merge(clinical[list(META_COLS)], on='SAMPLE_ID', how='inner')
    merged = merged.dropna(subset=['OS_MONTHS', 'OS_STATUS']).reset_index(drop=True)
    merged['OS_STATUS'] = merged['OS_STATUS'].astype(int)
    has_signal = (merged[nodes].values > 0).any(axis=1)
    merged = merged[has_signal].reset_index(drop=True)

    yt     = merged['OS_MONTHS'].values.astype(np.float32)
    ye     = merged['OS_STATUS'].values.astype(np.float32)
    X_all  = merged[nodes].values.astype(np.float32)
    clin_cols = [c for c in CLIN_FEATURES if c in merged.columns]

    print(f"  [{cancer}] Log-rank top-{K} selection ({len(nodes)} genes)...", flush=True)
    _, selected = _select_top_k(X_all, nodes, yt, ye, k=K)
    print(f"  [{cancer}] Selected {len(selected)} genes.", flush=True)

    result = {
        'df':        merged,
        'gene_cols': selected,
        'clin_cols': clin_cols,
    }
    _DATA_CACHE[cancer] = result
    return result

# Normalisation

def normalise_genomic(X_tr, *others):
    X_tr = np.log10(X_tr + 1e-12)
    mu   = X_tr.mean(axis=0, keepdims=True)
    sig  = X_tr.std( axis=0, keepdims=True)
    sig[sig < 1e-8] = 1.0
    X_tr = (X_tr - mu) / sig
    return [X_tr] + [(np.log10(X + 1e-12) - mu) / sig for X in others]

def normalise_clinical(cl_tr, *others):
    mu  = np.nanmean(cl_tr, axis=0, keepdims=True)
    sig = np.nanstd( cl_tr, axis=0, keepdims=True)
    sig[sig < 1e-8] = 1.0
    cl_tr = np.where(np.isnan(cl_tr), mu, cl_tr)
    return [(cl_tr - mu) / sig] + [
        ((np.where(np.isnan(c), mu, c)) - mu) / sig for c in others]

# Stratified split

def censoring_time_balanced_split(yt, ye, test_frac, val_frac, rng):
    yt = np.asarray(yt); ye = np.asarray(ye).astype(int)
    ev_idx = np.where(ye == 1)[0]
    if len(ev_idx) >= 8:
        q1, q2  = np.quantile(yt[ev_idx], [0.33, 0.66])
        ev_bins = np.digitize(yt[ev_idx], bins=[q1, q2], right=True)
    else:
        med     = np.median(yt[ev_idx]) if len(ev_idx) > 0 else np.median(yt)
        ev_bins = (yt[ev_idx] > med).astype(int)
    groups = [ev_idx[ev_bins == b] for b in np.unique(ev_bins)]
    groups.append(np.where(ye == 0)[0])
    tr_parts, va_parts, te_parts = [], [], []
    for g in groups:
        g = g.copy(); rng.shuffle(g)
        n    = len(g)
        n_te = max(1, int(round(n * test_frac))) if n >= 3 else max(0, n // 3)
        rem  = n - n_te
        n_va = max(1, int(round(rem * val_frac))) if rem >= 3 else max(0, rem // 3)
        te_parts.append(g[:n_te])
        va_parts.append(g[n_te:n_te + n_va])
        tr_parts.append(g[n_te + n_va:])
    tr = np.concatenate([x for x in tr_parts if len(x)])
    va = np.concatenate([x for x in va_parts if len(x)])
    te = np.concatenate([x for x in te_parts if len(x)])
    rng.shuffle(tr); rng.shuffle(va); rng.shuffle(te)
    return tr, va, te

# C-index (Harrell's C, pure numpy)

def cindex_numpy(yt, ye, risk):
    risk = np.asarray(risk, dtype=float).ravel()
    yt   = np.asarray(yt,   dtype=float).ravel()
    ye   = np.asarray(ye,   dtype=bool ).ravel()
    mask = np.isfinite(risk)
    yt, ye, risk = yt[mask], ye[mask], risk[mask]
    if ye.sum() < 1 or mask.sum() < 2:
        return np.nan
    concordant = discordant = 0
    for i in np.where(ye)[0]:
        diff_t    = yt[i] - yt
        diff_r    = risk[i] - risk
        comparable = (yt > yt[i]) | ye
        comparable[i] = False
        concordant += ((diff_t < 0) & (diff_r > 0) & comparable).sum()
        concordant += ((diff_t > 0) & (diff_r < 0) & comparable).sum()
        discordant += ((diff_t < 0) & (diff_r < 0) & comparable).sum()
        discordant += ((diff_t > 0) & (diff_r > 0) & comparable).sum()
    total = concordant + discordant
    return float(concordant / total) if total > 0 else np.nan

# Main API

def prepare_rep(cancer, rep):
    """
    Returns dict with 'tr', 'va', 'te', each a tuple:
        (X_genomic_norm, X_clinical_norm, yt, ye)

    Feature space is identical to RCoxNet AllPPI: RWR with best alpha,
    top-K genes by log-rank p-value on training patients.
    """
    data      = _build_rwr_df(cancer)
    df        = data['df']
    gene_cols = data['gene_cols']
    clin_cols = data['clin_cols']

    yt     = df['OS_MONTHS'].values.astype(np.float32)
    ye     = df['OS_STATUS'].values.astype(np.float32)
    X_gen  = df[gene_cols].values.astype(np.float32)
    X_clin = df[clin_cols].values.astype(np.float32)

    rng = np.random.RandomState(42 + rep)
    tr_idx, va_idx, te_idx = censoring_time_balanced_split(
        yt, ye, TEST_FRAC, VAL_FRAC, rng)

    Xg_tr, Xg_va, Xg_te = X_gen[tr_idx],  X_gen[va_idx],  X_gen[te_idx]
    Xc_tr, Xc_va, Xc_te = X_clin[tr_idx], X_clin[va_idx], X_clin[te_idx]
    yt_tr, yt_va, yt_te = yt[tr_idx],     yt[va_idx],     yt[te_idx]
    ye_tr, ye_va, ye_te = ye[tr_idx],     ye[va_idx],     ye[te_idx]

    Xg_tr, Xg_va, Xg_te = normalise_genomic( Xg_tr, Xg_va, Xg_te)
    Xc_tr, Xc_va, Xc_te = normalise_clinical(Xc_tr, Xc_va, Xc_te)

    return {
        'tr': (Xg_tr, Xc_tr, yt_tr, ye_tr),
        'va': (Xg_va, Xc_va, yt_va, ye_va),
        'te': (Xg_te, Xc_te, yt_te, ye_te),
    }
