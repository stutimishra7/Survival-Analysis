"""
RCoxNet end-to-end pipeline.

Performs a single joint grid search over all hyperparameters
(α × K × architecture × LR × L2), then trains the final model with
the best combination.  This is the approach described in the paper.

Usage
-----
    from rcoxnet.pipeline import run_pipeline
    result = run_pipeline('BRCA')
    print(result['test_ci_mean'], result['test_ci_std'])

Or from the command line:
    python run.py --cancer BRCA
"""

import os
import copy

import numpy as np
import pandas as pd
import torch
import torch.optim as optim

from rcoxnet.config import (
    CANCER_CONFIGS, CLIN_FEATURES,
    ALPHA_VALUES, GRID_ARCHITECTURES, GRID_L2, GRID_LR,
    N_REPEATS, TEST_FRAC, VAL_FRAC,
    GS_EPOCHS, GS_CHECKPOINT, GS_GRAD_CLIP,
)
from rcoxnet.clinical import extract_clinical
from rcoxnet.data_loader import load_raw_data
from rcoxnet.features import make_feature_df, logrank_pvalues
from rcoxnet.model import CoxPhRWRNet, cox_loss, c_index
from rcoxnet.mutations import parse_mutations
from rcoxnet.network import build_transition_matrix
from rcoxnet.preprocessing import (
    censoring_time_balanced_split,
    normalise_genomic,
    normalise_clinical,
    prep,
)
from rcoxnet.rwr import compute_rwr_scores
from rcoxnet.train import train_and_eval

# Grid search uses fewer repeats for speed; final evaluation uses N_REPEATS
GS_REPEATS = 3

def load_data(cancer: str, root: str) -> dict:
    """Load mutations, clinical info, and return a unified patient table."""
    cfg = CANCER_CONFIGS[cancer]
    mut_raw, clin_p_raw, clin_s_raw, ppi_path, _ = load_raw_data(root, cfg)
    seed_dict, mut_count_dict = parse_mutations(mut_raw)
    clinical, seed_by_patient, _ = extract_clinical(
        clin_p_raw, clin_s_raw, seed_dict, mut_count_dict, cfg)
    print(f'  Patients: {len(clinical)}  |  PPI file: {os.path.basename(ppi_path)}')
    return dict(clinical=clinical, seed_by_patient=seed_by_patient,
                ppi_path=ppi_path, cfg=cfg)

def build_network(ppi_path: str) -> dict:
    """Build the PPI transition matrix from the ConsensusPathDB file."""
    nodes, node_to_idx, W, _ = build_transition_matrix(ppi_path)
    print(f'  PPI genes: {len(nodes)}')
    return dict(nodes=nodes, node_to_idx=node_to_idx, W=W,
                all_idx=np.arange(len(nodes)))

def joint_grid_search(
    device: torch.device,
    dtype,
    clinical: pd.DataFrame,
    seed_by_patient: dict,
    network: dict,
) -> dict:
    """
    Joint grid search over α × K × architecture × LR × L2.

    RWR scores are computed once per α value and reused across all
    K / architecture / LR / L2 combinations, which keeps runtime manageable.

    Returns:
    best : dict   Keys: alpha, K, hidden1, hidden2, output_nodes, lr, l2, val_c
    rows : list   All evaluated combinations (saved to CSV by the caller)
    """
    nodes      = network['nodes']
    node_to_idx= network['node_to_idx']
    W          = network['W']
    all_idx    = network['all_idx']
    n_genes    = len(nodes)

    k_values   = _k_schedule(n_genes)
    n_combos   = (len(ALPHA_VALUES) * len(k_values)
                  * len(GRID_ARCHITECTURES) * len(GRID_L2) * len(GRID_LR))

    print(f'  α        : {ALPHA_VALUES}')
    print(f'  K values : {k_values}')
    print(f'  Arch     : {GRID_ARCHITECTURES}')
    print(f'  LR / L2  : {GRID_LR} / {GRID_L2}')
    print(f'  Combos   : {n_combos} × {GS_REPEATS} repeats = '
          f'{n_combos * GS_REPEATS} training runs')

    best, rows, combo = None, [], 0

    for alpha in ALPHA_VALUES:
        # --- compute RWR once per α, reuse for all K below ---
        rwr   = compute_rwr_scores(clinical, W, nodes, node_to_idx,
                                   seed_by_patient, feat_idx=all_idx, alpha=alpha)
        df_a  = make_feature_df(rwr, nodes, clinical)
        X_a   = df_a[nodes].values.astype(np.float32)
        yt_a  = df_a['OS_MONTHS'].values.astype(np.float32)
        ye_a  = df_a['OS_STATUS'].values.astype(np.float32)
        cl_a  = df_a[CLIN_FEATURES].values.astype(np.float32)

        # log-rank p-values computed once per α, reused for every K
        pvals  = logrank_pvalues(X_a, yt_a, ye_a)
        order  = np.argsort(pvals)

        for k in k_values:
            X_k = X_a[:, order[:k]] if k < n_genes else X_a

            for h1, h2, out in GRID_ARCHITECTURES:
                for l2 in GRID_L2:
                    for lr in GRID_LR:
                        combo += 1
                        val_c  = _eval_combo(device, dtype, X_k, cl_a,
                                             yt_a, ye_a, h1, h2, out,
                                             lr, l2, combo)
                        row = dict(alpha=alpha, K=k, hidden1=h1, hidden2=h2,
                                   output_nodes=out, LR=lr, L2=l2, val_c=val_c)
                        rows.append(row)
                        print(f'  [{combo:04d}/{n_combos}] '
                              f'α={alpha} K={k:4d} '
                              f'arch=({h1},{h2},{out}) '
                              f'L2={l2:.3f} LR={lr:.0e} | val_c={val_c:.4f}')
                        if best is None or val_c > best['val_c']:
                            best = row.copy()

    print(f'\n  Best → α={best["alpha"]}  K={best["K"]}  '
          f'arch=({best["hidden1"]},{best["hidden2"]},{best["output_nodes"]})  '
          f'LR={best["LR"]:.0e}  L2={best["L2"]:.3f}  val_c={best["val_c"]:.4f}')
    return best, rows

def train_final(
    device: torch.device,
    dtype,
    clinical: pd.DataFrame,
    seed_by_patient: dict,
    network: dict,
    best: dict,
) -> tuple:
    """
    Train the final model with the best hyperparameters over N_REPEATS splits.

    Returns test_cis, val_cis — one value per repeat.
    """
    nodes       = network['nodes']
    node_to_idx = network['node_to_idx']
    W           = network['W']
    all_idx     = network['all_idx']

    rwr    = compute_rwr_scores(clinical, W, nodes, node_to_idx,
                                seed_by_patient, feat_idx=all_idx,
                                alpha=best['alpha'])
    df     = make_feature_df(rwr, nodes, clinical)
    X      = df[nodes].values.astype(np.float32)
    yt     = df['OS_MONTHS'].values.astype(np.float32)
    ye     = df['OS_STATUS'].values.astype(np.float32)
    cl     = df[CLIN_FEATURES].values.astype(np.float32)

    pvals  = logrank_pvalues(X, yt, ye)
    order  = np.argsort(pvals)
    k      = best['K']
    X_top  = X[:, order[:k]] if k < len(nodes) else X

    model_kwargs = dict(
        input_nodes  = X_top.shape[1],
        hidden_nodes1= int(best['hidden1']),
        hidden_nodes2= int(best['hidden2']),
        output_nodes = int(best['output_nodes']),
        n_clin       = len(CLIN_FEATURES),
    )
    test_cis, val_cis = train_and_eval(
        device, dtype, X_top, cl, yt, ye,
        model_kwargs, lr=best['lr'], l2=best['l2'])
    return test_cis, val_cis

def run_pipeline(cancer: str, root: str = None) -> dict:
    """
    Run the full RCoxNet pipeline for one cancer type.

    Steps
    -----
    1. Load mutation + clinical data
    2. Build PPI transition matrix
    3. Joint grid search  (α × K × arch × LR × L2)
    4. Train final model  (20 independent random splits)
    5. Save results to Results/Pipeline_<CANCER>/

    Args:
    cancer : str   One of BRCA / LUNG / GBM / OV.
    root   : str   Repository root (default: parent of this file).

    Returns:
    dict with keys: cancer, best_alpha, best_K, hidden1, hidden2,
    output_nodes, lr, l2, test_ci_mean, test_ci_std, val_ci_mean, val_ci_std
    """
    assert cancer in CANCER_CONFIGS, \
        f'Unknown cancer: {cancer!r}.  Choose from {list(CANCER_CONFIGS)}.'

    root = root or os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    cfg  = CANCER_CONFIGS[cancer]

    results_dir = os.path.join(root, cfg['results_dir'])
    os.makedirs(results_dir, exist_ok=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dtype  = torch.FloatTensor

    _banner(cancer, device)

    # Step 1: data
    _section('Step 1 — Load data')
    data = load_data(cancer, root)

    # Step 2: network
    _section('Step 2 — Build PPI transition matrix')
    network = build_network(data['ppi_path'])

    # Step 3: joint grid search
    _section('Step 3 — Joint grid search  (α × K × arch × LR × L2)')
    best, rows = joint_grid_search(
        device, dtype, data['clinical'], data['seed_by_patient'], network)

    pd.DataFrame(rows) \
      .sort_values('val_c', ascending=False) \
      .reset_index(drop=True) \
      .to_csv(os.path.join(results_dir, 'joint_grid_results.csv'), index=False)

    # rename LR/L2 keys to lr/l2 for train_final
    best['lr'] = best.pop('LR')
    best['l2'] = best.pop('L2')

    # Step 4: final model
    _section(f'Step 4 — Final model  '
             f'(α={best["alpha"]}, K={best["K"]}, '
             f'arch=({best["hidden1"]},{best["hidden2"]},{best["output_nodes"]}))')
    test_cis, val_cis = train_final(
        device, dtype, data['clinical'], data['seed_by_patient'], network, best)

    pd.DataFrame({'test_ci': test_cis, 'val_ci': val_cis}) \
      .to_csv(os.path.join(results_dir, 'final_model_cindex.csv'), index=False)

    _summary(cancer, best, test_cis, val_cis, results_dir)

    return {
        'cancer':       cancer,
        'best_alpha':   best['alpha'],
        'best_K':       best['K'],
        'hidden1':      best['hidden1'],
        'hidden2':      best['hidden2'],
        'output_nodes': best['output_nodes'],
        'lr':           best['lr'],
        'l2':           best['l2'],
        'test_ci_mean': float(np.mean(test_cis)),
        'test_ci_std':  float(np.std(test_cis)),
        'val_ci_mean':  float(np.mean(val_cis)),
        'val_ci_std':   float(np.std(val_cis)),
    }

def _k_schedule(n_total: int) -> list:
    """100, 200, …, 800, 1000, 2000, … n_total."""
    vals = list(range(100, min(801, n_total + 1), 100))
    if n_total > 800:
        nxt = 1000
        while nxt < n_total:
            vals.append(nxt)
            nxt += 1000
    if not vals or vals[-1] != n_total:
        vals.append(n_total)
    return vals

def _eval_combo(device, dtype, X, cl, yt, ye,
                h1, h2, out, lr, l2, seed_offset):
    """Train one hyperparameter combo for GS_REPEATS splits, return mean val C-index."""
    n_clin = len(CLIN_FEATURES)
    model_kwargs = dict(input_nodes=X.shape[1], hidden_nodes1=h1,
                        hidden_nodes2=h2, output_nodes=out, n_clin=n_clin)
    val_cis = []
    for rep in range(GS_REPEATS):
        rng = np.random.RandomState(rep * 7 + 13 + seed_offset)
        tr, va, te = censoring_time_balanced_split(yt, ye, TEST_FRAC, VAL_FRAC, rng)

        X_tr, X_va, X_te    = normalise_genomic(X[tr], X[va], X[te])
        cl_tr, cl_va, cl_te = normalise_clinical(cl[tr], cl[va], cl[te])

        x_tr, cl_tr_t, yt_tr, ye_tr = prep(device, dtype, X_tr, cl_tr, yt[tr], ye[tr])
        x_va, cl_va_t, yt_va, ye_va = prep(device, dtype, X_va, cl_va, yt[va], ye[va])
        x_te, cl_te_t, yt_te, ye_te = prep(device, dtype, X_te, cl_te, yt[te], ye[te])

        torch.manual_seed(rep + seed_offset)
        np.random.seed(rep + seed_offset)
        net = CoxPhRWRNet(**model_kwargs).to(device)
        opt = optim.Adam(net.parameters(), lr=lr, weight_decay=l2)
        best_val, best_st = 0.0, copy.deepcopy(net.state_dict())

        for ep in range(GS_EPOCHS):
            net.train(); opt.zero_grad()
            loss = cox_loss(net(x_tr, cl_tr_t), yt_tr, ye_tr)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(net.parameters(), GS_GRAD_CLIP)
            opt.step()
            if ep % GS_CHECKPOINT == 0:
                net.eval()
                with torch.no_grad():
                    vc = c_index(net(x_va, cl_va_t), yt_va, ye_va).item()
                if vc > best_val:
                    best_val, best_st = vc, copy.deepcopy(net.state_dict())

        val_cis.append(best_val)
    return float(np.mean(val_cis))

def _banner(cancer, device):
    line = '=' * 58
    print(f'\n{line}')
    print(f'  RCoxNet — {cancer}')
    print(f'  Device         : {device}')
    print(f'  Clinical input : {CLIN_FEATURES}')
    print(f'  Grid search    : α × K × arch × LR × L2  (joint)')
    print(f'  GS repeats     : {GS_REPEATS}  |  Final repeats: {N_REPEATS}')
    print(f'{line}\n')

def _section(title):
    print(f'\n{title}')

def _summary(cancer, best, test_cis, val_cis, results_dir):
    line = '=' * 58
    print(f'\n{line}')
    print(f'  {cancer} — done')
    print(f'  α={best["alpha"]}  K={best["K"]}  '
          f'arch=({best["hidden1"]},{best["hidden2"]},{best["output_nodes"]})')
    print(f'  LR={best["lr"]:.0e}  L2={best["l2"]:.3f}')
    print(f'  Test C-index : {np.mean(test_cis):.4f} ± {np.std(test_cis):.4f}')
    print(f'  Val  C-index : {np.mean(val_cis):.4f} ± {np.std(val_cis):.4f}')
    print(f'  Saved to     : {results_dir}')
    print(f'{line}\n')
