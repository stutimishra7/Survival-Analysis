"""
Hyperparameter grid search and joint (alpha, K) grid search.

run_hyperparameter_grid_search
    Searches over architecture × learning-rate × L2 on a single train/val
    split and returns the best configuration.

run_alpha_k_grid_search
    For each combination of RWR restart probability (alpha) and number of
    COSMIC genes (K), trains and evaluates the model with 20 repeats and
    returns the (alpha, K) pair with the highest mean validation C-index.
"""

import copy

import numpy as np
import pandas as pd
import torch
import torch.optim as optim

from rcoxnet.config import (
    ALPHA_VALUES, CLIN_FEATURES,
    GRID_ARCHITECTURES, GRID_L2, GRID_LR,
    GS_CHECKPOINT, GS_EPOCHS, GS_GRAD_CLIP, GS_SEED,
    N_REPEATS, TEST_FRAC, VAL_FRAC,
)
from rcoxnet.features import logrank_pvalues, make_feature_df
from rcoxnet.model import CoxPhRWRNet, c_index, cox_loss
from rcoxnet.preprocessing import (
    censoring_time_balanced_split,
    normalise_clinical,
    normalise_genomic,
    prep,
)
from rcoxnet.rwr import compute_rwr_scores
from rcoxnet.train import train_and_eval

def run_hyperparameter_grid_search(
    device: torch.device,
    dtype,
    X: np.ndarray,
    clin_arr: np.ndarray,
    yt: np.ndarray,
    ye: np.ndarray,
    n_genes: int,
):
    """
    Search over architecture × LR × L2 on a single train/val split.

    Args:
    X        : (n_patients, n_genes)   RWR genomic features (top-K).
    clin_arr : (n_patients, 3)         [AGE, msi_score, tmb_score].
    yt, ye   : survival time and event indicators.
    n_genes  : int   number of genomic input features.

    Returns:
    best : dict   Keys: hidden1, hidden2, output_nodes, LR, L2, val_c.
    gs_df : pd.DataFrame   Full grid results sorted by val_c descending.
    """
    gs_rng = np.random.RandomState(GS_SEED)
    tr, va, _ = censoring_time_balanced_split(yt, ye, TEST_FRAC, VAL_FRAC, gs_rng)

    X_tr, X_va   = normalise_genomic(X[tr], X[va])
    cl_tr, cl_va = normalise_clinical(clin_arr[tr], clin_arr[va])
    tx_tr, tcl_tr, tyt_tr, tye_tr = prep(device, dtype, X_tr, cl_tr, yt[tr], ye[tr])
    tx_va, tcl_va, tyt_va, tye_va = prep(device, dtype, X_va, cl_va, yt[va], ye[va])

    total = len(GRID_ARCHITECTURES) * len(GRID_L2) * len(GRID_LR)
    n_clin = len(CLIN_FEATURES)
    print(f'Grid search: {total} combos  |  '
          f'n_genes={n_genes}  |  n_clin={n_clin} {CLIN_FEATURES}')
    print(f'Split: train={len(tr)}, val={len(va)}\n')

    results, best, idx = [], None, 0
    for h1, h2, out in GRID_ARCHITECTURES:
        for l2 in GRID_L2:
            for lr in GRID_LR:
                idx += 1
                torch.manual_seed(GS_SEED)
                np.random.seed(GS_SEED)

                net = CoxPhRWRNet(
                    input_nodes=n_genes, hidden_nodes1=h1,
                    hidden_nodes2=h2, output_nodes=out,
                    n_clin=n_clin,
                ).to(device)
                opt = optim.Adam(net.parameters(), lr=lr, weight_decay=l2)
                sch = optim.lr_scheduler.CosineAnnealingLR(
                    opt, T_max=GS_EPOCHS, eta_min=1e-6)
                best_val, best_st = 0.0, copy.deepcopy(net.state_dict())

                for ep in range(GS_EPOCHS):
                    net.train(); opt.zero_grad()
                    loss = cox_loss(net(tx_tr, tcl_tr), tyt_tr, tye_tr)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(net.parameters(), GS_GRAD_CLIP)
                    opt.step(); sch.step()
                    if ep % GS_CHECKPOINT == 0:
                        net.eval()
                        with torch.no_grad():
                            vc = c_index(net(tx_va, tcl_va), tyt_va, tye_va).item()
                        if vc > best_val:
                            best_val, best_st = vc, copy.deepcopy(net.state_dict())

                row = dict(hidden1=h1, hidden2=h2, output_nodes=out,
                           L2=l2, LR=lr, val_c=best_val)
                results.append(row)
                print(f'[{idx:02d}/{total}] arch=({h1},{h2},{out}) '
                      f'L2={l2:.3f} LR={lr:.0e} | val_c={best_val:.4f}')
                if best is None or best_val > best['val_c']:
                    best = row.copy()

    print(f'\nBest: arch=({best["hidden1"]},{best["hidden2"]},{best["output_nodes"]}) '
          f'L2={best["L2"]:.3f} LR={best["LR"]:.0e} | val_c={best["val_c"]:.4f}')
    gs_df = pd.DataFrame(results).sort_values('val_c', ascending=False)
    return best, gs_df

def run_alpha_k_grid_search(
    device: torch.device,
    dtype,
    clinical,
    W,
    nodes: list,
    node_to_idx: dict,
    seed_by_patient: dict,
    cosmic_in_net: list,
    cosmic_idx: np.ndarray,
    base_model_kwargs: dict,
    lr: float,
    l2: float,
    k_values: list = None,
):
    """
    Grid search over RWR restart probability (alpha) and number of genes (K).

    For each (alpha, K) pair:
      1. Re-compute RWR with the given alpha.
      2. Select the top-K genes by log-rank p-value.
      3. Train and evaluate with N_REPEATS random splits.

    Args:
    clinical         : pd.DataFrame   Patient table with SAMPLE_ID column.
    W                : sparse matrix  PPI transition matrix.
    nodes            : list[str]      All network gene names.
    node_to_idx      : dict           maps gene name to its index in W.
    seed_by_patient  : dict           maps patient ID to their list of mutated seed genes.
    cosmic_in_net    : list[str]      COSMIC genes present in the network.
    cosmic_idx       : np.ndarray     Their indices in W.
    base_model_kwargs: dict           Architecture kwargs (input_nodes updated per K).
    lr, l2           : float          Learning rate and weight decay from grid search.

    Returns:
    grid_results : dict  keys are (alpha, K) tuples, values have 'test' and 'val' lists
    best_alpha   : float
    best_K       : int
    k_values     : list[int]
    """
    n_cosmic_total  = len(cosmic_in_net)
    if k_values is None:
        k_values = list(range(100, n_cosmic_total, 100))
        if not k_values or k_values[-1] != n_cosmic_total:
            k_values.append(n_cosmic_total)
    else:
        k_values = [k for k in k_values if k <= n_cosmic_total]
        if not k_values or k_values[-1] != n_cosmic_total:
            k_values.append(n_cosmic_total)

    base_sample_ids = set(clinical['SAMPLE_ID'])
    print(f'Grid: alpha={ALPHA_VALUES}')
    print(f'Grid: K={k_values}  (n_cosmic_in_net={n_cosmic_total})')
    print(f'Total: {len(ALPHA_VALUES) * len(k_values)} combos × {N_REPEATS} repeats\n')

    grid_results = {}
    for alpha in ALPHA_VALUES:
        # Re-compute RWR once per alpha (reused across all K values)
        rwr_a = compute_rwr_scores(
            clinical, W, nodes, node_to_idx,
            seed_by_patient, cosmic_idx, alpha=alpha)
        df_a  = make_feature_df(rwr_a, cosmic_in_net, clinical)
        df_a  = df_a[df_a['SAMPLE_ID'].isin(base_sample_ids)].reset_index(drop=True)

        gene_arr_a = df_a[cosmic_in_net].values.astype(np.float32)
        yt_a  = df_a['OS_MONTHS'].values.astype(np.float32)
        ye_a  = df_a['OS_STATUS'].values.astype(np.float32)
        cl_a  = df_a[CLIN_FEATURES].values.astype(np.float32)  # AGE, msi_score, tmb_score

        # Log-rank p-values computed once per alpha, reused across K
        pvals_a    = logrank_pvalues(gene_arr_a, yt_a, ye_a)
        sorted_idx = np.argsort(pvals_a)

        for k in k_values:
            X_a = (gene_arr_a[:, sorted_idx[:k]]
                   if k < n_cosmic_total else gene_arr_a)
            mk = {**base_model_kwargs, 'input_nodes': X_a.shape[1]}

            test_cis, val_cis = train_and_eval(
                device, dtype, X_a, cl_a, yt_a, ye_a, mk, lr=lr, l2=l2)
            grid_results[(alpha, k)] = {'test': test_cis, 'val': val_cis}

            print(f'alpha={alpha}  K={k:4d}  '
                  f'val={np.mean(val_cis):.4f}±{np.std(val_cis):.4f}  '
                  f'test={np.mean(test_cis):.4f}±{np.std(test_cis):.4f}')

    best_alpha, best_K = max(
        grid_results, key=lambda ak: np.mean(grid_results[ak]['val']))
    print(f'\nBest: alpha={best_alpha}  K={best_K}  '
          f'val={np.mean(grid_results[(best_alpha, best_K)]["val"]):.4f}  '
          f'test={np.mean(grid_results[(best_alpha, best_K)]["test"]):.4f}')
    return grid_results, best_alpha, best_K, k_values
