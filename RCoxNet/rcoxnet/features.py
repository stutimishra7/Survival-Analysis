import numpy as np
import pandas as pd
from lifelines.statistics import logrank_test

from rcoxnet.config import META_COLS


def make_feature_df(rwr_mat, feat_names, clinical):
    """
    Merge the RWR score matrix with clinical metadata.
    Drops patients with missing survival info or no RWR signal.
    """
    rwr_df = pd.DataFrame(rwr_mat, columns=feat_names)
    rwr_df.insert(0, 'SAMPLE_ID', clinical['SAMPLE_ID'].values)

    merged = rwr_df.merge(clinical[META_COLS], on='SAMPLE_ID', how='inner')
    merged = (merged
              .dropna(subset=['OS_MONTHS', 'OS_STATUS'])
              .reset_index(drop=True))
    merged['OS_STATUS'] = merged['OS_STATUS'].astype(int)

    has_signal = (merged[feat_names].values > 0).any(axis=1)
    return merged[has_signal].reset_index(drop=True)


def logrank_pvalues(gene_arr, yt, ye):
    """
    Log-rank p-value for each gene using a median split of RWR scores.
    Genes where either group has fewer than 5 patients get p=1.
    """
    pvals = []
    for j in range(gene_arr.shape[1]):
        col  = gene_arr[:, j]
        mask = col > np.median(col)
        if mask.sum() < 5 or (~mask).sum() < 5:
            pvals.append(1.0)
            continue
        r = logrank_test(yt[mask], yt[~mask], ye[mask], ye[~mask])
        pvals.append(r.p_value)
    return np.array(pvals)


def select_top_k(gene_arr, gene_names, yt, ye, k):
    """Select the top-K genes ranked by log-rank p-value."""
    pvals      = logrank_pvalues(gene_arr, yt, ye)
    sorted_idx = np.argsort(pvals)
    top_idx    = sorted_idx[:k] if k < len(gene_names) else sorted_idx
    return (gene_arr[:, top_idx],
            [gene_names[i] for i in top_idx],
            pvals,
            sorted_idx)
