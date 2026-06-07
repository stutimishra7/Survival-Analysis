"""
Builds the column-stochastic PPI transition matrix from a ConsensusPathDB
tissue-specific interaction file.

File format (tab-separated, 2 header lines):
    gene1  gene2  IntConf  tpm_1  tpm_2
"""

import numpy as np
import pandas as pd
from scipy import sparse


def build_transition_matrix(txt_path, int_conf_thresh=0, tpm_thresh=0):
    edges = pd.read_csv(
        txt_path, sep='\t', skiprows=2, header=None,
        names=['src', 'dst', 'IntConf', 'tpm_1', 'tpm_2'],
    )
    edges = edges.dropna()
    for col in ['IntConf', 'tpm_1', 'tpm_2']:
        edges[col] = pd.to_numeric(edges[col], errors='coerce')
    edges = edges.dropna(subset=['IntConf', 'tpm_1', 'tpm_2'])

    before = len(edges)
    edges  = edges[
        (edges['IntConf'] > int_conf_thresh) &
        (edges['tpm_1']   >= tpm_thresh) &
        (edges['tpm_2']   >= tpm_thresh)
    ].copy()
    print(f'  Edges before filter: {before:,}  after: {len(edges):,}  '
          f'(IntConf>{int_conf_thresh}, tpm>={tpm_thresh})')

    edges = edges[edges['src'] != edges['dst']].copy()
    edges['weight'] = edges['IntConf']

    nodes       = sorted(set(edges['src']) | set(edges['dst']))
    node_to_idx = {g: i for i, g in enumerate(nodes)}
    n           = len(nodes)

    si = edges['src'].map(node_to_idx).values
    di = edges['dst'].map(node_to_idx).values
    w  = edges['weight'].values.astype(np.float64)

    # symmetrise (undirected graph)
    A = sparse.coo_matrix(
        (np.concatenate([w, w]),
         (np.concatenate([si, di]), np.concatenate([di, si]))),
        shape=(n, n),
    ).tocsr()

    # column-normalise to get a column-stochastic transition matrix
    col_sums = np.asarray(A.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    W = A @ sparse.diags(1.0 / col_sums)

    edge_df = edges[['src', 'dst', 'weight']].reset_index(drop=True)
    print(f'  Network genes: {n:,}  |  Non-zero entries: {W.nnz:,}')
    return nodes, node_to_idx, W, edge_df
