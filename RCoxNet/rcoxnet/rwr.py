"""
Random Walk with Restart (RWR) over the PPI network.

Each patient's somatic mutations seed the walk; the resulting
stationary distribution is used as the genomic feature vector.
"""

import numpy as np
import pandas as pd

from rcoxnet.config import RWR_MAX, RWR_TOL


def run_rwr(W, seed_vec, alpha=0.3):
    """Single RWR run — iterates until convergence or RWR_MAX steps."""
    p = seed_vec.copy()
    for _ in range(RWR_MAX):
        p_new = (1.0 - alpha) * W.dot(p) + alpha * seed_vec
        if np.linalg.norm(p_new - p, ord=1) < RWR_TOL:
            break
        p = p_new
    return p


def compute_rwr_scores(clinical, W, nodes, node_to_idx,
                       seed_by_patient, feat_idx, alpha=0.3,
                       weighted=False, mut_count_by_patient=None):
    """
    Compute the RWR feature matrix for all patients.

    Returns an (n_patients, len(feat_idx)) array of RWR scores.
    Patients with no mutated seeds in the network get an all-zero row.
    """
    n_pts  = len(clinical)
    n_fea  = len(feat_idx)
    scores = np.zeros((n_pts, n_fea), dtype=np.float32)

    for i, row in clinical.iterrows():
        pid   = row['PATIENT_ID']
        seeds = seed_by_patient.get(pid, [])
        valid = [g for g in seeds if g in node_to_idx]
        if not valid:
            continue

        s = np.zeros(len(nodes), dtype=np.float64)
        for g in valid:
            if weighted and mut_count_by_patient:
                cnt = mut_count_by_patient.get(pid, {}).get(g, 1)
                s[node_to_idx[g]] = float(cnt)
            else:
                s[node_to_idx[g]] = 1.0
        s /= s.sum()

        scores[i] = run_rwr(W, s, alpha=alpha)[feat_idx].astype(np.float32)

    return scores
