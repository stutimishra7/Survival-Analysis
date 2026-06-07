"""
Data normalisation and train/val/test splitting.

Splitting is stratified by event status and survival-time tertile so that
each fold has a similar event rate and time distribution.
"""

import numpy as np
import torch


def normalise_genomic(X_tr, *others):
    """
    Log10 + z-score normalisation for the RWR feature matrix.
    Mean and standard deviation are computed on the training split only
    and then applied to the validation and test splits.
    """
    X_tr = np.log10(X_tr + 1e-12)
    mu   = X_tr.mean(axis=0, keepdims=True)
    sig  = X_tr.std(axis=0, keepdims=True)
    sig[sig < 1e-8] = 1.0
    X_tr = (X_tr - mu) / sig
    return [X_tr] + [(np.log10(X + 1e-12) - mu) / sig for X in others]


def normalise_clinical(cl_tr, *others):
    """
    Z-score normalisation for clinical features (AGE, msi_score, tmb_score).
    Missing values are imputed with the training-set mean before scaling.
    Mean and std are fit on the training split and applied to the others.
    """
    mu  = np.nanmean(cl_tr, axis=0, keepdims=True)
    sig = np.nanstd(cl_tr, axis=0, keepdims=True)
    sig[sig < 1e-8] = 1.0
    cl_tr = np.where(np.isnan(cl_tr), mu, cl_tr)
    return [(cl_tr - mu) / sig] + [
        ((np.where(np.isnan(c), mu, c)) - mu) / sig for c in others
    ]


def sort_by_time_desc(yt, *arrs):
    # Cox loss requires descending survival time order
    idx = np.argsort(yt)[::-1].copy()
    return [a[idx] for a in (yt,) + arrs]


def to_tensor(device, dtype, *arrays):
    return [
        torch.from_numpy(a.reshape(-1, 1) if a.ndim == 1 else a)
        .type(dtype).to(device)
        for a in arrays
    ]


def prep(device, dtype, X_norm, cl_norm, yt, ye):
    # sort by descending time, then convert to tensors
    yt_s, X_s, cl_s, ye_s = sort_by_time_desc(yt, X_norm, cl_norm, ye)
    return to_tensor(device, dtype,
                     X_s, cl_s,
                     yt_s.astype(np.float32),
                     ye_s.astype(np.float32))


def censoring_time_balanced_split(yt, ye, test_frac, val_frac, rng):
    """
    Train/val/test split that keeps event rate stable across folds.

    Event patients are split into tertiles by survival time; censored patients
    form a fourth group. Each group is split proportionally into train/val/test.
    """
    yt = np.asarray(yt)
    ye = np.asarray(ye).astype(int)

    ev_idx = np.where(ye == 1)[0]
    if len(ev_idx) >= 8:
        q1, q2  = np.quantile(yt[ev_idx], [0.33, 0.66])
        ev_bins = np.digitize(yt[ev_idx], bins=[q1, q2], right=True)
    else:
        med     = np.median(yt[ev_idx]) if len(ev_idx) > 0 else np.median(yt)
        ev_bins = (yt[ev_idx] > med).astype(int)

    # three event-time tertiles plus one group for censored patients
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

    tr = np.concatenate([x for x in tr_parts if len(x) > 0])
    va = np.concatenate([x for x in va_parts if len(x) > 0])
    te = np.concatenate([x for x in te_parts if len(x) > 0])
    rng.shuffle(tr); rng.shuffle(va); rng.shuffle(te)
    return tr, va, te
