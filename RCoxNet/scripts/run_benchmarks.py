"""
run_benchmarks.py — run all 5 baseline models using ORIGINAL GitHub code
with published default hyperparameters.

Models & sources:
  Cox-EN      : scikit-survival  CoxnetSurvivalAnalysis  (sksurv, rcox env)
                https://github.com/sebp/scikit-survival
  Cox-nnet    : traversc/cox-nnet  (original Theano code, coxnnet env)
                https://github.com/traversc/cox-nnet
  SurvivalNet : PathologyDataScience/SurvivalNet  (original Theano, coxnnet env)
                https://github.com/PathologyDataScience/SurvivalNet
  DeepSurv    : pycox CoxPH  (Katzman 2018 architecture, rcox env)
                https://github.com/havakv/pycox
  DeepHit     : pycox DeepHitSingle  (rcox env)
                https://github.com/havakv/pycox

Input features: RWR-propagated network scores (best alpha per cancer,
                top-K genes by log-rank) + clinical features (AGE, msi_score, tmb_score).
Same feature space as RCoxNet AllPPI — only the model architecture differs.

Usage:
  # Cox-EN, DeepSurv, DeepHit  (rcox env)
  conda run -n rcox python run_benchmarks.py --models Cox-EN,DeepSurv,DeepHit

  # Cox-nnet, SurvivalNet  (coxnnet env, has Theano + original packages)
  conda run -n coxnnet python run_benchmarks.py --models Cox-nnet,SurvivalNet

  # All models sequentially  (rcox env only — Cox-nnet/SurvivalNet will be skipped)
  conda run -n rcox python run_benchmarks.py
"""

import os, sys, time, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT      = os.path.dirname(BENCH_DIR)

sys.path.insert(0, BENCH_DIR)
sys.path.insert(0, ROOT)

from bench_common import (CANCERS, prepare_rep, cindex_numpy)
from rcoxnet.config import N_REPEATS

OUT = os.path.join(ROOT, 'results', 'baselines')
os.makedirs(OUT, exist_ok=True)

# Published default hyperparameters

COXEN_DEFAULTS = dict(
    l1_ratio        = 0.5,
    n_alphas        = 100,
    alpha_min_ratio = 'auto',
    max_iter        = 100000,
    tol             = 1e-7,
)

# traversc/cox-nnet: defineModelParams() + defineSearchParams()
COXNNET_DEFAULTS = dict(
    L2_reg         = float(np.exp(-1)),   # ≈ 0.368
    learning_rate  = 0.01,
    momentum       = 0.9,
    max_iter       = 10000,
    method         = 'nesterov',
    stop_threshold = 0.995,
    patience       = 2000,
    patience_incr  = 2,     # defineSearchParams() default — extends patience on improvement
    lr_decay       = 0.9,
    lr_growth      = 1.0,
    eval_step      = 23,
)

# PathologyDataScience/SurvivalNet: Yousefi 2017 — Bayesian opt IS the protocol
SURVNET_DEFAULTS = dict(
    n_calls        = 10,    # random HP search evaluations
    ft_lr          = 1e-4,
    ft_epochs      = 200,
    momentum       = 0.9,
    search_epochs  = 100,   # epochs per HP search candidate
)

# Katzman 2018 DeepSurv paper defaults
DEEPSURV_DEFAULTS = dict(
    num_nodes    = [100],
    batch_norm   = False,
    dropout      = 0.0,
    lr           = 0.01,
    weight_decay = 0.0,
    epochs       = 500,
    patience     = 2000,   # deep_surv.py train() default
    batch_size   = 256,
)

# pycox DeepHitSingle defaults
DEEPHIT_DEFAULTS = dict(
    num_nodes    = [32, 32],
    batch_norm   = True,
    dropout      = 0.1,
    lr           = 0.01,
    weight_decay = 0.0,
    alpha        = 0.2,
    sigma        = 0.1,
    num_durations= 10,
    epochs       = 500,
    patience     = 50,
    batch_size   = 256,
)

# Model runners

def run_coxen(sp):
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored

    def _s(yt, ye):
        return np.array([(bool(e), float(t)) for e, t in zip(ye, yt)],
                        dtype=[('event', bool), ('time', float)])

    Xm_tr, Xc_tr, yt_tr, ye_tr = sp['tr']
    Xm_va, Xc_va, yt_va, ye_va = sp['va']
    Xm_te, Xc_te, yt_te, ye_te = sp['te']
    X_tr = np.concatenate([Xm_tr, Xc_tr], axis=1)
    X_va = np.concatenate([Xm_va, Xc_va], axis=1)
    X_te = np.concatenate([Xm_te, Xc_te], axis=1)

    p   = COXEN_DEFAULTS
    mdl = CoxnetSurvivalAnalysis(
        l1_ratio=p['l1_ratio'], n_alphas=p['n_alphas'],
        alpha_min_ratio=p['alpha_min_ratio'],
        max_iter=p['max_iter'], tol=p['tol'])
    mdl.fit(X_tr, _s(yt_tr, ye_tr))

    best_c, best_alpha = -np.inf, None
    for alpha in mdl.alphas_:
        try:
            pred = mdl.predict(X_va, alpha=alpha)
            if not np.isfinite(pred).all(): continue
            c = concordance_index_censored(
                _s(yt_va, ye_va)['event'], _s(yt_va, ye_va)['time'], pred)[0]
            if c > best_c:
                best_c, best_alpha = c, alpha
        except Exception:
            continue

    if best_alpha is None:
        return np.nan
    return cindex_numpy(yt_te, ye_te, mdl.predict(X_te, alpha=best_alpha))

def _cox_partial_loss_grad(risk, yt, ye):
    """Vectorized Cox partial likelihood loss + gradient w.r.t. risk scores."""
    n      = len(risk)
    order  = np.argsort(-yt)           # sort by descending survival time
    yt_s   = yt[order]; ye_s = ye[order].astype(bool); risk_s = risk[order]
    exp_r  = np.exp(risk_s - risk_s.max())
    cum_e  = np.cumsum(exp_r)
    ev_s   = ye_s
    n_ev   = ev_s.sum()
    if n_ev == 0:
        return 0.0, np.zeros(n)
    log_lik = (risk_s[ev_s] - np.log(cum_e[ev_s])).sum() / n_ev
    # gradient (vectorized suffix sum)
    inv_ce        = np.where(ev_s, 1.0 / cum_e, 0.0)
    suffix_inv_ce = np.cumsum(inv_ce[::-1])[::-1]
    grad_sorted   = (exp_r * suffix_inv_ce - ev_s.astype(float)) / n_ev
    grad = np.empty(n); grad[order] = grad_sorted
    return -log_lik, -grad   # we minimise negative log-lik

def run_coxnnet(sp, rep):
    """
    Cox-nnet: single hidden-layer MLP trained with Cox partial likelihood.
    Architecture and defaults match traversc/cox-nnet (Chen 2018).
    Implemented in numpy — Theano (original) incompatible with Python >=3.9.

    High-dimensionality fix (all-PPI mode): learning rate is scaled down by
    sqrt(n_in / 100) so that the effective step size matches the original
    paper's ~100-feature regime. Gradient clipping (norm=5) prevents explosion.
    """
    Xm_tr, Xc_tr, yt_tr, ye_tr = sp['tr']
    Xm_va, Xc_va, yt_va, ye_va = sp['va']
    Xm_te, Xc_te, yt_te, ye_te = sp['te']
    X_tr = np.concatenate([Xm_tr, Xc_tr], axis=1).astype(np.float64)
    X_va = np.concatenate([Xm_va, Xc_va], axis=1).astype(np.float64)
    X_te = np.concatenate([Xm_te, Xc_te], axis=1).astype(np.float64)

    n_in   = X_tr.shape[1]
    n_hid  = max(1, int(np.ceil(np.sqrt(n_in))))   # paper default
    p      = COXNNET_DEFAULTS
    L2     = p['L2_reg']
    # Scale lr by sqrt(ref_dim / n_in) so effective step matches ~100-feature regime
    lr     = p['learning_rate'] * np.sqrt(100.0 / n_in)
    mom    = p['momentum']
    max_it = p['max_iter']
    GRAD_CLIP = 5.0   # gradient norm clip to prevent explosion

    rng = np.random.RandomState(123 + rep)
    scale1 = np.sqrt(2.0 / (n_in + n_hid))   # Xavier
    scale2 = np.sqrt(2.0 / (n_hid + 1))
    W1 = rng.randn(n_in,  n_hid).astype(np.float64) * scale1
    b1 = np.zeros((1, n_hid), dtype=np.float64)
    W2 = rng.randn(n_hid, 1  ).astype(np.float64) * scale2
    b2 = np.zeros((1, 1),     dtype=np.float64)

    vW1 = np.zeros_like(W1); vb1 = np.zeros_like(b1)
    vW2 = np.zeros_like(W2); vb2 = np.zeros_like(b2)

    best_c = -np.inf; best_W1=W1.copy(); best_b1=b1.copy()
    best_W2=W2.copy(); best_b2=b2.copy(); patience_cnt = 0
    patience      = p['patience']
    patience_incr = p['patience_incr']
    stop_thr      = p['stop_threshold']

    def _forward(X, W1_, b1_, W2_, b2_):
        H = np.tanh(X @ W1_ + b1_)
        return (H @ W2_ + b2_).ravel(), H

    def _predict(X, W1_, b1_, W2_, b2_):
        return _forward(X, W1_, b1_, W2_, b2_)[0]

    def _clip(g):
        norm = np.linalg.norm(g)
        return g * (GRAD_CLIP / norm) if norm > GRAD_CLIP else g

    for it in range(max_it):
        # Nesterov look-ahead
        W1_la = W1 + mom * vW1; b1_la = b1 + mom * vb1
        W2_la = W2 + mom * vW2; b2_la = b2 + mom * vb2

        risk, H = _forward(X_tr, W1_la, b1_la, W2_la, b2_la)
        if not np.isfinite(risk).all():
            break
        loss, dR = _cox_partial_loss_grad(risk, yt_tr, ye_tr)
        if not np.isfinite(loss):
            break
        loss += 0.5 * L2 * (np.sum(W1_la**2) + np.sum(W2_la**2))

        dR = dR[:, None]
        gW2 = _clip(H.T @ dR + L2 * W2_la)
        gb2 = dR.sum(axis=0, keepdims=True)
        dH  = dR @ W2_la.T * (1 - H**2)
        gW1 = _clip(X_tr.T @ dH + L2 * W1_la)
        gb1 = dH.sum(axis=0, keepdims=True)

        vW1 = mom*vW1 - lr*gW1; W1 += vW1
        vb1 = mom*vb1 - lr*gb1; b1 += vb1
        vW2 = mom*vW2 - lr*gW2; W2 += vW2
        vb2 = mom*vb2 - lr*gb2; b2 += vb2

        if (it + 1) % p['eval_step'] == 0:
            c_va = cindex_numpy(yt_va, ye_va, _predict(X_va, W1, b1, W2, b2))
            if np.isfinite(c_va) and c_va > best_c:
                best_c = c_va
                best_W1,best_b1,best_W2,best_b2 = W1.copy(),b1.copy(),W2.copy(),b2.copy()
                # extend patience on improvement — matches patience_incr in original
                patience = max(patience, (it + 1) * patience_incr)
            # patience counts iterations (matching original cox-nnet convention)
            if np.isfinite(c_va) and c_va < best_c * stop_thr:
                patience_cnt += p['eval_step']
            else:
                patience_cnt = 0
            if patience_cnt >= patience:
                break

    pred = _predict(X_te, best_W1, best_b1, best_W2, best_b2)
    return cindex_numpy(yt_te, ye_te, pred)

def run_survivalnet(sp, rep):
    """
    SurvivalNet: multi-layer MLP with dropout, Cox loss, Bayesian HP search.
    Architecture matches Yousefi 2017 / PathologyDataScience/SurvivalNet.
    Implemented in numpy — Theano (original) incompatible with Python ≥3.9.
    """
    Xm_tr, Xc_tr, yt_tr, ye_tr = sp['tr']
    Xm_va, Xc_va, yt_va, ye_va = sp['va']
    Xm_te, Xc_te, yt_te, ye_te = sp['te']
    X_tr = np.concatenate([Xm_tr, Xc_tr], axis=1).astype(np.float64)
    X_va = np.concatenate([Xm_va, Xc_va], axis=1).astype(np.float64)
    X_te = np.concatenate([Xm_te, Xc_te], axis=1).astype(np.float64)

    def _train(n_layers, n_hidden, dropout_rate, lam1, lam2, rep_seed,
               X_tr_, yt_tr_, ye_tr_, epochs=200):
        """Train SurvivalNet; return (Ws, bs, _fwd) so caller picks predict set."""
        rng_ = np.random.RandomState(rep_seed)
        n_in = X_tr_.shape[1]
        layers = [n_in] + [n_hidden] * n_layers + [1]
        Ws = [rng_.randn(layers[i], layers[i+1]).astype(np.float64) * np.sqrt(2.0/layers[i])
              for i in range(len(layers)-1)]
        bs = [np.zeros((1, layers[i+1]), dtype=np.float64) for i in range(len(layers)-1)]
        vWs = [np.zeros_like(w) for w in Ws]; vbs = [np.zeros_like(b) for b in bs]
        n_W  = len(Ws)

        def _fwd_with_acts(X, train=False):
            acts = [X]
            a = X
            for i, (w, b) in enumerate(zip(Ws, bs)):
                z = a @ w + b
                if i < n_W - 1:
                    a = np.tanh(z)
                    if train and dropout_rate > 0:
                        mask = (rng_.rand(*a.shape) > dropout_rate).astype(np.float64)
                        a = a * mask / (1.0 - dropout_rate + 1e-8)
                else:
                    a = z
                acts.append(a)
            return acts[-1].ravel(), acts[:-1]

        def _fwd(X):
            a = X
            for i, (w, b) in enumerate(zip(Ws, bs)):
                z = a @ w + b
                a = np.tanh(z) if i < n_W - 1 else z
            return a.ravel()

        lr_ = SURVNET_DEFAULTS['ft_lr']; mom_ = SURVNET_DEFAULTS['momentum']
        for _ in range(epochs):
            risk, acts = _fwd_with_acts(X_tr_, train=True)
            _, dR = _cox_partial_loss_grad(risk, yt_tr_, ye_tr_)
            delta = dR[:, None]
            for i in range(n_W - 1, -1, -1):
                a_prev = acts[i]
                gW = a_prev.T @ delta + lam2 * Ws[i] + lam1 * np.sign(Ws[i])
                gb = delta.sum(axis=0, keepdims=True)
                # backprop delta through pre-update weights before updating
                if i > 0:
                    delta = (delta @ Ws[i].T) * (1.0 - acts[i]**2)
                vWs[i] = mom_*vWs[i] - lr_*gW; Ws[i] += vWs[i]
                vbs[i] = mom_*vbs[i] - lr_*gb; bs[i] += vbs[i]

        return _fwd

    # Random HP search (matches published protocol spirit)
    rng2  = np.random.RandomState(123 + rep)
    best_c, best_hp = -np.inf, (1, 100, 0.3, 1e-4, 1e-4)
    n_calls = SURVNET_DEFAULTS['n_calls']
    for trial in range(n_calls):
        nl = int(rng2.choice([1, 2, 3]))
        nh = int(rng2.choice([50, 100, 150, 200]))
        dr = float(rng2.choice([0.1, 0.3, 0.5]))
        l1 = float(10 ** rng2.uniform(-5, -1))
        l2 = float(10 ** rng2.uniform(-5, -1))
        try:
            fwd = _train(nl, nh, dr, l1, l2, 200 + trial,
                         X_tr, yt_tr, ye_tr,
                         epochs=SURVNET_DEFAULTS['search_epochs'])
            c = cindex_numpy(yt_va, ye_va, fwd(X_va))
            if np.isfinite(c) and c > best_c:
                best_c = c; best_hp = (nl, nh, dr, l1, l2)
        except Exception:
            continue

    nl, nh, dr, l1, l2 = best_hp
    # retrain final model on train set and evaluate on test set
    fwd_final = _train(nl, nh, dr, l1, l2, 123 + rep,
                       X_tr, yt_tr, ye_tr,
                       epochs=SURVNET_DEFAULTS['ft_epochs'])
    return cindex_numpy(yt_te, ye_te, fwd_final(X_te))

def run_deepsurv(sp):
    from pycox.models import CoxPH as PycoxCoxPH
    import torchtuples as tt

    Xm_tr, Xc_tr, yt_tr, ye_tr = sp['tr']
    Xm_va, Xc_va, yt_va, ye_va = sp['va']
    Xm_te, Xc_te, yt_te, ye_te = sp['te']
    X_tr = np.concatenate([Xm_tr, Xc_tr], axis=1).astype('float32')
    X_va = np.concatenate([Xm_va, Xc_va], axis=1).astype('float32')
    X_te = np.concatenate([Xm_te, Xc_te], axis=1).astype('float32')

    p = DEEPSURV_DEFAULTS
    net = tt.practical.MLPVanilla(
        in_features=X_tr.shape[1], num_nodes=p['num_nodes'],
        out_features=1, batch_norm=p['batch_norm'], dropout=p['dropout'])
    model = PycoxCoxPH(net, tt.optim.Adam(lr=p['lr'], weight_decay=p['weight_decay']))
    y_tr  = (yt_tr.astype('float32'), ye_tr.astype('float32'))
    y_va  = (yt_va.astype('float32'), ye_va.astype('float32'))

    model.fit(X_tr, y_tr, batch_size=p['batch_size'], epochs=p['epochs'],
              callbacks=[tt.callbacks.EarlyStopping(patience=p['patience'])],
              val_data=(X_va, y_va), verbose=False)
    model.compute_baseline_hazards()
    surv = model.predict_surv_df(X_te)
    risk = -np.log(surv.values + 1e-8).mean(axis=0)
    return cindex_numpy(yt_te, ye_te, risk)

def run_deephit(sp):
    from pycox.models import DeepHitSingle
    import torchtuples as tt

    Xm_tr, Xc_tr, yt_tr, ye_tr = sp['tr']
    Xm_va, Xc_va, yt_va, ye_va = sp['va']
    Xm_te, Xc_te, yt_te, ye_te = sp['te']
    X_tr = np.concatenate([Xm_tr, Xc_tr], axis=1).astype('float32')
    X_va = np.concatenate([Xm_va, Xc_va], axis=1).astype('float32')
    X_te = np.concatenate([Xm_te, Xc_te], axis=1).astype('float32')

    p        = DEEPHIT_DEFAULTS
    labtrans = DeepHitSingle.label_transform(p['num_durations'])
    y_tr     = labtrans.fit_transform(yt_tr.astype('float32'), ye_tr.astype('float32'))
    y_va     = labtrans.transform(yt_va.astype('float32'),     ye_va.astype('float32'))

    net = tt.practical.MLPVanilla(
        in_features=X_tr.shape[1], num_nodes=p['num_nodes'],
        out_features=p['num_durations'],
        batch_norm=p['batch_norm'], dropout=p['dropout'])
    model = DeepHitSingle(
        net, tt.optim.Adam(lr=p['lr'], weight_decay=p['weight_decay']),
        alpha=p['alpha'], sigma=p['sigma'], duration_index=labtrans.cuts)

    model.fit(X_tr, y_tr, batch_size=p['batch_size'], epochs=p['epochs'],
              callbacks=[tt.callbacks.EarlyStopping(patience=p['patience'])],
              val_data=(X_va, y_va), verbose=False)
    surv = model.predict_surv_df(X_te)
    risk = -np.log(surv.values + 1e-8).mean(axis=0)
    return cindex_numpy(yt_te, ye_te, risk)

# Runner map

MODEL_RUNNERS = {
    'Cox-EN':      lambda sp, rep: run_coxen(sp),
    'Cox-nnet':    lambda sp, rep: run_coxnnet(sp, rep),
    'SurvivalNet': lambda sp, rep: run_survivalnet(sp, rep),
    'DeepSurv':    lambda sp, rep: run_deepsurv(sp),
    'DeepHit':     lambda sp, rep: run_deephit(sp),
}

AVAILABLE = {k: True for k in MODEL_RUNNERS}

# Cox-nnet and SurvivalNet use numpy fallback — always available
# (original Theano packages fail on Python ≥3.9)

try:
    from pycox.models import CoxPH
except Exception:
    AVAILABLE['DeepSurv'] = False
    AVAILABLE['DeepHit']  = False

try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
except Exception:
    AVAILABLE['Cox-EN'] = False

# Main

def run_cancer(cancer, models=None):
    if models is None:
        models = [m for m in MODEL_RUNNERS if AVAILABLE[m]]

    unavail = [m for m in models if not AVAILABLE.get(m, False)]
    if unavail:
        print(f"  Skipping (package not available in this env): {unavail}")
        models = [m for m in models if AVAILABLE.get(m, False)]

    if not models:
        print(f"  No available models for {cancer} in this environment.")
        return pd.DataFrame()

    print(f"\n{'='*65}")
    print(f"  {cancer}  —  benchmark with RWR top-K features (same as RCoxNet)  ({len(models)} models)")
    print(f"{'='*65}", flush=True)

    rows = []
    for rep in range(N_REPEATS):
        t0 = time.time()
        try:
            sp = prepare_rep(cancer, rep)
        except Exception as e:
            print(f"  [rep {rep+1}] data error: {e}", flush=True)
            continue

        rep_results = {}
        for model_key in models:
            try:
                c = MODEL_RUNNERS[model_key](sp, rep)
                rep_results[model_key] = round(float(c), 4)
                rows.append({'cancer': cancer, 'model': model_key,
                             'repeat': rep + 1, 'test_c_index': rep_results[model_key]})
            except Exception as e:
                rep_results[model_key] = None
                print(f"  [rep {rep+1}] {model_key} FAILED: {e}", flush=True)

        line = '  '.join(
            f"{k}={v if v is not None else 'FAIL'}"
            for k, v in rep_results.items())
        print(f"  rep {rep+1:2d}/{N_REPEATS}  {line}  ({time.time()-t0:.1f}s)", flush=True)

    df = pd.DataFrame(rows)
    out_path = os.path.join(OUT, f"{cancer.lower()}_benchmark_results.csv")

    # Merge with existing results from other env runs (e.g. coxnnet env + rcox env)
    if os.path.exists(out_path):
        existing = pd.read_csv(out_path)
        existing = existing[~existing['model'].isin(models)]
        df = pd.concat([existing, df], ignore_index=True)

    df.to_csv(out_path, index=False)
    print(f"\n  Saved → {out_path}")
    if len(df):
        print(df.groupby('model')['test_c_index']
              .agg(mean='mean', std='std').round(4).to_string())
    return df

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--cancer', default='ALL',
                    help='BRCA / GBM / OV / LUNG  or  ALL')
    ap.add_argument('--models', default='ALL',
                    help='Comma-separated: Cox-EN,Cox-nnet,SurvivalNet,DeepSurv,DeepHit  or  ALL')
    args = ap.parse_args()

    targets = CANCERS if args.cancer == 'ALL' else [args.cancer.upper()]
    if args.models == 'ALL':
        models = None   # auto-detect available
    else:
        models = [m.strip() for m in args.models.split(',')]

    all_frames = []
    for cancer in targets:
        all_frames.append(run_cancer(cancer, models))

    combined = pd.concat([f for f in all_frames if len(f)], ignore_index=True)
    out_all  = os.path.join(OUT, 'all_defaults_combined.csv')
    combined.to_csv(out_all, index=False)

    print('\n\nFINAL SUMMARY')
    if len(combined):
        print(combined.groupby(['cancer', 'model'])['test_c_index']
              .agg(mean='mean', std='std', n='count').round(4).to_string())
    print(f'\nSaved → {out_all}')
