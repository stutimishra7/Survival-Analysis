"""
Training loop and repeated evaluation over independent random splits.
"""

import copy

import numpy as np
import torch
import torch.optim as optim

from rcoxnet.config import (
    N_REPEATS, TEST_FRAC, VAL_FRAC, EPOCHS, CHECKPOINT, GRAD_CLIP,
)
from rcoxnet.model import CoxPhRWRNet, cox_loss, c_index
from rcoxnet.preprocessing import (
    censoring_time_balanced_split,
    normalise_genomic,
    normalise_clinical,
    prep,
)


def train_model(device, dtype, model_kwargs,
                x_tr, cl_tr, yt_tr, ye_tr,
                x_va, cl_va, yt_va, ye_va,
                x_te, cl_te, yt_te, ye_te,
                seed, lr, l2,
                epochs=EPOCHS, checkpoint=CHECKPOINT):
    """
    Train one model run. Saves the checkpoint with the best validation C-index
    and returns test C-index and best validation C-index.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    net = CoxPhRWRNet(**model_kwargs).to(device)
    opt = optim.Adam(net.parameters(), lr=lr, weight_decay=l2)
    best_val, best_state = 0.0, copy.deepcopy(net.state_dict())

    for ep in range(epochs):
        net.train()
        opt.zero_grad()
        loss = cox_loss(net(x_tr, cl_tr), yt_tr, ye_tr)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(net.parameters(), GRAD_CLIP)
        opt.step()

        if ep % checkpoint == 0:
            net.eval()
            with torch.no_grad():
                vc = c_index(net(x_va, cl_va), yt_va, ye_va).item()
            if vc > best_val:
                best_val, best_state = vc, copy.deepcopy(net.state_dict())

    net.load_state_dict(best_state)
    net.eval()
    with torch.no_grad():
        tc = c_index(net(x_te, cl_te), yt_te, ye_te).item()
    return tc, best_val


def train_and_eval(device, dtype, X_data, clin_data, yt, ye,
                   model_kwargs, lr, l2, n_repeats=N_REPEATS):
    """
    Run training over n_repeats independent random splits and collect
    test and validation C-index for each repeat.
    """
    test_cis, val_cis = [], []

    for rep in range(n_repeats):
        rng = np.random.RandomState(rep * 7 + 13)
        tr, va, te = censoring_time_balanced_split(yt, ye, TEST_FRAC, VAL_FRAC, rng)

        X_trn, X_van, X_ten    = normalise_genomic(X_data[tr], X_data[va], X_data[te])
        cl_trn, cl_van, cl_ten = normalise_clinical(clin_data[tr], clin_data[va], clin_data[te])

        x_tr_t, cl_tr_t, yt_tr_t, ye_tr_t = prep(device, dtype, X_trn, cl_trn, yt[tr], ye[tr])
        x_va_t, cl_va_t, yt_va_t, ye_va_t = prep(device, dtype, X_van, cl_van, yt[va], ye[va])
        x_te_t, cl_te_t, yt_te_t, ye_te_t = prep(device, dtype, X_ten, cl_ten, yt[te], ye[te])

        tc, vc = train_model(
            device, dtype, model_kwargs,
            x_tr_t, cl_tr_t, yt_tr_t, ye_tr_t,
            x_va_t, cl_va_t, yt_va_t, ye_va_t,
            x_te_t, cl_te_t, yt_te_t, ye_te_t,
            seed=rep, lr=lr, l2=l2,
        )
        test_cis.append(tc)
        val_cis.append(vc)

    return test_cis, val_cis
