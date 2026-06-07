"""
RCoxNet — Random-Walk-with-Restart Cox Proportional-Hazards Network
====================================================================
Package layout
--------------
rcoxnet/
    config.py        — cancer configs, shared constants
    data_loader.py   — raw file loading
    mutations.py     — functional-variant parsing
    clinical.py      — clinical feature extraction
    network.py       — PPI transition-matrix construction
    rwr.py           — Random Walk with Restart
    features.py      — feature-matrix assembly, log-rank selection
    preprocessing.py — normalisation, train/val/test splitting
    model.py         — CoxPhRWRNet, cox_loss, c_index
    train.py         — single-run and repeated training
    grid_search.py   — hyperparameter and (alpha, K) grid search
    pipeline.py      — end-to-end orchestrator
"""

try:
    from rcoxnet.pipeline import run_pipeline
    __all__ = ["run_pipeline"]
except ModuleNotFoundError:
    # torch not available in this environment (e.g. coxnnet, survnet)
    pass
