"""
RCoxNet — main entry point.

Usage
-----
# Reproduce paper results (uses best_params.json, skips grid search):
    python run.py --reproduce --cancer BRCA
    python run.py --reproduce --cancer ALL

# Run full pipeline with joint grid search (α × K × arch × LR × L2):
    python run.py --cancer BRCA
    python run.py --cancer ALL
"""

import argparse
import json
import os

import numpy as np
import torch

CANCERS = ['BRCA', 'LUNG', 'GBM', 'OV']

def reproduce(cancer: str, params: dict, root: str):
    """Run 20-repeat evaluation using paper-reported hyperparameters."""
    from rcoxnet.pipeline import load_data, build_network, train_final

    print(f'\n{cancer}  (reproduce mode)')
    print(f'  alpha={params["alpha"]}  K={params["K"]}  '
          f'arch=({params["hidden1"]},{params["hidden2"]},{params["output_nodes"]})  '
          f'lr={params["lr"]}  l2={params["l2"]}')

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dtype  = torch.FloatTensor

    data    = load_data(cancer, root)
    network = build_network(data['ppi_path'])

    best = dict(
        alpha       = params['alpha'],
        K           = params['K'],
        hidden1     = params['hidden1'],
        hidden2     = params['hidden2'],
        output_nodes= params['output_nodes'],
        lr          = params['lr'],
        l2          = params['l2'],
    )
    test_cis, val_cis = train_final(
        device, dtype, data['clinical'], data['seed_by_patient'], network, best)

    print(f'\n  {cancer}  C-index: {np.mean(test_cis):.3f} ± {np.std(test_cis):.3f}')
    return test_cis

def full_search(cancer: str, root: str):
    """Run full pipeline with joint hyperparameter grid search."""
    from rcoxnet.pipeline import run_pipeline
    res = run_pipeline(cancer, root=root)
    print(f'\n  {cancer}  C-index: {res["test_ci_mean"]:.3f} ± {res["test_ci_std"]:.3f}')
    return res

def main():
    parser = argparse.ArgumentParser(
        description='RCoxNet — cancer survival prediction')
    parser.add_argument('--cancer', default='ALL',
                        choices=CANCERS + ['ALL'],
                        help='cancer type to run  (default: ALL)')
    parser.add_argument('--reproduce', action='store_true',
                        help='skip grid search; use hyperparameters from best_params.json')
    parser.add_argument('--root', default=None,
                        help='repository root directory  (default: location of this file)')
    args = parser.parse_args()

    root    = args.root or os.path.dirname(os.path.abspath(__file__))
    targets = CANCERS if args.cancer == 'ALL' else [args.cancer]

    if args.reproduce:
        params_path = os.path.join(root, 'best_params.json')
        with open(params_path) as f:
            best_params = json.load(f)
        for cancer in targets:
            reproduce(cancer, best_params[cancer], root)
    else:
        for cancer in targets:
            full_search(cancer, root)

if __name__ == '__main__':
    main()
