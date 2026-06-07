"""
CoxPhRWRNet: Cox PH network that takes RWR genomic features and clinical
covariates (AGE, MSI, TMB) as input and outputs a scalar risk score.
"""

import torch
import torch.nn as nn


def _tril_indicator(x):
    n = x.size(0)
    return torch.tril(torch.ones(n, n, device=x.device))


def cox_loss(pred, ytime, yevent):
    # normalised Cox partial log-likelihood; inputs sorted by descending time
    n_obs    = yevent.sum()
    risk_set = _tril_indicator(ytime)
    shift    = pred.detach().max()
    log_risk = torch.log(risk_set.mm(torch.exp(pred - shift))) + shift
    diff     = pred - log_risk
    sum_obs  = diff.t().mm(yevent)
    raw_loss = -(sum_obs / n_obs)
    baseline = (torch.log(risk_set.sum(dim=1, keepdim=True)) * yevent).sum() / n_obs
    return raw_loss / baseline


def c_index(pred, ytime, yevent):
    # Harrell's C-index; censored patients excluded from numerator
    yt_mat = _tril_indicator(ytime)
    yt_mat = yt_mat - torch.diag(torch.diag(yt_mat))
    censor = (yevent.squeeze() == 0).nonzero(as_tuple=False).squeeze()
    yt_mat[censor, :] = 0

    pred_col = pred.view(-1, 1)
    pred_row = pred.view(1, -1)
    pred_mat = ((pred_col > pred_row).float() +
                0.5 * (pred_col == pred_row).float())
    return torch.div(torch.sum(pred_mat.mul(yt_mat)), torch.sum(yt_mat))


class CoxPhRWRNet(nn.Module):
    """
    Three-layer network for Cox PH survival prediction.

    Genomic branch: input -> hidden1 -> hidden2 -> latent (all Tanh).
    Final layer: concat(latent, clinical) -> risk score (no bias, no activation).
    """

    def __init__(self, input_nodes, hidden_nodes1, hidden_nodes2,
                 output_nodes, n_clin=3):
        super().__init__()
        self.tanh          = nn.Tanh()
        self.rwr_layer     = nn.Linear(input_nodes,   hidden_nodes1)
        self.hidden_layer1 = nn.Linear(hidden_nodes1, hidden_nodes2)
        self.hidden_layer2 = nn.Linear(hidden_nodes2, output_nodes, bias=False)
        self.cox_layer     = nn.Linear(output_nodes + n_clin, 1, bias=False)
        self.cox_layer.weight.data.uniform_(-0.001, 0.001)

    def forward(self, x_genomic, x_clin):
        x = self.tanh(self.rwr_layer(x_genomic))
        x = self.tanh(self.hidden_layer1(x))
        x = self.tanh(self.hidden_layer2(x))
        return self.cox_layer(torch.cat([x, x_clin], dim=1))
