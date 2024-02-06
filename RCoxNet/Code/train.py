
# %%writefile CostFunc_CIndex.pyimport torch
from CostFunc_CIndex import create_indicator_matrix, concordance_index,  negative_partial_log_likelihood
from model import  CoxPhRWRNet
import torch
import torch.optim as optim
import copy
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd

dtype = torch.FloatTensor


import torch
import torch.optim as optim
import copy
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd

dtype = torch.FloatTensor

def train_cox_rwrnnet(train_x, train_age, train_ytime, train_yevent, train_msi, train_tmb, eval_x, eval_age, eval_ytime, eval_yevent, eval_msi, eval_tmb, In_Nodes, hidden_nodes1, hidden_nodes2, Out_Nodes, Learning_Rate, L2, Num_Epochs):

    net = CoxPhRWRNet(In_Nodes,hidden_nodes1, hidden_nodes2, Out_Nodes)
    ###if gpu is being used
    if torch.cuda.is_available():
        net.cuda()
    ###
    ###optimizer
    opt = optim.Adam(net.parameters(), lr=Learning_Rate, weight_decay=L2)

    for epoch in range(Num_Epochs + 1):
        net.train()
        opt.zero_grad()  ###reset gradients to zeros

        pred = net(train_x, train_age, train_msi, train_tmb)  ###Forward
        loss = negative_partial_log_likelihood(pred, train_ytime, train_yevent)  ###calculate loss
        loss.backward()  ###calculate gradients
        opt.step()  ###update weights and biases

        if epoch % 200 == 0:
            net.train()
            train_pred = net(train_x, train_age, train_msi, train_tmb)
            train_loss =negative_partial_log_likelihood(train_pred, train_ytime, train_yevent).view(1,)

            net.eval()
            eval_pred = net(eval_x, eval_age, eval_msi, eval_tmb)
            eval_loss = negative_partial_log_likelihood(eval_pred, eval_ytime, eval_yevent).view(1,)

            train_cindex = concordance_index(train_pred, train_ytime, train_yevent)
            eval_cindex = concordance_index(eval_pred, eval_ytime, eval_yevent)
            print("Loss in Train: ", train_loss)

    return (train_loss, eval_loss, train_cindex, eval_cindex)
