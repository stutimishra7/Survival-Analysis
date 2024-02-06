import torch
from Run_train import *
from model import CoxPhRWRNet
import torch.optim as optim
import copy
from scipy.interpolate import interp1d
import numpy as np

dtype = torch.FloatTensor

def InterpretCoxPhRWRNet(x, age, msi, tmb, ytime, yevent,
                        In_Nodes, hidden_nodes1 , hidden_nodes2, Out_Nodes,
                        Learning_Rate, L2, Num_Epochs, outpath):

    net = CoxPhRWRNet(In_Nodes,hidden_nodes1 , hidden_nodes2, Out_Nodes)
    ###if gpu is being used
    if torch.cuda.is_available():
        net.cuda()
    ###
    ###optimizer
    opt = optim.Adam(net.parameters(), lr=Learning_Rate, weight_decay=L2)

    for epoch in range(Num_Epochs + 1):
        net.train()
        opt.zero_grad()  ###reset gradients to zeros
        ###Randomize dropout masks
        # net.do_m1 = dropout_mask(860, Dropout_Rate[0])
        # net.do_m2 = dropout_mask(Hidden_Nodes, Dropout_Rate[1])

        pred = net(x, age, msi, tmb)  ###Forward
        loss = negative_partial_log_likelihood(pred, ytime, yevent)  ###calculate loss
        loss.backward()  ###calculate gradients
        opt.step()  ###update weights and biases
        print( loss)


        net_state_dict = net.state_dict()


    torch.save(net.state_dict(), outpath)

    return
