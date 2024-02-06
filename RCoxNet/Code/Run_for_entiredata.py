
from Run_train import *

from Train_for_entiredata import InterpretCoxPhRWRNet
import pandas as pd 
from model import CoxPhRWRNet
import torch
import numpy as np
# Optimal L2:  0.005 Optimal LR:  0.01 GBM
# C-index in Test:  tensor(0.6500, device='cuda:0')
dtype = torch.FloatTensor
''' Net Settings'''
In_Nodes = 549 ###number of genes
hidden_nodes1 = 860 ###number of pathways
hidden_nodes2 = 150 ###number of hidden nodes
Out_Nodes = 30 ###number of hidden nodes in the last hidden layer
''' Initialize '''
Initial_Learning_Rate = opt_lr_loss
L2_Lambda = opt_l2_loss
Num_EPOCHS = 3

###sub-network setup
Dropout_Rate = [0.7, 0.5]

''' load data and pathway '''
# pathway_mask = load_pathway("pathway_mask_tcga_ov.csv", dtype)
x, ytime, yevent, age,msi,tmb = load_sorted_data("entire_data_tcga_Brca_all_cosmic_tmb.csv", dtype)

outpath = "InterpretCoxPASNet_lin_pred_brca_nopath.pt"
'''train Cox-PASNet for model interpretation'''
InterpretCoxPhRWRNet(x, age,msi,tmb, ytime, yevent, \
                   In_Nodes, hidden_nodes1 , hidden_nodes2, Out_Nodes, \
					Initial_Learning_Rate, L2_Lambda, Num_EPOCHS, outpath)

'''load trained Cox-PASNet'''
net = CoxPhRWRNet(In_Nodes,hidden_nodes1, hidden_nodes2, Out_Nodes)
net.load_state_dict(torch.load(outpath))
###if gpu is being used
if torch.cuda.is_available():
	net.cuda()
###
'''save weights and node values into files individually'''
w_sc1 = net.rwr_layer.weight.data.cpu().detach().numpy()
w_sc2 = net.hidden_layer1.weight.data.cpu().detach().numpy()
w_sc3 = net.hidden_layer2 .weight.data.cpu().detach().numpy()
w_sc4 = net.cox_layer.weight.data.cpu().detach().numpy()
np.savetxt("w_sc1_lin_pred_brca_nopath_150.csv", w_sc1, delimiter = ",")
np.savetxt("w_sc2_lin_pred_brca_nopath_150.csv", w_sc2, delimiter = ",")
np.savetxt("w_sc3_lin_pred_brca_nopath_150.csv", w_sc3, delimiter = ",")
np.savetxt("w_sc4_lin_pred_brca_nopath_150.csv", w_sc4, delimiter = ",")

#  self.rwr_layer = nn.Linear(input_nodes, hidden_nodes1)
#         self.hidden_layer1 = nn.Linear(hidden_nodes1, hidden_nodes2)
#         self.hidden_layer2 = nn.Linear(hidden_nodes2, output_nodes, bias=False)
#         self.cox_layer = nn.Linear(output_nodes + 3, 1, bias=False)  # Additional features (age, msi, tmb)
#         self.cox_layer.weight.data.uniform_(-0.001, 0.001)
pathway_node = net.tanh(net.rwr_layer(x))
hidden_node = net.tanh(net.hidden_layer1(pathway_node))
hidden_2_node = net.tanh(net.hidden_layer2(hidden_node))
x_cat = torch.cat((hidden_2_node,age,msi,tmb), 1)
lin_pred = net.cox_layer(x_cat)
output_data = np.column_stack((lin_pred.cpu().detach().numpy(), ytime.cpu().detach().numpy(), yevent.cpu().detach().numpy(), age.cpu().detach().numpy(), msi.cpu().detach().numpy(), tmb.cpu().detach().numpy()))

# Save the concatenated data to lin_pred_OV_AN.csv
np.savetxt("lin_pred_BRCA_TEST.csv", output_data, delimiter=",")
lin_pred = net.cox_layer(x_cat)
np.savetxt("pathway_lin_pred_brca_nopath_150.csv", pathway_node.cpu().detach().numpy(), delimiter = ",")
np.savetxt("hidden_node_lin_pred_brca_nopath_150.csv", hidden_node.cpu().detach().numpy(), delimiter = ",")
np.savetxt("hidden_2_lin_pred_brca_nopath_150.csv", x_cat.cpu().detach().numpy(), delimiter = ",")
# np.savetxt("lin_pred_OV_AN.csv", lin_pred.cpu().detach().numpy(), delimiter = ",")
columns = ['PI','OS_MONTHS','OS_EVENT','AGE','tmb_score',	'msi_score']
surv_data = pd.DataFrame(output_data , columns=columns)
