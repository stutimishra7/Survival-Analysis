
import torch
import torch.nn as nn

class CoxPhRWRNet(nn.Module):
    def __init__(self, input_nodes, hidden_nodes1, hidden_nodes2, output_nodes):
        super(CoxPhRWRNet, self).__init__()
        self.tanh = nn.Tanh()

        self.rwr_layer = nn.Linear(input_nodes, hidden_nodes1)
        self.hidden_layer1 = nn.Linear(hidden_nodes1, hidden_nodes2)
        self.hidden_layer2 = nn.Linear(hidden_nodes2, output_nodes, bias=False)
        self.cox_layer = nn.Linear(output_nodes + 3, 1, bias=False)  # Additional features (age, msi, tmb)
        self.cox_layer.weight.data.uniform_(-0.001, 0.001)

    def forward(self, x_genomic, x_msi, x_tmb, x_age):
        x_rwr = self.tanh(self.rwr_layer(x_genomic))
        x_hidden1 = self.tanh(self.hidden_layer1(x_rwr))
        x_hidden2 = self.tanh(self.hidden_layer2(x_hidden1))

        # Combine with additional features
        x_combined = torch.cat((x_hidden2, x_msi, x_tmb, x_age), 1)
        cox_output = self.cox_layer(x_combined)

        return cox_output
