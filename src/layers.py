import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv

class MeanAct(nn.Module):
    def __init__(self):
        super(MeanAct, self).__init__()

    def forward(self, x):
        return torch.clamp(torch.exp(x), min=1e-5, max=1e6)

class DispAct(nn.Module):
    def __init__(self):
        super(DispAct, self).__init__()

    def forward(self, x):
        return torch.clamp(F.softplus(x), min=1e-4, max=1e4)
    
class Encoder(nn.Module):
    def __init__(self, layer_config):
        super().__init__()
        self.stacked_gnn = nn.ModuleList(
            [GCNConv(layer_config[i - 1], layer_config[i]) for i in range(1, len(layer_config))])
        self.stacked_bns = nn.ModuleList(
            [nn.BatchNorm1d(layer_config[i], momentum=0.01) for i in range(1, len(layer_config))])
        self.stacked_prelus = nn.ModuleList([nn.PReLU() for _ in range(1, len(layer_config))])

    def forward(self, x, edge_index, edge_weight=None):
        for i, gnn in enumerate(self.stacked_gnn):
            x = gnn(x, edge_index, edge_weight)
            x = self.stacked_bns[i](x)
            x = self.stacked_prelus[i](x)

        return x
    
def buildNetwork(layers, activation="elu"):
    net = nn.Sequential()
    for i in range(1, len(layers)):
        net.add_module('linear%d'%i, nn.Linear(layers[i-1], layers[i]))
        net.add_module('batchnorm%d'%i, nn.BatchNorm1d(layers[i], affine=True))
        if activation=="relu":
            net.add_module('relu%d'%i, nn.ReLU())
        elif activation=="selu":
            net.add_module('selu%d'%i, nn.SELU())
        elif activation=="sigmoid":
            net.add_module('sigmoid%d'%i, nn.Sigmoid())
        elif activation=="elu":
            net.add_module('elu%d'%i, nn.ELU())
    return net

class Discriminator(nn.Module):
    def __init__(self, latent_dim, hidden_dim):
        super(Discriminator, self).__init__()
        self.fc1 = nn.Linear(latent_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, 1)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = torch.sigmoid(self.fc2(x))
        return x