import torch
import torch.nn as nn
import torch.nn.functional as F

tau = 1
alpha = 1
EPS = 1e-15

# NB Loss
class NBLoss(nn.Module):
    def __init__(self):
        super(NBLoss, self).__init__()

    def forward(self, x, mean, disp, scale_factor=1.0):
        eps = 1e-10
        #scale_factor = scale_factor[:, None]
        mean = mean * scale_factor
        
        t1 = torch.lgamma(disp+eps) + torch.lgamma(x+1.0) - torch.lgamma(x+disp+eps)
        t2 = (disp+x) * torch.log(1.0 + (mean/(disp+eps))) + (x * (torch.log(disp+eps) - torch.log(mean+eps)))
        result = t1 + t2

        result = torch.mean(result)
        return result

# ZINB Loss
class ZINBLoss(nn.Module):
    def __init__(self):
        super(ZINBLoss, self).__init__()

    def forward(self, x, mean, disp, pi, scale_factor=1.0, ridge_lambda=0.0):
        eps = 1e-10
        #scale_factor = scale_factor[:, None]
        mean = mean * scale_factor

        t1 = torch.lgamma(disp + eps) + torch.lgamma(x + 1.0) - torch.lgamma(x + disp + eps)
        t2 = (disp + x) * torch.log(1.0 + (mean / (disp + eps))) + (x * (torch.log(disp + eps) - torch.log(mean + eps)))
        nb_final = t1 + t2

        nb_case = nb_final - torch.log(1.0 - pi + eps)
        zero_nb = torch.pow(disp / (disp + mean + eps), disp)
        zero_case = -torch.log(pi + ((1.0 - pi) * zero_nb) + eps)
        result = torch.where(torch.le(x, 1e-8), zero_case, nb_case)

        if ridge_lambda > 0:
            ridge = ridge_lambda * torch.square(pi)
            result += ridge

        result = torch.mean(result)
        return result   

def soft_assign(z, mu):
    q = 1.0 / (1.0 + torch.sum((z.unsqueeze(1) - mu) ** 2, dim=2) / alpha)
    q = q ** ((alpha + 1.0) / 2.0)
    q = (q.t() / torch.sum(q, dim=1)).t()
    return q

def target_distribution(q):
    p = q ** 2 / q.sum(0)
    return (p.t() / p.sum(1)).t() 

# Cluster Loss
def cluster_loss(p, q):
    def kld(target, pred):
        return torch.mean(torch.sum(target * torch.log(target / (pred + 1e-6)), dim=-1))
    kldloss = kld(p, q)
    return kldloss
    
# Adversarial Loss
def adversarial_loss(real_output, fake_output):
    real_loss = F.binary_cross_entropy(real_output, torch.ones_like(real_output))
    fake_loss = F.binary_cross_entropy(fake_output, torch.zeros_like(fake_output))
    total_loss = real_loss + fake_loss
    return total_loss