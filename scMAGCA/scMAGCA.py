import numpy as np
import scipy.sparse as sp
import os

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import adjusted_mutual_info_score,normalized_mutual_info_score,adjusted_rand_score,silhouette_score,davies_bouldin_score,calinski_harabasz_score

import torch
import torch.nn as nn
from torch.nn import Parameter
import torch.nn.functional as F
import torch.optim as optim
from torch_sparse import SparseTensor
from torch_geometric.utils import add_self_loops,remove_self_loops,to_undirected
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

from scMAGCA.loss import *
from scMAGCA.layers import *
from scMAGCA.utils import *
from scMAGCA.Dataset import MyDataset

def dopca(X, dim=50):
    pcaten = PCA(n_components=dim,random_state=42)
    X_pca = pcaten.fit_transform(X)
    return X_pca

def get_adj(count, k=30, pca=100):
    if pca:
        countp = dopca(count, dim=pca)
    else:
        countp = count
    A = kneighbors_graph(countp, k, metric="cosine", include_self=False).toarray()
    adj = sp.coo_matrix(A)
    return adj  

class scMultiCluster(nn.Module):
    def __init__(self, input_dim1, input_dim2, alpha, beta, gama, device="cuda"):
        super(scMultiCluster, self).__init__()
        self.input_dim1 = input_dim1
        self.input_dim2 = input_dim2
        self.alpha = alpha
        self.beta = beta
        self.gama = gama
        self.device = device
        self.z_dim = 32
        # Graph Encoder
        self.encoder = Encoder(layer_config=[input_dim1+input_dim2,1024,256,64,self.z_dim])
        # Linear Decoder
        self.decoder = nn.Sequential(
            nn.Linear(in_features=self.z_dim, out_features=512),nn.BatchNorm1d(512),nn.PReLU(),
            nn.Linear(in_features=512, out_features=1024),nn.BatchNorm1d(1024),nn.PReLU(),
            nn.Linear(in_features=1024, out_features=input_dim1+input_dim2),
        ).to(self.device)
        self.dec_mean = nn.Sequential(nn.Linear(self.z_dim, 256),nn.Linear(256, 512),nn.Linear(512,input_dim1+input_dim2), MeanAct())
        self.dec_disp = nn.Sequential(nn.Linear(self.z_dim, 256),nn.Linear(256, 512),nn.Linear(512,input_dim1+input_dim2), DispAct())
        self.dec_pi = nn.Sequential(nn.Linear(self.z_dim, 256),nn.Linear(256, 512),nn.Linear(512,input_dim1+input_dim2), nn.Sigmoid())
        self.zinb_loss = ZINBLoss()

    def pretrain_autoencoder(self, X1, X2, X1_raw, X2_raw, epochs=400, file=None, ad_out=16):
        print("Pretraining stage")
        self.num = X1.shape[0]
        X1, X2 = torch.tensor(X1), torch.tensor(X2)
        self.X = torch.cat([X1, X2], dim=-1)
        X1_raw, X2_raw = torch.tensor(X1_raw), torch.tensor(X2_raw)
        self.X_raw1, self.X_raw2 = X1_raw.to(self.device), X2_raw.to(self.device)
        self.X_raw = torch.cat([X1_raw, X2_raw,], dim=-1).to(self.device)

        # 图邻接矩阵A(细胞连接图)
        adj = get_adj(self.X)
        edge_index = torch.tensor(np.array([adj.row, adj.col]), dtype=torch.long)
        edge_index, _ = remove_self_loops(edge_index)
        edge_index = to_undirected(edge_index, self.num)
        x = torch.tensor(self.X, dtype=torch.float32)
        y = None
        data = Data(x=x, edge_index=edge_index, y=y)
        data.edge_attr = torch.ones(data.edge_index.shape[1])
        nodes = torch.tensor(np.arange(data.num_nodes), dtype=torch.long)
        edge_index, edge_attr = add_self_loops(data.edge_index, data.edge_attr)
        data = Data(nodes=nodes, edge_index=data.edge_index, edge_attr=data.edge_attr, x=data.x, y=data.y,
                    num_nodes=data.num_nodes, neighbor_index=edge_index, neighbor_attr=edge_attr)

        data.to(self.device)
        if not os.path.exists('../datasets/'+file+'/raw'):
            os.makedirs('../datasets/'+file+'/raw')
        if not os.path.exists('../datasets/'+file+'/processed'):
            os.makedirs('../datasets/'+file+'/processed')       
        torch.save(data, '../datasets/'+file+'/raw/raw.pt')
        
        dataset = MyDataset(file=file)
        dataloader = DataLoader(dataset, batch_size=128, shuffle=True)
        optimizer = optim.Adam(filter(lambda p: p.requires_grad, self.parameters()), lr=1e-3, amsgrad=True)

        self.discriminator = Discriminator(self.z_dim,ad_out).to(self.device) #CITE-seq:16; SMAGE-seq:32
        self.real_distribution = torch.randn(self.num, self.z_dim).to(self.device)

        for epoch in range(epochs):
            for batch_data in dataloader:
                z = self.encoder(x=batch_data.x, 
                                edge_index=SparseTensor(row=batch_data.edge_index[0], col=batch_data.edge_index[1]), 
                                edge_weight=batch_data.edge_attr).to(self.device)

                # reconstruct_loss
                recon_x = self.decoder(z)
                recon_loss = F.mse_loss(recon_x, batch_data.x)

                # ZINB loss
                mean,disp,pi = self.dec_mean(z),self.dec_disp(z),self.dec_pi(z)
                zinb_loss = self.zinb_loss(x=self.X_raw, mean=mean, disp=disp, pi=pi)

                # adversarial loss
                self.dis_real_logit = self.discriminator(self.real_distribution)
                self.dis_fake_logit = self.discriminator(z)
                discriminator_loss = adversarial_loss(self.dis_real_logit,self.dis_fake_logit)

                loss = self.alpha*recon_loss + self.beta*zinb_loss + self.gama*discriminator_loss

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            print('Pretrain epoch {}, recon_loss:{:.6f}, zinb_loss:{:.6f}, adversial_loss:{:.6f}'.format(epoch+1, recon_loss.item(), zinb_loss.item(), discriminator_loss.item()))

        data.to(self.device)
        return z

    def fit(self, y=None, n_clusters=-1, num_epochs=400, file=None, tol=1e-3, pretrain_latent=None, resolution=.08):
        print("Clustering stage")
        dataset = MyDataset(file=file)
        dataloader = DataLoader(dataset, batch_size=128, shuffle=True)
        optimizer = optim.Adadelta(filter(lambda p: p.requires_grad, self.parameters()), lr=1, rho=.95)

        if y is not None:
            print("Initializing cluster centers with kmeans.")
            kmeans = KMeans(n_clusters=n_clusters)
            self.mu = Parameter(torch.Tensor(n_clusters, self.z_dim), requires_grad=True).to(self.device)
            for _, batch_data in enumerate(dataloader):
                z = self.encoder(x=batch_data.x, edge_index=SparseTensor(row=batch_data.edge_index[0], col=batch_data.edge_index[1]), edge_weight=batch_data.edge_attr)
            self.y_pred = kmeans.fit_predict(z.data.cpu().numpy())
            self.y_pred_last = self.y_pred
            self.mu.data.copy_(torch.Tensor(kmeans.cluster_centers_))
            
            ami = np.round(adjusted_mutual_info_score(y, self.y_pred), 5)
            nmi = np.round(normalized_mutual_info_score(y, self.y_pred), 5)
            ari = np.round(adjusted_rand_score(y, self.y_pred), 5)
            acc = np.round(cluster_acc(y, self.y_pred), 5)

            print('Initializing k-means: AMI= %.4f, NMI= %.4f, ARI= %.4f, ACC= %.4f' % (ami,nmi,ari,acc))
            self.train()
        
            final_nmi, final_ami, final_ari, final_acc= 0, 0, 0, 0
            for epoch in range(num_epochs):
                for _, batch_data in enumerate(dataloader):
                    z = self.encoder(x=batch_data.x, 
                                    edge_index=SparseTensor(row=batch_data.edge_index[0], col=batch_data.edge_index[1]), 
                                    edge_weight=batch_data.edge_attr)
                  
                    # reconstruct loss
                    recon_x = self.decoder(z)
                    recon_loss = F.mse_loss(recon_x, batch_data.x)

                    # zinb loss
                    mean,disp,pi = self.dec_mean(z),self.dec_disp(z),self.dec_pi(z)
                    zinb_loss = self.zinb_loss(x=self.X_raw, mean=mean, disp=disp, pi=pi)

                    # cluster loss
                    q = soft_assign(z, self.mu)
                    p = target_distribution(q).data
                    self.y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
                    cluste_loss = cluster_loss(p,q)

                    # overall loss
                    loss = self.alpha*recon_loss + self.beta*zinb_loss + cluste_loss

                    optimizer.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(self.mu, 1)
                    optimizer.step()

                self.y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
                final_ami = ami = np.round(adjusted_mutual_info_score(y, self.y_pred), 5)
                final_nmi = nmi = np.round(normalized_mutual_info_score(y, self.y_pred), 5)
                final_ari = ari = np.round(adjusted_rand_score(y, self.y_pred), 5)
                final_acc = acc = np.round(cluster_acc(y, self.y_pred), 5)

                print('Training epoch {}, recon_loss:{:.6f}, zinb_loss:{:.6f}, cluster_loss:{:.6f}'.format(epoch+1,recon_loss.item(),zinb_loss.item(),cluste_loss.item()))
                print('Clustering   %d: AMI= %.4f, NMI= %.4f, ARI= %.4f, ACC= %.4f' % (epoch+1,ami,nmi,ari,acc))
                     
                delta_label = np.sum(self.y_pred != self.y_pred_last).astype(np.float32) / self.num
                print(delta_label)
                self.y_pred_last, self.embedding = self.y_pred, z
                if epoch>0 and delta_label < tol:
                    print('delta_label ', delta_label, '< tol ', tol)
                    print("Reach tolerance threshold. Stopping training.")
                    break

            self.y_pred_last, self.embedding = self.y_pred, z
            print('Final Result : AMI= %.4f, NMI= %.4f, ARI= %.4f, ACC= %.4f'% (final_ami,final_nmi,final_ari,final_acc))         
                
            return self.y_pred_last, self.embedding

        else:
            n_clusters = GetCluster(pretrain_latent.cpu().detach().numpy(),n=20,res=resolution)
            print("Initializing cluster centers with kmeans.")
            kmeans = KMeans(n_clusters=n_clusters)
            self.mu = Parameter(torch.Tensor(n_clusters, self.z_dim), requires_grad=True).to(self.device)

            for _, batch_data in enumerate(dataloader):
                z = self.encoder(x=batch_data.x, edge_index=SparseTensor(row=batch_data.edge_index[0], col=batch_data.edge_index[1]), edge_weight=batch_data.edge_attr)
            self.y_pred = kmeans.fit_predict(z.data.cpu().numpy())
            self.y_pred_last = self.y_pred
            self.mu.data.copy_(torch.Tensor(kmeans.cluster_centers_))
            asw = np.round(silhouette_score(z.data.cpu().numpy(), self.y_pred), 5)
            db = np.round(davies_bouldin_score(z.data.cpu().numpy(), self.y_pred), 5)
            ch = np.round(calinski_harabasz_score(z.data.cpu().numpy(), self.y_pred), 5)

            print('Initializing k-means: ASW= %.4f, DB= %.4f, CH= %.4f' % (asw,db,ch))
            self.train()
        
            final_asw, final_db, final_ch = 0, 0, 0
            for epoch in range(num_epochs):
                for _, batch_data in enumerate(dataloader):
                    z = self.encoder(x=batch_data.x, 
                                    edge_index=SparseTensor(row=batch_data.edge_index[0], col=batch_data.edge_index[1]), 
                                    edge_weight=batch_data.edge_attr)
                  
                    # reconstruct loss
                    recon_x = self.decoder(z)
                    recon_loss = F.mse_loss(recon_x, batch_data.x)

                    # zinb loss
                    mean,disp,pi = self.dec_mean(z),self.dec_disp(z),self.dec_pi(z)
                    zinb_loss = self.zinb_loss(x=self.X_raw, mean=mean, disp=disp, pi=pi)
                    
                    # cluster loss
                    q = soft_assign(z, self.mu)
                    p = target_distribution(q).data
                    self.y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
                    cluste_loss = cluster_loss(p,q)

                    # overall loss
                    loss = self.alpha*recon_loss + self.beta*zinb_loss + cluste_loss

                    optimizer.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(self.mu, 1)
                    optimizer.step()
            
                self.y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
                final_asw = asw = np.round(silhouette_score(z.data.cpu().numpy(), self.y_pred), 5)
                final_db = db = np.round(davies_bouldin_score(z.data.cpu().numpy(), self.y_pred), 5)
                final_ch = ch = np.round(calinski_harabasz_score(z.data.cpu().numpy(), self.y_pred), 5)

                print('Training epoch {}, recon_loss:{:.6f}, zinb_loss:{:.6f}, cluster_loss:{:.6f}'.format(epoch+1,recon_loss.item(),zinb_loss.item(),cluste_loss.item()))
                print('Clustering   %d: ASW= %.4f, DB= %.4f, CH= %.4f' % (epoch+1,asw,db,ch))

            self.y_pred_last, self.embedding = self.y_pred, z
            print('Final Result : ASW= %.4f, DB= %.4f, CH= %.4f'% (final_asw,final_db,final_ch))  


            return self.y_pred_last, self.embedding

