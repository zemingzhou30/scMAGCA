import numpy as np
import pandas as pd
from time import time
import torch
import scanpy as sc
import random
import os

from scMAGCA.preprocess import read_dataset, preprocess_dataset
from scMAGCA.utils import *
from scMAGCA.scMAGCA_3omics import scMultiCluster

if __name__ == "__main__":

    # setting the hyper parameters
    import argparse
    parser = argparse.ArgumentParser(description='train',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_clusters', default=-1, type=int)
    parser.add_argument('--data_file', required=True)
    parser.add_argument('--pretrain_epochs', default=400, type=int)
    parser.add_argument('--maxiter', default=2000, type=int)
    parser.add_argument('--save_dir', default='../results/')
    parser.add_argument('--embedding_file', action='store_true', default=True)
    parser.add_argument('--prediction_file', action='store_true', default=True)
    parser.add_argument('--ad_out', type=int, default=16, help='Output dim of discriminator')
    parser.add_argument('--tol', default=0.001, type=float, help='Criterion to stop the model')
    parser.add_argument('--resolution', default=0.1, type=float, help='Resolution parameter to estimate K')
    parser.add_argument('--filter1', action='store_true', default=False, help='Do ADT selection')
    parser.add_argument('--filter2', action='store_true', default=False, help='Do ATAC selection')
    parser.add_argument('--filter3', action='store_true', default=False, help='Do mRNA selection')
    parser.add_argument('--f1', default=2000, type=int, help='Number of ADT after feature selection')
    parser.add_argument('--f2', default=2000, type=int, help='Number of ATAC after feature selection')
    parser.add_argument('--f3', default=2000, type=int, help='Number of mRNA after feature selection')
    parser.add_argument('--alpha', default=0.2, type=float, help='Weight of reconstruction loss')
    parser.add_argument('--beta', default=0.8, type=float, help='Weight of ZINB loss')
    parser.add_argument('--gama', default=0.01, type=float, help='Weight of adversarial loss')
    parser.add_argument('--device', default='cuda')
    args = parser.parse_args()
    print(args)

    # set seed
    random.seed(3407)
    np.random.seed(3407)
    torch.manual_seed(3407)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False

    ### Read dataset
    x1 = sc.read_h5ad('../datasets/GSE158013/GSE158013_adt.h5ad').layers['counts']
    feature1 = sc.read_h5ad('../datasets/GSE158013/GSE158013_adt.h5ad').var.index
    x2 = sc.read_h5ad('../datasets/GSE158013/GSE158013_atac.h5ad').layers['counts'].A
    feature2 = sc.read_h5ad('../datasets/GSE158013/GSE158013_atac.h5ad').var.index
    x3 = sc.read_h5ad('../datasets/GSE158013/GSE158013_rna.h5ad').layers['counts'].A
    feature3 = sc.read_h5ad('../datasets/GSE158013/GSE158013_rna.h5ad').var.index
    y = None

    # Gene filter
    if args.filter1:
        importantGenes = geneSelection(x1, n=args.f1)
        x1 = x1[:, importantGenes]
        feature1 = feature1[importantGenes]
    if args.filter2:
        importantGenes = geneSelection(x2, n=args.f2)
        x2 = x2[:, importantGenes]
        feature2 = feature2[importantGenes]
    if args.filter3:
        importantGenes = geneSelection(x3, n=args.f3)
        x3 = x3[:, importantGenes]
        feature3 = feature3[importantGenes]
        
    # preprocessing scRNA-seq read counts matrix
    adata1 = sc.AnnData(x1)
    adata1 = read_dataset(adata1, copy=True)
    adata1 = preprocess_dataset(adata1, normalize_input=True, logtrans_input=True)
    adata1.var['importantGenes'] = feature1
    print(adata1)

    adata2 = sc.AnnData(x2)
    adata2 = read_dataset(adata2, copy=True)
    adata2 = preprocess_dataset(adata2, normalize_input=True, logtrans_input=True)
    adata2.var['importantGenes'] = feature2
    print(adata2)

    adata3 = sc.AnnData(x3)
    adata3 = read_dataset(adata3, copy=True)
    adata3 = preprocess_dataset(adata3, normalize_input=True, logtrans_input=True)
    adata3.var['importantGenes'] = feature3
    print(adata3)
    
    input_size1 = adata1.n_vars
    input_size2 = adata2.n_vars
    input_size3 = adata3.n_vars

    model = scMultiCluster(input_dim1=input_size1,input_dim2=input_size2,input_dim3=input_size3,
                           alpha=args.alpha,beta=args.beta,gama=args.gama,device=args.device).to(args.device)
    print(str(model))

    if not os.path.exists(args.save_dir):
        print("Create save_dir")
        os.makedirs(args.save_dir)

    t0 = time()
    pretrain_latent = model.pretrain_autoencoder(
                        X1=adata1.X, X2=adata2.X, X3=adata3.X, X1_raw=adata1.raw.X, X2_raw=adata2.raw.X, X3_raw=adata3.raw.X, 
                        epochs=args.pretrain_epochs, file=args.data_file, ad_out=args.ad_out)
    print('Pretraining time: %d seconds.' % int(time() - t0))

    if y is not None:   
        y_pred, final_latent = model.fit(y=y, n_clusters=args.n_clusters, num_epochs=args.maxiter, file=args.data_file, tol=args.tol)
        if args.prediction_file:
            y_pred_ = best_map(y, y_pred)
            y_pred_ = y_pred_.astype(int)
            np.savetxt(args.save_dir + str(args.data_file) + "_pred.csv", y_pred, delimiter=",")
    else:
        y_pred, final_latent = model.fit(y=y, n_clusters=args.n_clusters, num_epochs=args.maxiter, file=args.data_file, tol=args.tol, pretrain_latent=pretrain_latent, resolution=args.resolution)
        if args.prediction_file:
            np.savetxt(args.save_dir + str(args.data_file) +  "_pred.csv", y_pred, delimiter=",")

    if args.embedding_file:
       final_latent = final_latent.cpu().detach().numpy()
       np.savetxt(args.save_dir + str(args.data_file) + "_embedding.csv", final_latent, delimiter=",")

    print('Total time: %d seconds.' % int(time() - t0))