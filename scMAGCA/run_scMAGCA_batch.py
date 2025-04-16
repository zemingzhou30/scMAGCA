import numpy as np
import pandas as pd
from time import time
import torch
import scanpy as sc
import random

from sklearn.preprocessing import OneHotEncoder
from scMAGCA.preprocess import read_dataset, preprocess_dataset
from scMAGCA.utils import *
from scMAGCA.scMAGCA_batch import scMultiCluster

if __name__ == "__main__":
    # setting the hyper parameters
    import argparse
    parser = argparse.ArgumentParser(description='train',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_clusters', type=int)
    parser.add_argument('--n_batches', type=int)
    parser.add_argument('--data_file', required=True)
    parser.add_argument('--pretrain_epochs', default=400, type=int)
    parser.add_argument('--maxiter', default=2000, type=int)
    parser.add_argument('--save_dir', default='../results/')
    parser.add_argument('--embedding_file', action='store_true', default=True)
    parser.add_argument('--prediction_file', action='store_true', default=True)
    parser.add_argument('--ad_out', type=int, default=16, help='Output dim of discriminator')
    parser.add_argument('--tol', default=0.001, type=float, help='Criterion to stop the model')
    parser.add_argument('--resolution', default=0.08, type=float, help='Resolution parameter to estimate K')
    parser.add_argument('--filter1', action='store_true', default=False, help='Do ADT/ATAC selection')
    parser.add_argument('--filter2', action='store_true', default=False, help='Do mRNA selection')
    parser.add_argument('--f1', default=2000, type=int, help='Number of ADT/ATAC after feature selection')
    parser.add_argument('--f2', default=2000, type=int, help='Number of mRNA after feature selection')
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

    ### read multiBatch dataset
    x1 = np.array(sc.read_h5ad('../datasets/GSE164378/GSE164378_adt.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/GSE164378/GSE164378_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/GSE164378/GSE164378_label.csv')['Cluster']).astype('float32')
    b = np.array(pd.read_csv('../datasets/GSE164378/GSE164378_batch.csv')['Batch']).astype('float32')
    enc = OneHotEncoder()
    enc.fit(b.reshape(-1, 1))
    B = enc.transform(b.reshape(-1, 1)).toarray()

    # Gene filter
    if args.filter1:
        importantGenes = geneSelection(x1, n=args.f1)
        x1 = x1[:, importantGenes]
    if args.filter2:
        importantGenes = geneSelection(x2, n=args.f2)
        x2 = x2[:, importantGenes]
        
    # preprocessing scRNA-seq read counts matrix
    adata1 = sc.AnnData(x1)
    adata1 = read_dataset(adata1, copy=True)
    adata1 = preprocess_dataset(adata1, normalize_input=True, logtrans_input=True)
    print(adata1)

    adata2 = sc.AnnData(x2)
    adata2 = read_dataset(adata2, copy=True)
    adata2 = preprocess_dataset(adata2, normalize_input=True, logtrans_input=True)
    print(adata2)
    
    input_size1 = adata1.n_vars
    input_size2 = adata2.n_vars

    model = scMultiCluster(input_dim1=input_size1,input_dim2=input_size2,n_batch=args.n_batches, 
                           alpha=args.alpha,beta=args.beta,gama=args.gama,device=args.device).to(args.device)
    print(str(model))

    t0 = time()
    pretrain_latent = model.pretrain_autoencoder(
                        X1=adata1.X, X2=adata2.X, X1_raw=adata1.raw.X, X2_raw=adata2.raw.X, B=B,
                        epochs=args.pretrain_epochs, file=args.data_file)
    print('Pretraining time: %d seconds.' % int(time() - t0))
 
    y_pred, final_latent = model.fit(y=y, n_clusters=args.n_clusters, num_epochs=args.maxiter, file=args.data_file)
    if args.prediction_file:
        y_pred_ = best_map(y, y_pred)
        y_pred_ = y_pred_.astype(int)
        np.savetxt(args.save_dir + str(args.data_file) + "_pred.csv", y_pred, delimiter=",")
    if args.embedding_file:
       final_latent = final_latent.cpu().detach().numpy()
       np.savetxt(args.save_dir + str(args.data_file) + "_embedding.csv", final_latent, delimiter=",")

    print('Total time: %d seconds.' % int(time() - t0))