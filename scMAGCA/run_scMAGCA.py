import numpy as np
import pandas as pd
from time import time
import torch
import scanpy as sc
import random
import os

from scMAGCA.preprocess import read_dataset, preprocess_dataset
from scMAGCA.utils import *
from scMAGCA.scMAGCA import scMultiCluster

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

    ### read CITE-Seq datasets
    x1 = np.array(sc.read_h5ad('../datasets/10X1kpbmc/10X1kpbmc_adt.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/10X1kpbmc/10X1kpbmc_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/10X1kpbmc/10X1kpbmc_label.csv')['Cluster']).astype('float32')
    '''x1 = np.array(sc.read_h5ad('../datasets/10X5kpbmc/10X5kpbmc_adt.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/10X5kpbmc/10X5kpbmc_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/10X5kpbmc/10X5kpbmc_label.csv')['Cluster']).astype('float32')'''
    '''x1 = np.array(sc.read_h5ad('../datasets/10X5kpbmcTotalSeq/10X5kpbmcTotalSeq_adt.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/10X5kpbmcTotalSeq/10X5kpbmcTotalSeq_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/10X5kpbmcTotalSeq/10X5kpbmcTotalSeq_label.csv')['Cluster']).astype('float32')'''
    '''x1 = np.array(sc.read_h5ad('../datasets/10Xmalt/10Xmalt_adt.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/10Xmalt/10Xmalt_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/10Xmalt/10Xmalt_label.csv')['Cluster']).astype('float32')'''
    '''data_mat = h5py.File('../datasets/GSE128639/GSE128639.h5')
    x1, x2, y = np.array(data_mat['X2']).astype('float32'), np.array(data_mat['X1']).astype('float32'), np.array(data_mat['Y']).astype('float32')
    data_mat.close()'''
    '''data_mat = h5py.File('../datasets/spector/spector.h5')
    x1, x2, y = np.array(data_mat['X2']).astype('float32'), np.array(data_mat['X1']).astype('float32'), np.array(data_mat['Y']).astype('float32')'''

    ### read SMAGE-Seq datasets
    '''x1 = np.array(sc.read_h5ad('../datasets/human_brain_3k/human_brain_3k_atac.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/human_brain_3k/human_brain_3k_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/human_brain_3k/human_brain_3k_label.csv')['Cluster']).astype('float32')'''
    '''x1 = np.array(sc.read_h5ad('../datasets/human_pbmc_3k/human_pbmc_3k_atac.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/human_pbmc_3k/human_pbmc_3k_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/human_pbmc_3k/human_pbmc_3k_label_a.csv')['Cluster']).astype('float32')'''
    '''x1 = np.array(sc.read_h5ad('../datasets/mouse_brain_5k/mouse_brain_5k_atac.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/mouse_brain_5k/mouse_brain_5k_rna.h5ad').to_df()).astype('float32')
    y = np.array(pd.read_csv('../datasets/mouse_brain_5k/mouse_brain_5k_label_a.csv')['Cluster']).astype('float32')'''
    ### read no-label SMAGE-Seq datasets
    '''x1 = np.array(sc.read_h5ad('../datasets/GSM4949911/GSM4949911_atac.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/GSM4949911/GSM4949911_rna.h5ad').to_df()).astype('float32')
    y = None'''
    '''x1 = np.array(sc.read_h5ad('../datasets/pbmc_10X_public/pbmc_10X_atac_public.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/pbmc_10X_public/pbmc_10X_rna_public.h5ad').to_df()).astype('float32')
    y = None'''
    '''x1 = np.array(sc.read_h5ad('../datasets/10X-Multiome-Pbmc10k/10X-Multiome-Pbmc10k-ATAC.h5ad').to_df()).astype('float32')
    x2 = np.array(sc.read_h5ad('../datasets/10X-Multiome-Pbmc10k/10X-Multiome-Pbmc10k-RNA.h5ad').to_df()).astype('float32')
    y = None'''

    # gene filter
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

    model = scMultiCluster(input_dim1=input_size1,input_dim2=input_size2, 
                           alpha=args.alpha,beta=args.beta,gama=args.gama,device=args.device).to(args.device)
    print(str(model))

    if not os.path.exists(args.save_dir):
        print("Create save_dir")
        os.makedirs(args.save_dir)

    t0 = time()
    pretrain_latent = model.pretrain_autoencoder(
                        X1=adata1.X, X2=adata2.X, X1_raw=adata1.raw.X, X2_raw=adata2.raw.X, 
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