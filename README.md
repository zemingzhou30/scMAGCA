# scMAGCA

Single Cell Multi-omics adversarial graph convolutional autoencoder

We develop a multimodal graph convolutional model based on adversarial learning (scMAGCA) to jointly analyze scMulti-omics efficiently. scMAGCA effectively integrates scMulti-omics data and learns the joint representation.  Extensive experiments demonstrate that scMAGCA exhibits superior clustering performance, outperforming existing multi-modal and single-modal clustering techniques across various single-cell multimodal datasets including CITE-seq datasets and SMAGE-seq datasets.

## Table of contents

- [Dependencies](#Dependencies)
- [Usage](#Usage)
- [Output](#Output)
- [Arguments](#Arguments)

## <a name="Dependencies">Dependencies</a>

Python 3.9.7

Pytorch 2.1.0

Pytorch Geometric 2.5.2

Scanpy 1.9.6

Scipy 1.11.3

Sklearn 0.22.1

Numpy 1.22.0

Pandas 1.4.4

All experiments of scMAGCA in this study are conducted on Nvidia 4090 (24G) GPU. We suggest to install the dependencies in a conda environment (conda create -n scMAGCA).

## <a name="Dependencies">Usage</a>

Required objects in h5ad or h5 file for running scMAGCA

1. X1: ADT/ATAC count matrix

2. X2: mRNA count matrix

3. Y: True labels (if exist)

4. Batch: batch indicator (for multi-batch analysis)

Other objects in the h5ad or h5 files:

1. ADT: feature names in ADT count matirx (only in CITE-seq data)

2. GenesFromPeaks: feature names in the gene-to-cell matrix mapped from scATAC-seq (only in SMAGE-seq data)

3. Genes: feature names in mRNA count matrix

4. Cell types: cell type of each cell (if exist).



1. Prepare the input data in h5ad or h5 format (All used datassets are in "datasets" floder and the generated cell-cell graph structure will be saved in the "raw" folder corresponding to the dataset). 

2. Run scMAGCA according to the running script in "scripts" folder (run_scMAGCA.sh for ADT/ATAC+mRNA data and run_scMAGCA_batch.sh for multi-batch data clustering). For example, you can execute scMAGCA on multi-batch dataset "GSE164378"  by shell command: "bash run_scMAGCA_batch.sh".

## <a name="Dependencies">Output</a>

1. scMAGCA outputs latent representation of data which can be used for further downstream analyses and visualized by Umap; 
2. scMAGCA outputs predicted label of data which can be used for calculating related evaluation indexs.
3. Multi-batch scMAGCA outputs a latent representation of integrated datasets on which the batch effects are corrected. 

## <a name="Dependencies">Arguments</a>

Structure: X1(ADT or ATAC), X2(RNA), Y(label, if exit), Batch (Batch indicator for multi-batch data clustering).

--n_cluster: number of clusters (K); scMAGCA will estimate K if this arguments is set to -1.

--data_file: data input and cell-cell graph structure saving path.

--pretrain_epochs: number of epochs for pre-training. Default: 400.

--maxiter: maximum epochs of training. Default: 2000.

--save_dir: the directory to store the outputs. 

--embedding_file: if save embedding file. Default: Yes.

--prediction_file: if save prediction file. Default: Yes.

--ad_out: the output dim of discriminator. Default: 16 for CITE-Seq; 32 for SMAGE-Seq.

--tol: the criterion to stop the model, which is a percentage of changed labels. Default: 0.001.

--resolution: the resolution parameter to estimate K. Default: 0.08.

--filter1: if do feature selection on ADT/ATAC. Default: No.

--filter2: if do feature selection on Genes. Default: No.		

--f1: number of high variable genes from ATAC(in X1) used for clustering if doing the featue selection. Default: 2000.

--f2: number of high variable genes from Genes(in X2) used for clustering if doing the featue selection. Default: 2000.	

--device: training device. Default: cuda.

 *We denote  antibody-derived tags (ADTs) + Gene Expression as CITE-Seq and 10X Single-Cell Multiome ATAC + Gene Expression technology as SMAGE-seq for convenience. 
