# Example real datasets

Example datasets used in the article can be downloaded from https://pan.baidu.com/s/12eSdzQLEDml_0nXflYgF1w?pwd=exsk 

Required objects in h5/h5ad file for running scMDC
1) X1: ADT/ATAC count matrix
2) X2: mRNA count matrix
3) Y: True labels (if exist)
4) Batch: batch indicator (for multi-batch analysis)

Other objects in the h5/h5ad files:
1) ADT: feature names in ADT count matirx (only in CITE-seq data)
2) GenesFromPeaks: feature names in the gene-to-cell matrix mapped from scATAC-seq (only in SMAGE-seq data)
3) Genes: feature names in mRNA count matrix
4) Cell types: cell type of each cell (if exist)
5) Barcodes: cell barcodes (if exits)
