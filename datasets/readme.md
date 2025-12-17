# Example real datasets

Example datasets used in the article can be downloaded from [https://doi.org/10.6084/m9.figshare.30164773.v1](https://doi.org/10.6084/m9.figshare.30899447)

Required objects in h5/h5ad file for running scMAGCA
1) X1: ADT/ATAC count matrix
2) X2: mRNA count matrix
3) Y: True labels (if exist)
4) Batch: batch indicator (for multi-batch analysis)

Other objects in the h5/h5ad files:
1) ADT: feature names in ADT count matirx (only in RNA+ADT data)
2) GenesFromPeaks: feature names in the gene-to-cell matrix mapped from scATAC-seq (only in RNA+ATAC data)
3) Genes: feature names in mRNA count matrix
4) Cell types: cell type of each cell (if exist)
5) Barcodes: cell barcodes (if exits)

