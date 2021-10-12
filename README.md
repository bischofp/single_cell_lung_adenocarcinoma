# Single-cell RNA sequencing of lung adenocarcinomas
Here you find complementary code used for our paper "Single-cell RNA sequencing reveals distinct tumor microenvironmental patterns in lung adenocarcinoma".
The main code used to generate the main figures can be found here: https://codeocean.com/capsule/8321305/tree/v1.

## Infer copy number aberrations
1) Run the preprocessing and cell type annotation steps as described here: https://codeocean.com/capsule/8321305/tree/v1
2) Use the seurat object "epi_anno" and apply the inferCNV algorithm as described here: infercnv_nsclc.R
3) Proceed with further epithelial analyses as described here: https://codeocean.com/capsule/8321305/tree/v1

## Validation using the dataset from Kim et al. 2020
1) Gene expression and clinical data can be downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907
2) Data was preprocessed as described here: lung_kim_et_al_preprocessing.Rmd
3) Cell type labels were transfered from our dataset to the Kim et al. dataset as described here: lung_ingest.ipynb
4) Cell type labels were added to Kim et al. seurat objects as described here: lung_kim_et_al_preprocessing.Rmd
5) For further analyses of stromal and immune cells and correlation analyses of different cell types, see: https://codeocean.com/capsule/8321305/tree/v1

## CellPhoneDB
1) Run the preprocessing and cell type annotation steps as described here: https://codeocean.com/capsule/8321305/tree/v1
2) Fetch and normalize gene expression data as described here: lung_CellPhoneDB_normalization.Rmd
3) Run CellPhoneDB with default parameters and generate figures as described here: lung_CellPhoneDB.R
