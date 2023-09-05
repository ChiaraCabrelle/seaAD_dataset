#setwd("Documents/postdoc_GiorgiLab/scripps_Sanna/waterMaze_scRNAseq/watermaze/seaAD_comparison/")
setwd( "C:/Users/ccabrelle/OneDrive - Scripps Research/Documents/seaAD_comparison/")
library(DESeq2)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(magick)
library(aws.s3)

# Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD)
# https://registry.opendata.aws/allen-sea-ad-atlas/
# Single cell profiling (transcriptomics and epigenomics) data files are in  a public bucket:
# Amazon Resource Name (ARN): arn:aws:s3:::sea-ad-single-cell-profiling
# AWS region: us-west-2

# Finding the bucket
bucket_exists(
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2"
)

# Listing bucket contents
sea_trans_epi_file <- get_bucket_df(
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2", 
  #max = 20000
  ) %>% 
  as_tibble()

save_object(
  object = "MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.csv",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/Reference_MTG_RNAseq_all-nuclei.2022-06-07.csv"
)
read_csv("data/Reference_MTG_RNAseq_all-nuclei.2022-06-07.csv", n_max = 100) 
# samples in rows, genes in columns

#another file that might be useful to download
#MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds

save_object(
  object = "MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
)

SeuratDisk::Convert("C:/Users/ccabrelle/OneDrive - Scripps Research/Documents/seaAD_comparison/data/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad", dest = "Reference_MTG_RNAseq_final-nuclei.h5seurat" ,overwrite = FALSE)

save_object(
  object = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad"
)


### dorsolateral prefrontal cortex samples
save_object(
  object = "DLPFC/RNAseq/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad"
)

setwd( "C:/Users/ccabrelle/OneDrive - Scripps Research/Documents/seaAD_comparison/data/")
library(Seurat)
library(SeuratDisk)
SeuratDisk::Convert(source = "SEAAD_MTG_RNAseq_final-nuclei_2023-05-05.h5ad", dest = "ciao.h5seurat", overwrite = FALSE)
ciao <- LoadH5Seurat("ciao.h5seurat")
install.packages("anndata")
library(anndata)
sea = anndata::read_h5ad("SEAAD_MTG_RNAseq_final-nuclei_2023-05-05.h5ad")
library(reticulate)
ad <- import("anndata", convert = F)
sea_ad <- ad$read_h5ad("SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad")
sea <- Convert(sea_ad, to = "seurat")

setwd("C:/Users/ccabrelle/Documents/seaAD_dataset/data/")
SeuratDisk::Convert("1027251029-raw_feature_bc_matrix.h5", dest = "h5seurat", overwrite = FALSE)
seu_1027251029 <- LoadH5Seurat("1027251029-raw_feature_bc_matrix.h5seurat")

#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
#>  2 dimensional reductions calculated: pca, umap
# S3 method for H5File
Convert(
  source = "1027251029-raw_feature_bc_matrix.h5",
  dest = "h5seurat",
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

install.packages("BiocManager")
BiocManager::install("rhdf5")
library(rhdf5)
h5createFile("myhdf5file.h5")
H5Fopen("")
library(SCP)
library(reticulate)
sc <- import("scanpy")


setwd("C:/Users/ccabrelle/Desktop/")
library(aws.s3)
library(tibble)
library(SeuratDisk)
# Finding the bucket
bucket_exists(
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2"
)

# Listing bucket contents
sea_trans_epi_file <- get_bucket_df(
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2", 
  #max = 20000
) %>% 
  as_tibble()

save_object(
  object = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "SEAAD_MTG_RNAseq_finalnuclei_2023.h5ad"
)

SeuratDisk::Convert("SEAAD_MTG_RNAseq_finalnuclei_2023.h5ad", dest = "SEAAD_MTG_RNAseq_finalnuclei_2023.h5seurat", assay ="RNA" ,overwrite = FALSE)

ad <- anndata::read_h5ad('Results/celltype_assigned_raw.h5ad')
sceasy::convertFormat(ad, from="anndata", to="seurat", outFile='file.rds')

install.packages("renv")
library(renv)
renv::init()
renv::install("reticulate")
renv::use_python() 
# "~/.virtualenvs/r-reticulate/Scripts/python.exe"
#1 virtualenv r-reticulate
#install R packages
pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scater",
  "multtest"
)

# If you are using an {renv} environment
renv::install(pkgs)

# Otherwise do it the normal way
#install.packages(pkgs)

# Install Bioconductor packages
BiocManager::install(bioc_pkgs, update = FALSE)

py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)

reticulate::py_install(py_pkgs)
renv::snapshot()

suppressPackageStartupMessages({
  library("reticulate")
  library("ggplot2")
  library("SingleCellExperiment")
  library("scater")
  library("Seurat")
})
sc <- import("scanpy")
adata <- sc$read_h5ad("SEAAD_MTG_RNAseq_finalnuclei_2023.h5ad")
adata <- sc$read_h5ad("C:/Users/ccabrelle/OneDrive - Scripps Research/Documents/seaAD_comparison/data/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
srt <- scp::adata_to_srt(adata)
srt

head(adata$obs)
head(adata$var)
adata$X[1:5, 1:5]
ggplot(adata$obs, aes(x = n_counts, y = n_genes, colour = louvain)) +
  geom_point()
py_run_string("draw_graph = r.adata.uns['draw_graph']")
py$draw_graph
