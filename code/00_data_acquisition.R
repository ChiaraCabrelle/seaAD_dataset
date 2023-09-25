# renv::init()
# renv::status()
library(tibble)
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


# # middle temporal gyrus
# MTG/RNAseq/Supplementary Information/Supplementary Table 5.csv
#
# MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2023-05-05.h5ad
# MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad
# 
# MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds
# MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad
# MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad

save_object(
  object = "MTG/RNAseq/Supplementary Information/Supplementary Table 5.csv",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/MTG_RNAseq_Supplementary_Table5.csv"
)


save_object(
  object = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/SEAAD_MTG_RNAseq_finalnuclei_20230505.h5ad"
)


save_object(
  object = "MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2023-05-05.h5ad",
  bucket = "s3://sea-ad-single-cell-profiling/", 
  region = "us-west-2",
  file = "data/SEAAD_MTG_RNAseq_allnuclei_20230505.h5ad"
)

