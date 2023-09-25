import anndata as ad
import numpy as np
import pandas as pd
from scipy import io
import gc

# read the data with backed mode
adata = ad.read_h5ad("data/SEAAD_MTG_RNAseq_finalnuclei_20230505.h5ad", backed = 'r+')
# 1378211 cells and 36601 genes sparse matrix
# Compressed Sparse Row format

## data exploration
adata
#adata.X
# samples
adata.obs
obs_df = pd.DataFrame(adata.obs)
# obs_df.to_csv("data/SEA_AD_MTG_obs.csv", index=False, header=False) #1378211 rows x 134 columns
del obs_df
obs_names = adata.obs.keys()
obs_names = pd.DataFrame(obs_names) #134 
# obs_names.to_csv("data/SEA_AD_MTG_obsnames.csv", index=False, header=False)
del obs_names

#genes
adata.var
var_df = pd.DataFrame(adata.var)
# var_df.to_csv("data/SEA_AD_MTG_var_geneIDs.csv", index=False, header=False) # 36601 rows x 1 columns
del var_df
gc.collect()

# conditions
adata.obs['Braak']
adata[adata.obs.Braak == "Braak VI"] # 185466 × 36601
adata[adata.obs.Braak == "Reference"] # 137303 × 36601
# group by condition
groups = adata.obs.groupby("Braak").indices

# extract cells from 'Reference' and 'Braak VI' samples
braak6_ref = ad.concat([adata[inds] for key, inds in groups.items() if key in ['Reference', 'Braak VI']], merge="same")
# currently SeuratDisk Convert and LoadH5Seurat don't work on all the h5ad formats I produced 
# braak6_ref.write_h5ad("data/SEAAD_MTG_braak6.h5ad")
io.mmwrite("data/mtg_braak6/matrix.mtx", braak6_ref.layers['UMIs'].T)
braak6_ref.obs_names
with open('data/mtg_braak6/barcodes.tsv', 'w') as f:
  for item in braak6_ref.obs_names:
    f.write(item + '\n')
braak6_ref.var_names
with open('data/mtg_braak6/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in braak6_ref.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_braak6/matrix', braak6_ref.X.T)
!ls data/mtg_braak6/
!gzip data/mtg_braak6/*
!ls data/mtg_braak6/
columns_to_keep = ['sample_id','sample_name', 'Subclass','Subclass confidence','Donor ID',  'Braak']
braak6_ref.obs = braak6_ref.obs[columns_to_keep]
braak6_ref.write_csvs("data/mtg_braak6/metadata")
del braak6_ref
gc.collect()


# extract cells from 'Reference' samples
ref = adata[adata.obs.Braak == "Reference"].copy() 
io.mmwrite("data/mtg_ref/matrix.mtx", ref.layers['UMIs'].T)
ref.obs_names
with open('data/mtg_ref/barcodes.tsv', 'w') as f:
  for item in ref.obs_names:
    f.write(item + '\n')
ref.var_names
with open('data/mtg_ref/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in ref.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_ref/matrix', ref.X.T)
!ls data/mtg_ref/
!gzip data/mtg_ref/*
!ls data/mtg_ref/
ref.obs = ref.obs[columns_to_keep]
ref.write_csvs("data/mtg_ref/metadata"
del ref
gc.collect()

# extract cells from 'Braak V' samples
braak5 = adata[adata.obs.Braak == 'Braak V'].copy() 
# braak5 = ad.concat([adata[inds] for key, inds in groups.items() if key in ['Reference', 'Braak V']], merge="same")
# braak5.write_h5ad("data/SEAAD_MTG_braak5.h5ad")
!mkdir data/mtg_braak5/
io.mmwrite("data/mtg_braak5/matrix.mtx", braak5.layers['UMIs'].T)
braak5.obs_names
with open('data/mtg_braak5/barcodes.tsv', 'w') as f:
  for item in braak5.obs_names:
    f.write(item + '\n')
    
braak5.var_names
with open('data/mtg_braak5/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in braak5.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_braak5/matrix', braak5.X.T)
!ls data/mtg_braak5/
!gzip data/mtg_braak5/*
!ls data/mtg_braak5/
braak5.obs = braak5.obs[columns_to_keep]
braak5.write_csvs("data/mtg_braak5/metadata"
del braak5
gc.collect()

# Extract cells 'Braak IV' samples
braak4 = adata[adata.obs.Braak == 'Braak IV'].copy() 
# braak4 = ad.concat([adata[inds] for key, inds in groups.items() if key in ['Reference', 'Braak IV']], merge="same")
# braak4.write_h5ad("data/SEAAD_MTG_braak4.h5ad")
io.mmwrite("data/mtg_braak4/matrix.mtx", braak4.layers['UMIs'].T)
braak4.obs_names
with open('data/mtg_braak4/barcodes.tsv', 'w') as f:
  for item in braak4.obs_names:
    f.write(item + '\n')
    
braak4.var_names
with open('data/mtg_braak4/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in braak4.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_braak4/matrix', braak4.X.T)
!ls data/mtg_braak4/
!gzip data/mtg_braak4/*
!ls data/mtg_braak4/
braak4.obs = braak4.obs[columns_to_keep]
braak4.write_csvs("data/mtg_braak4/metadata"
del braak4
gc.collect()

# Extract cells from 'Braak III' samples
braak3 = adata[adata.obs.Braak == 'Braak III'].copy()
# braak3 = ad.concat([adata[inds] for key, inds in groups.items() if key in ['Reference', 'Braak III']], merge="same")
# braak3.write_h5ad("data/SEAAD_MTG_braak3.h5ad")
io.mmwrite("data/mtg_braak3/matrix.mtx", braak3.layers['UMIs'].T)
braak3.obs_names
with open('data/mtg_braak3/barcodes.tsv', 'w') as f:
  for item in braak3.obs_names:
    f.write(item + '\n')
    
braak3.var_names
with open('data/mtg_braak3/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in braak3.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_braak3/matrix', braak3.X.T)
!ls data/mtg_braak3/
!gzip data/mtg_braak3/*
!ls data/mtg_braak3/
braak3.obs = braak3.obs[columns_to_keep]
braak3.write_csvs("data/mtg_braak3/metadata"
del braak3
gc.collect()

# Extract cells from 'Braak II' samples
braak2 = adata[adata.obs.Braak == 'Braak II'].copy()
# braak2 = ad.concat([adata[inds] for key, inds in groups.items() if key in ['Reference', 'Braak II']], merge="same")
# braak2.write_h5ad("data/SEAAD_MTG_braak2.h5ad")
io.mmwrite("data/mtg_braak2/matrix.mtx", braak2.layers['UMIs'].T)
braak2.obs_names
with open('data/mtg_braak2/barcodes.tsv', 'w') as f:
  for item in braak2.obs_names:
    f.write(item + '\n')
    
braak2.var_names
with open('data/mtg_braak2/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in braak2.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_braak2/matrix', braak2.X.T)
!ls data/mtg_braak2/
!gzip data/mtg_braak2/*
!ls data/mtg_braak2/
braak2.obs = braak2.obs[columns_to_keep]
braak2.write_csvs("data/mtg_braak2/metadata"
del braak2
gc.collect()

# Extract cells from 'Braak I' samples
braak1 = adata[adata.obs.Braak == 'Braak I'].copy()
# braak1 = ad.concat([adata[inds] for key, inds in groups.items() if key in ['Reference', 'Braak I']], merge="same")
# braak1.write_h5ad("data/SEAAD_MTG_braak1.h5ad")
io.mmwrite("data/mtg_braak1/matrix.mtx", braak1.layers['UMIs'].T)
braak1.obs_names
with open('data/mtg_braak1/barcodes.tsv', 'w') as f:
  for item in braak1.obs_names:
    f.write(item + '\n')
    
braak1.var_names
with open('data/mtg_braak1/features.tsv', 'w') as f:
  for item in ['\t'.join([x,x,'Gene Expression']) for x in braak1.var_names]:
    f.write(item + '\n')
io.mmwrite('data/mtg_braak1/matrix', braak1.X.T)
!ls data/mtg_braak1/
!gzip data/mtg_braak1/*
!ls data/mtg_braak1/
braak1.obs = braak1.obs[columns_to_keep]
braak1.write_csvs("data/mtg_braak1/metadata"
del braak1
gc.collect()
