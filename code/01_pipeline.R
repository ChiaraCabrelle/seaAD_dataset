library(Seurat)
library(dplyr)
# library(ggplot2)

h5_mtg <-
  ReadMtx(
    'data/mtg_braak6/SEAAD_MTG_braak6_UMIs.mtx.gz',
    features = 'data/mtg_braak6/features.tsv',
    cells = 'data/mtg_braak6/barcodes.tsv'
  )
seuset = CreateSeuratObject(
  h5_mtg,
  project = 'sea_mtg',
  min.cells = 3,
  min.features = 200
) # 35163  genes * 322769 samples
save(seuset, file = 'results/01_seuset.rda')

# Visualize QC metrics as a violin plot
seuset[["percent.mt"]] <- PercentageFeatureSet(seuset, pattern = "^MT-")
VlnPlot(seuset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# ok, data were preprocessed as expected

# load(file = 'results/01_seuset.rda')
obs <- read.csv("data/mtg_braak6/metadata/obs.csv") # 322769      7
# genes <- read.csv('results/mtg_braak6/var.csv')
all(colnames(seuset) == obs$exp_component_name) # TRUE
unique(obs$Braak) # "Reference" "Braak VI"

ref <-
  unique(obs[obs$Braak == 'Reference', ]$Donor.ID) # 5 Reference samples
br6 <-
  unique(obs[obs$Braak == 'Braak VI', ]$Donor.ID) # 15 Braak VI samples
obs <- obs %>%
  mutate(ID_name = case_when(
    Braak == 'Reference' ~ paste0('REF', match(Donor.ID, ref)),
    Braak == 'Braak VI' ~ paste0('BRVI', match(Donor.ID, br6)),
    TRUE ~ Donor.ID  # Keep the original sample_name if none of the conditions match
  ))
rownames(obs) <-
  paste(obs$ID_name, obs$exp_component_name, sep =  '_')
colnames(seuset) <- rownames(obs)

seuset@meta.data$orig.ident <- obs$sample_name
seuset@meta.data$ID_name <- obs$ID_name
seuset@meta.data$condition <- obs$Braak
seuset@meta.data$celltype <- obs$Subclass

# How many reads can be mapped to azimuth
which(is.na(seuset@meta.data$celltype)) # 0
length(unique(seuset@meta.data$celltype)) # 24 cell types
# saveRDS(seuset, 'results/01_seuset_celltype.rds')
save(seuset, file = 'results/01_seuset_Wannotation.rda')
metadata <- seuset@meta.data
save(metadata, file = 'results/01_metadata.rda')

## Sanity check: cells per sample
ref_samples <-subset(seuset, subset = condition == 'Reference')
braak6_samples <-subset(seuset, subset = condition == 'Braak VI')
# distribution in Braak VI and Reference is the same
VlnPlot(ref_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(braak6_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


### 2 - Create a rawcounts matrix by cell type and sample type ----
# seuset <- readRDS('results/01_seuset_celltype.rds')
# obs <- read.csv("results/mtg_braak6/obs.csv") # 322769      7
celltypes <-  setNames(seuset@meta.data$celltype, colnames(seuset))
ID_name <-  setNames(seuset@meta.data$ID_name, colnames(seuset))
condition <- setNames(seuset@meta.data$condition, colnames(seuset))
rawcounts <- seuset@assays$RNA@counts
# There are 20 individuals, 24 cell types
20 * 24 #  480 cellcounts ID_name
combos <- c()
for (s in unique(ID_name)) {
  for (c in unique(celltypes)) {
    new <- paste0(s, "_", c)
    combos <- c(combos, new)
  }
}
cellcounts <- matrix(0, nrow = nrow(rawcounts), ncol = length(combos))
rownames(cellcounts) <- rownames(rawcounts)
colnames(cellcounts) <- combos
metadata
# Let's fill it
design <- matrix(NA, nrow = length(combos), ncol = 3)
colnames(design) <- c("celltype", "sample", "condition")
rownames(design) <- combos
for (s in unique(ID_name)) {
  message("Doing ", s)
  sampcells <- rownames(metadata[metadata$ID_name == s,])
  for (c in unique(celltypes)) {
    sampname <- paste0(s, "_", c)
    design[sampname, ] <- 
      c(c, s, metadata[metadata$ID_name == s,]$condition[1])
    herecells <- rownames(metadata[metadata$celltype == c,])
    herecells <- intersect(sampcells, herecells)
    if (length(herecells) > 1) {
      subcounts <- rawcounts[, herecells]
      sampvector <- apply(subcounts, 1, sum)
      cellcounts[, sampname] <- sampvector
    }
  }
}
save(cellcounts, design, file = "results/01_cellcounts.rda")


### 4 - Differential expression analysis ----
library(DESeq2)
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx48G"))
library(xlsx)

de <- list()
xlsxfile <- "results/01_de.xlsx"
unlink(xlsxfile)
celltypes <- unique(design[, "celltype"])

celltypes <- c("pseudobulk", celltypes)

for (celltype in celltypes) {
  message("Doing ", celltype)
  if (celltype == "pseudobulk") {
    subcounts <- cellcounts
    subdesign <- design
  } else{
    subcounts <- cellcounts[, design[, "celltype"] == celltype]
    subdesign <- design[design[, "celltype"] == celltype, ]
  }
  nonzerogenes <- which(apply(subcounts, 1, sum) != 0)
  nonzerosamples <- which(apply(subcounts, 2, sum) != 0)
  subcounts <- subcounts[nonzerogenes, ]
  subcounts <- subcounts[, nonzerosamples]
  subdesign <- subdesign[nonzerosamples, ]
  sbd_trt <- table(subdesign[, 'condition'])
  print(sbd_trt)
  if (sbd_trt['Braak VI'] > 1 & sbd_trt['Reference'] > 1) {
    dds <- DESeqDataSetFromMatrix(
      countData = subcounts,
      colData = subdesign,
      design =  ~ condition
    )
    dds$condition <- relevel(dds$condition, ref = "Reference")
    dds <- DESeq(dds)
    #
    res <- results(dds, contrast = c("condition", "Braak VI", "Reference"))
    res <- res[order(res$pvalue), ]
    contrast <- paste0("br6VSctrl_", celltype)
    contrastxlsx <- gsub("/", "_", contrast)
    de[[contrast]] <- res
    write.xlsx2(res,
                file = xlsxfile,
                sheetName = contrastxlsx,
                append = TRUE)
  }
}
length(de)
save(de, file = 'results/01_de.rda')

# sanity check
load(file = 'results/01_de.rda')
library(EnhancedVolcano)
EnhancedVolcano(de$br6VSctrl_pseudobulk,
                lab = rownames(de$br6VSctrl_pseudobulk),
                x = 'log2FoldChange',
                y = 'pvalue')

