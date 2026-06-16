# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(Seurat) 
library(dplyr) 
library(tidyverse) 
library(RColorBrewer) 
library(viridis) 
library(ggthemes) 
library(ggplot2) 
library(cowplot) 
library(data.table) 
library(gridExtra) 

# ── Data mining & creating seurat object ────────────────────────────────────── 

# Read in the feature-barcode matrices 

samples = c("LGG-04-1" , "LGG-04-2" , "LGG-04-3", "LGG-03", "SF11136" , "SF12017" , "SF11949") 
IDHmut.data = list()

IDHmut.data[[1]] <- ReadMtx(mtx = "GSM5518630_LGG-04-1_matrix.mtx.gz", cells = "GSM5518630_LGG-04-1_barcodes.tsv.gz", features = "GSM5518630_LGG-04-1_features.tsv.gz", feature.column = 1) 
IDHmut.data[[2]] <- ReadMtx(mtx = "GSM5518631_LGG-04-2_matrix.mtx.gz", cells = "GSM5518631_LGG-04-2_barcodes.tsv.gz", features = "GSM5518631_LGG-04-2_features.tsv.gz", feature.column = 1) 
IDHmut.data[[3]] <- ReadMtx(mtx = "GSM5518632_LGG-04-3_matrix.mtx.gz", cells = "GSM5518632_LGG-04-3_barcodes.tsv.gz", features = "GSM5518632_LGG-04-3_features.tsv.gz", feature.column = 1) 
IDHmut.data[[4]] <- ReadMtx(mtx = "GSM5518638_LGG-03_matrix.mtx.gz", cells = "GSM5518638_LGG-03_barcodes.tsv.gz", features = "GSM5518638_LGG-03_features.tsv.gz", feature.column = 1) 
IDHmut.data[[5]] <- ReadMtx(mtx = "GSM4119535_SF11136_matrix.mtx.gz", cells = "GSM4119535_SF11136_barcodes.tsv.gz", features = "GSM4119535_SF11136_features.tsv.gz", feature.column = 1) 
IDHmut.data[[6]] <- ReadMtx(mtx = "GSM4119536_SF12017_matrix.mtx.gz", cells = "GSM4119536_SF12017_barcodes.tsv.gz", features = "GSM4119536_SF12017_features.tsv.gz", feature.column = 1) 
IDHmut.data[[7]] <- ReadMtx(mtx = "GSM4119538_SF11949_matrix.mtx.gz", cells = "GSM4119538_SF11949_barcodes.tsv.gz", features = "GSM4119538_SF11949_features.tsv.gz", feature.column = 1) 

# Convert each feature-barcode matrix to a seurat object 

scrna.list = list() 
scrna.list[[1]] = CreateSeuratObject(counts = IDHmut.data[[1]], min.cells=100, min.features=700, project = samples[1]) 
scrna.list[[1]][["DataSet"]] = samples[1] 
scrna.list[[2]] = CreateSeuratObject(counts = IDHmut.data[[2]], min.cells=100, min.features=700, project = samples[2]) 
scrna.list[[2]][["DataSet"]] = samples[2] 
scrna.list[[3]] = CreateSeuratObject(counts = IDHmut.data[[3]], min.cells=100, min.features=700, project = samples[3]) 
scrna.list[[3]][["DataSet"]] = samples[3] 
scrna.list[[4]] = CreateSeuratObject(counts = IDHmut.data[[4]], min.cells=100, min.features=700, project = samples[4]) 
scrna.list[[4]][["DataSet"]] = samples[4] 
scrna.list[[5]] = CreateSeuratObject(counts = IDHmut.data[[5]], min.cells=100, min.features=700, project = samples[5]) 
scrna.list[[5]][["DataSet"]] = samples[5] 
scrna.list[[6]] = CreateSeuratObject(counts = IDHmut.data[[6]], min.cells=100, min.features=700, project = samples[6]) 
scrna.list[[6]][["DataSet"]] = samples[6] 
scrna.list[[7]] = CreateSeuratObject(counts = IDHmut.data[[7]], min.cells=100, min.features=700, project = samples[7]) 
scrna.list[[7]][["DataSet"]] = samples[7] 

# Merge the Seurat objects into a single object 

scrna <- merge(x = scrna.list[[1]], y = c(scrna.list[[2]], scrna.list[[3]], scrna.list[[4]], scrna.list[[5]], scrna.list[[6]], scrna.list[[7]]), add.cell.ids = c("LGG-04-1" , "LGG-04-2" , "LGG-04-3", "LGG-03", "SF11136" , "SF12017", "SF11949"), project = "mIDH") 

# ── Add metadata ────────────────────────────────────────────────────────────── 

# Classify into astrocytoma and oligodendroglioma based on metadata 

datasets <- FetchData(object = scrna, vars = "DataSet") 
datasets %>% dplyr::count(DataSet) 
datasets <- datasets %>% dplyr::mutate(Types = ifelse(DataSet %in% c("SF11136", "SF12017", "LGG-03"), "Astrocytoma", "Oligodendroglioma"))
datasets <- datasets %>% dplyr::mutate(Grade = ifelse(DataSet %in% c("SF11949"), "Grade_3", "Grade_2"))

scrna <- AddMetaData(object = scrna, metadata = datasets$Types, col.name = "Types") 
scrna <- AddMetaData(object = scrna, metadata = datasets$Grade, col.name = "Grade") 

# Add cell type labels based on the original studies 

metadata_org <- read.csv("Single-cell metadata.csv", row.names = 1)
metadata_sub <- subset(metadata_org, rownames(metadata_org) %in% colnames(scrna)) 
all(rownames(metadata_sub) %in% colnames(scrna)) 
all(colnames(scrna) %in% rownames(metadata_sub)) 

scrna <- AddMetaData(object = scrna, metadata = metadata_sub$celltypes, col.name = 'immunecells') 
scrna <- AddMetaData(object = scrna, metadata = metadata_sub$Sex, col.name = 'Sex') 

# ── Quality control ─────────────────────────────────────────────────────────── 

head(scrna@meta.data) 

# Mitochondrial and ribosomal percentages 

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna), value = TRUE) 
percent.rb <- Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'));
scrna[['percent.rb']] <- percent.rb 
scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-") 

# nCount RNA 

vln_nCount <- VlnPlot(object = scrna, features = "nCount_RNA", pt.size=0, group.by = "DataSet") 
vln_nCount 
ggsave("01. QC metrics (nCount).pdf", plot = vln_nCount, width = 5, height = 3, units = "in")

# nFeature RNA 

vln_nFeatr <- VlnPlot(object = scrna, features = "nFeature_RNA", pt.size=0, group.by = "DataSet") 
vln_nFeatr
ggsave("01. QC metrics (nFeature).pdf", plot = vln_nFeatr, width = 5, height = 3, units = "in")

# Feature-feature relationships (scatter plot) 

sca_ctmt <- FeatureScatter(scrna, feature1 = "nCount_RNA", feature2 = "percent.mt")
sca_ctft <- FeatureScatter(scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sca_comb <- sca_ctmt + sca_ctft 
ggsave("02. F-F scatters.pdf", plot = sca_comb, width = 6, height = 4, units = "in")

# Subset the main data according to QC analysis 

scrna <- subset(scrna, subset = nFeature_RNA > 700 & nCount_RNA < 25000 & percent.mt < 10 & percent.rb < 40)

# ── Integration ───────────────────────────────────────────────────────────────  

# Split the datasets into a list 

sc.list <- SplitObject(scrna, split.by = "DataSet")

# Normalize and identify variable features for each dataset independently

sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1500)
}) 

# Select features and perform integration 

i.feature <- SelectIntegrationFeatures(object.list = sc.list) 
i.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = i.feature)
i.combine <- IntegrateData(anchorset = i.anchors) 

# original unmodified data still resides in the 'RNA' assay

DefaultAssay(i.combine) <- "integrated"

# Run the standard workflow for visualization and clustering

i.combine <- ScaleData(i.combine, verbose = FALSE) %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) 

saveRDS(i.combine, file = "mIDH (CCA integrated Seurat 4.3.0).RDS")

# ── Plot everything ─────────────────────────────────────────────────────────── 

# Samples  

umap_samples <- DimPlot(object = i.combine, reduction = "umap", group.by = "DataSet", pt.size = 0.1, label = TRUE, repel = TRUE)
umap_samples 
tsne_samples <- DimPlot(object = i.combine, reduction = "tsne", group.by = "DataSet", pt.size = 0.1, label = TRUE, repel = TRUE)
tsne_samples
ggsave("03. Samples Umap.pdf", plot = umap_samples, width = 4.7, height = 3.5, units = "in")
ggsave("03. Samples tSNE.pdf", plot = tsne_samples, width = 4.7, height = 3.5, units = "in")

# Cell categories 

umap_cellpop <- DimPlot(object = i.combine, reduction = "umap", group.by = "immunecells", pt.size = 0.1, label = TRUE, repel = TRUE)
umap_cellpop 
tsne_cellpop <- DimPlot(object = i.combine, reduction = "tsne", group.by = "immunecells", pt.size = 0.1, label = TRUE, repel = TRUE)
tsne_cellpop
ggsave("04. Immune cells Umap.pdf", plot = umap_cellpop, width = 4, height = 3, units = "in")
ggsave("04. Immune cells tSNE.pdf", plot = tsne_cellpop, width = 4, height = 3, units = "in")

# Seurat clsuters  

umap_seuclus <- DimPlot(object = i.combine, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.2, label = TRUE, repel = TRUE) 
umap_seuclus 
tsne_seuclus <- DimPlot(object = i.combine, reduction = "tsne", group.by = "seurat_clusters", pt.size = 1.2, label = TRUE, repel = TRUE) 
tsne_seuclus
ggsave("05. Seurat clusters Umap.pdf", plot = umap_seuclus, width = 5.6, height = 5, units = "in") 
ggsave("05. Seurat clusters tSNE.pdf", plot = tsne_seuclus, width = 4.6, height = 4, units = "in") 

