# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse)
library(openxlsx)
library(Seurat) 
library(garnett) 
library(org.Hs.eg.db)

# ── Prepare data ────────────────────────────────────────────────────────────── 

sc_tot <- readRDS("mIDH total (CCA integrated Seurat 4.3.0).RDS") 
sc_tot@meta.data %>% dplyr::count(ann_basic)
sc_mye <- subset(sc_tot, subset = ann_basic %in% c("Myeloid"))

# Perform integration 

sc_list <- SplitObject(sc_mye, split.by = "DataSet") 

for (i in 1:length(sc_list)) {
  DefaultAssay(sc_list[[i]]) <- "RNA" 
  sc_list[[i]] <- SCTransform(sc_list[[i]], assay = "RNA", verbose = FALSE)
}

feature <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = 3000)
sc_list <- PrepSCTIntegration(object.list = sc_list, anchor.features = feature)

anchors <- FindIntegrationAnchors(object.list = sc_list,
                                  normalization.method = "SCT",
                                  anchor.features = feature,
                                  reduction = "cca",
                                  dims = 1:10) 

sc_imye <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(sc_imye) <- "integrated"

sc_imye <- sc_imye %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

saveRDS(sc_imye, file = "mIDH myeloid (CCA integrated Seurat 4.3.0).RDS") 

# Plot metadata

p_clust <- DimPlot(sc_imye, reduction = 'umap', group.by = 'seurat_clusters', label = T)
p_immun <- DimPlot(sc_imye, reduction = 'umap', group.by = 'subcat.immunecells', label = T)
p_types <- DimPlot(sc_imye, reduction = 'umap', group.by = 'Types', label = T)
p_dataa <- DimPlot(sc_imye, reduction = 'umap', group.by = 'DataSet', label = T)

p_basic <- p_clust + p_immun + p_types + p_dataa
p_basic
ggsave(p_basic, filename = "1. Myeloid classification metadata (Basic).pdf", width = 10, height = 7, units = c("in"))

# ── Build and train own classifier ────────────────────────────────────────────

# Create CDS object

sc_imye <- readRDS("mIDH myeloid (CCA integrated Seurat 4.3.0).RDS")

expression_matx <- sc_imye@assays$RNA@counts
pheno_metadatas <- sc_imye@meta.data

if (all.equal(colnames(expression_matx), rownames(pheno_metadatas))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matx), nrow(pheno_metadatas)))
  print("If the counts are equal, sort differences will throw this error")
} 

gene_annotation <- data.frame(gene_short_name = rownames(sc_imye@assays$RNA), row.names = rownames(sc_imye@assays$RNA))

if (all.equal(rownames(expression_matx), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matx), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

pheno_metadatas <- new("AnnotatedDataFrame", data = pheno_metadatas)
gene_annotation <- new("AnnotatedDataFrame", data = gene_annotation)

cds.mye_obj <- newCellDataSet(as(expression_matx, "dgCMatrix"),
                              phenoData = pheno_metadatas, 
                              featureData = gene_annotation) 

cds.mye_obj <- estimateSizeFactors(cds.mye_obj)

# Constructing a marker file 

set.seed(123) 
marker_direc <- "C:\\01. Research\\01. mIDH (Hypoxia)\\02. Analysis\\02-01. scRNA-seq (Tirosh)\\06. Annotation (myeloid)\\Abdelfattah_natcomm_myeloid_markers.txt"
marker_check <- check_markers(cds.mye_obj, 
                              marker_direc,
                              db = org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

p_mark <- plot_markers(marker_check)
p_mark
ggsave(p_mark, filename = "2. Myeloid subtypes marker genes (Garnett).pdf", width = 5, height = 50, units = c("in"), limitsize = FALSE)

# Train the classifier and get 

mye_classifier <- train_cell_classifier(cds = cds.mye_obj,
                                        marker_file = marker_direc,
                                        db = org.Hs.eg.db,
                                        cds_gene_id_type = "SYMBOL",
                                        num_unknown = 5,
                                        marker_file_gene_id_type = "SYMBOL", 
                                        classifier_gene_id_type = "SYMBOL") 

saveRDS(mye_classifier, file = "mIDH myeloid Garnett classifier.RDS")

feature_geness <- get_feature_genes(mye_classifier, 
                                    node = "root",
                                    db = org.Hs.eg.db)
head(feature_geness) 

# Classify cells 

cds.mye_obj <- classify_cells(cds = cds.mye_obj, 
                              classifier = mye_classifier, 
                              db = org.Hs.eg.db, 
                              cluster_extend = TRUE, 
                              cds_gene_id_type = "SYMBOL")

saveRDS(cds.mye_obj, file = "mIDH myeloid CDS.RDS")

# ── Fix cell types ──────────────────────────────────────────────────────────── 

# Prioritize annotation from original study first, then label unknown cells according to garnett classifier 

meta_cds <- pData(cds.mye_obj) %>% as.data.frame()
meta_cds$subcat.immunecells[is.na(meta_cds$subcat.immunecells)] <- "Unknown"
meta_cds %>% dplyr::count(subcat.immunecells)

# Broad classification

meta_cds <- meta_cds %>% mutate(ann_mye_details = ifelse(subcat.immunecells %in% "Unknown", meta_cds$cell_type, subcat.immunecells))
meta_cds %>% dplyr::count(ann_mye_details)

meta_cds$ann_mye_details <- factor(meta_cds$ann_mye_details,
                                   levels = c('MC01_i_microglia', 'MC02_h_microglia', 'MC03_s_mac1', 'MC04_MDSC', 'MC05_s_mac2', 'MC06_AP_microglia', 'MC07_a_microglia', 'MC08_DCs', 'MC09_Proliferating', 'Proliferating'), 
                                   labels = c('i-microglia', 'h-microglia', 's-mac 1', 'MDSC', 's-mac 2', 'AP-microglia', 'a-microglia', 'DCs', 'Proliferating mac', 'Proliferating mac'))

meta_cds %>% dplyr::count(ann_mye_details)

# Short classification

meta_cds$ann_mye_general <- factor(meta_cds$ann_mye_details,
                                   levels = c('i-microglia','h-microglia','s-mac 1','s-mac 2','AP-microglia','a-microglia','Proliferating mac','DCs','MDSC'), 
                                   labels = c('Microglia','Microglia','Macrophage','Macrophage','Microglia','Microglia','Macrophage','DCs','MDSC')) 

meta_cds$ann_mye_general[is.na(meta_cds$ann_mye_general)] <- 'Unknown'
meta_cds %>% dplyr::count(ann_mye_general)

# Add classifications to seurat object

sc_imye <- AddMetaData(sc_imye, 
                       meta_cds[, c("cell_type","ann_mye_general","ann_mye_details")], c("garnett_cell_type","ann_mye_general","ann_mye_details")) 

p_garnt <- DimPlot(sc_imye, group.by = c("garnett_cell_type"), reduction = "umap")
p_broad <- DimPlot(sc_imye, group.by = c("ann_mye_general"), reduction = "umap")
p_short <- DimPlot(sc_imye, group.by = c("ann_mye_details"), reduction = "umap")

p_class <- p_garnt + p_broad + p_short
p_class
ggsave(p_class, filename = "3. Myeloid classification metadata (Garnett).pdf", width = 16, height = 7.25, units = c("in"))

# Save myeloid seurat object

saveRDS(sc_imye, file = "mIDH myeloid (CCA integrated Seurat 4.3.0).RDS") 

# ── Include classifications in main seurat ────────────────────────────────────

meta.tot <- sc_tot@meta.data
meta.mye <- sc_imye@meta.data

cell_bar <- rownames(meta.mye)

meta.tot <- meta.tot %>% 
  mutate(garnett_cell_type = ifelse(rownames(meta.tot) %in% cell_bar, meta.mye$garnett_cell_type[match(rownames(meta.tot), cell_bar)], immunecells), 
         ann_mye_general = ifelse(rownames(meta.tot) %in% cell_bar, meta.mye$ann_mye_general[match(rownames(meta.tot), cell_bar)], immunecells), 
         ann_mye_details = ifelse(rownames(meta.tot) %in% cell_bar, meta.mye$ann_mye_details[match(rownames(meta.tot), cell_bar)], immunecells))

# Add metadata to seurat object 

assign_metadata <- function(seurat_obj, meta.tot) {
  res_long_columns <- colnames(meta.tot) 
  for (col in res_long_columns) {
    seurat_obj@meta.data[[col]] <- meta.tot[[col]]
  }
  return(seurat_obj)
}

sc_tot <- assign_metadata(seurat_obj = sc_tot, meta.tot = meta.tot)
saveRDS(sc_tot, file = "mIDH total (CCA integrated Seurat 4.3.0).RDS") 
