# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(tidyr) 
library(openxlsx) 
library(Seurat) 
library(patchwork)
library(NMF)
library(pheatmap) 
library(RColorBrewer)
library(BiocGenerics) 

# ── Prepare data ────────────────────────────────────────────────────────────── 

# Re-integrate and re-cluster glioma population 

sc_object <- readRDS("mIDH (CCA integrated Seurat 4.3.0).RDS") 
sc_glioma <- subset(sc_object, subset = immunecells == "Glioma") 

# Remove smaller ones and merge if they exist 

DefaultAssay(sc_glioma) <- "RNA"

sc_list <- SplitObject(sc_glioma, split.by = "DataSet")
lrge_dt <- sc_list[sapply(sc_list, ncol) >= 30]
smal_dt <- sc_list[sapply(sc_list, ncol) < 30]

if (length(smal_dt) > 1) {
  mergsml <- merge(smal_dt[[1]], y = smal_dt[-1])
  sc_list <- c(lrge_dt, list(mergsml)) 
} else if (length(smal_dt) == 1) {
  sc_list <- c(lrge_dt, smal_dt)
} else {
  sc_list <- lrge_dt
}

# Perform integration 

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

sc_inew <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(sc_inew) <- "integrated"

sc_inew <- sc_inew %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

saveRDS(sc_inew, file = "Glioma-new (CCA integrated Seurat 4.3.0).rds") 

# Plot metadata 

sc_inew@meta.data %>% dplyr::count(DataSet) 

sn_metavars <- c("Types", "Grade", "DataSet", "seurat_clusters") 
type_colors <- c("Astrocytoma" = "#31b5ff", "Oligodendroglioma" = "#ff009e")  
grad_colors <- c("II" = "#BBC863", "III" = "#850E35")
data_colors <- c("LGG-03" = "#F4F754", "LGG-04-1" = "#0BA6DF", "LGG-04-2" = "#6495ED", "LGG-04-3" = "#F75270", "SF11949" = "#FFB300", "SF12017" = "#D62728")

umap_plots <- lapply(sn_metavars, function(group) {
  if(group == "seurat_clusters") {
    DimPlot(sc_inew, reduction = "umap", label = T, repel = T, raster = F) 
  } else {
    cols <- switch(group, "Types" = type_colors, "Grade" = grad_colors, "DataSet" = data_colors, NULL)
    DimPlot(sc_inew, group.by = group, reduction = "umap", label = T, repel = T, raster = F, cols = cols) + ggtitle(group)
  }
}) 

AO_umap <- wrap_plots(umap_plots, ncol = 2) 
AO_umap
ggsave("01. Metadata umaps (only glioma).pdf", plot = AO_umap, width = 8.25, height = 6, units = "in")  

# ── Run NMF ───────────────────────────────────────────────────────────────────  

# Prepare data 

sc_inew <- readRDS("Glioma-new (CCA integrated Seurat 4.3.0).rds")
sc_inew <- ScaleData(sc_inew, assay = "SCT", verbose = FALSE)

allgene <- rownames(sc_inew)
mt_gene <- grep("^MT", rownames(sc_inew@assays$SCT@scale.data), value = TRUE)
rb_gene <- grep("^RPS|^RPL", rownames(sc_inew@assays$SCT@scale.data), value = TRUE)
keep_gn <- setdiff(allgene, c(mt_gene, rb_gene))

# Run sample-wise NMF

nmf_file <- "nmf_new_results/"
dir.create(nmf_file, recursive = TRUE) 

nmf_list <- list() 

for (samp in unique(sc_inew$orig.ident)) {
  cat("Processing sample:", samp, "\n")
  tum <- sc_inew[, sc_inew$orig.ident == samp]
  tum <- ScaleData(tum)
  cts <- tum@assays$SCT@scale.data[keep_gn, ]
  cts[cts < 0] <- 0
  cts <- cts[rowSums(cts) > 0, ]
  nmf_res <- nmf(cts, rank = 2:10, nrun = 10, seed = 123, method = "snmf/r")
  save(nmf_res, cts, file = file.path(nmf_file, paste0("sample_", samp, "_nmf2to10.RData"))) 
  nmf_list[[samp]] <- nmf_res
} 

# ── Identify robust/refined programs ────────────────────────────────────────── 

top_n <- 50
silh.filter <- 0.35

top_gene_list <- list()
basis_vectors <- list()

extract_elements <- function(nmf_res, sample_name, top_n = top_n) {
  if (!"measures" %in% names(nmf_res)) {
    stop(paste("NMF object for", sample_name, "contains no 'measures' slot"))
  }
  m <- nmf_res$measures 
  if (!"silhouette.basis" %in% names(m)) {
    stop("No silhouette.basis found in measures")  # silhouette.basis must exist
  }
  good_k.ft <- m$rank[m$silhouette.basis > silh.filter] 
  fit_lists <- nmf_res$fit
  out_genes <- list()
  out_basis <- list()
  
  for (k in good_k.ft) {
    k_chr <- as.character(k)
    if (!k_chr %in% names(fit_lists)) next
    fit_obj <- fit_lists[[k_chr]]
    if (is.null(fit_obj)) next
    W <- basis(fit_obj) 
    colnames(W) <- paste0(sample_name, "_k", k, "_prog", seq_len(ncol(W)))
    for (p in colnames(W)) {
      out_basis[[p]] <- W[, p]
    }
    top_genes <- apply(W, 2, function(vec) {
      names(sort(vec, decreasing = TRUE)[1:top_n])
    })
    out_genes[[paste0(sample_name, "_k", k)]] <- top_genes
  } 
  return(list(top_genes = out_genes, basis = out_basis))
}

# Run on all samples 

for (samp in names(nmf_list)) {
  cat("Processing:", samp, "\n")
  res <- extract_elements(nmf_list[[samp]], samp, top_n)
  top_gene_list[[samp]] <- res$top_genes
  basis_vectors <- c(basis_vectors, res$basis)
}

length(basis_vectors)
head(names(basis_vectors))

# Build binary program matrix 

all_genes <- keep_gn
all_progs <- names(basis_vectors)

program_matrix <- matrix(0,
                         nrow = length(all_genes),
                         ncol = length(all_progs),
                         dimnames = list(all_genes, all_progs))

for (prog in all_progs) {
  top_g <- names(sort(basis_vectors[[prog]], decreasing = TRUE)[1:top_n])
  program_matrix[top_g, prog] <- 1
}

# ──  Metaprogram clustering programs ──────────────────────────────────────────

# Correlation 

program_corr <- cor(program_matrix, method = "pearson")

# Cluster annotation  

k_meta <- 4
hc_cal <- hclust(as.dist(1 - program_corr), method = "ward.D2")
clustr <- cutree(hc_cal, k = k_meta) 

# Sample annotation  

samp_p <- sub("_k.*", "", all_progs)
ann_df <- data.frame(Sample = factor(samp_p))
rownames(ann_df) <- all_progs

samp_c <- c("LGG-03" = "#228b22", "LGG-04-1" = "#ffa500", "LGG-04-2" = "#ff0000", "LGG-04-3" = "#a020f0", "SF11949" = "#0055aa", "SF12017" = "#e7298a")
ann_cl <- list(Sample = samp_c)

pht_mp <- pheatmap(program_corr,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   cluster_rows = hc_cal,
                   cluster_cols = hc_cal,
                   annotation_colors = ann_cl,
                   annotation_col = data.frame(MetaProgram = factor(clustr), Sample = samp_p),
                   annotation_row = data.frame(MetaProgram = factor(clustr)),
                   show_rownames = FALSE,
                   show_colnames = FALSE) 
pht_mp
ggsave("02. Identified MPs.pdf", plot = pht_mp, width = 4.5, height = 4, units = "in")

# ── Metaprogram specific genesets ───────────────────────────────────────────── 

# Organize programs by cluster (list of program names per meta-program) 

program_cluster_list <- split(names(clustr), clustr) 

top_genes_per_cluster <- list()

for (clust_id in names(program_cluster_list)) {
  progs <- program_cluster_list[[clust_id]] 
  mat <- do.call(cbind, basis_vectors[progs]) # subset basis vectors 
  avg_vec <- rowMeans(mat) 
  top_genes <- names(sort(avg_vec, decreasing = TRUE))[1:top_n]
  top_genes_per_cluster[[paste0("MetaProg", clust_id)]] <- top_genes
}

top_genes_per_cluster

top_genes_df <- bind_rows(lapply(names(top_genes_per_cluster), function(mp) {
  data.frame(MetaProgram = mp, Gene = top_genes_per_cluster[[mp]], stringsAsFactors = FALSE)
  })
)

# Determine maximum number of top genes

max_genes <- max(sapply(top_genes_per_cluster, length))

# Pad each list to the same length

top_genes_wd <- lapply(top_genes_per_cluster, function(g) {
  length(g) <- max_genes
  g
})

top_genes_df <- as.data.frame(top_genes_wd, stringsAsFactors = F)
write.csv(top_genes_df, file = "NMF-MP-4 geneset.csv", row.names = F)

# ── Scoring (not too fancy) ─────────────────────────────────────────────────── 

gs.comp <- list(MetaProg1 = as.character(top_genes_df$MetaProg1),
                MetaProg2 = as.character(top_genes_df$MetaProg2),
                MetaProg3 = as.character(top_genes_df$MetaProg3), 
                MetaProg4 = as.character(top_genes_df$MetaProg4)) 

gs.list <- lapply(gs.comp, function(x) list(as.character(x))) 

DefaultAssay(sc_inew) <- "SCT"

for (nm in names(gs.list)) {
  sc_inew <- AddModuleScore(sc_inew, features = gs.list[[nm]], name = nm) 
  oldname <- grep(paste0("^", nm), colnames(sc_inew@meta.data), value = TRUE)
  colnames(sc_inew@meta.data)[colnames(sc_inew@meta.data) == oldname] <- nm
}

mps_ft <- FeaturePlot(sc_inew, reduction = "umap", features = c("MetaProg1", "MetaProg2", "MetaProg3", "MetaProg4"))  
mps_ft 
ggsave("03. MPs scoring.pdf", plot = mps_ft, width = 8, height = 6, units = "in")

