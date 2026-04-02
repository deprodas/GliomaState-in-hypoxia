# Author : Depro Das, 3DBM Lab, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(Seurat)
library(openxlsx)
library(tidyverse)
library(decoupleR)
library(ggplot2)

# ── Data input ──────────────────────────────────────────────────────────────── 

sample_dir <- "C:\\01. Research\\01. mIDH (hypoxia)\\02. Analysis\\03. Analysis (stRNA-seq)\\02. Spatial enrichment"
output_dir <- file.path(sample_dir, "Res_dec_2025") 
dir.create(output_dir, showWarnings = FALSE)

strds_files <- list.files(sample_dir, pattern = "\\.RDS$", full.names = TRUE)
seurat_list <- lapply(strds_files, readRDS)
names(seurat_list) <- tools::file_path_sans_ext(basename(strds_files))

idh_geneset <- read.xlsx(file.path(sample_dir, "mIDH programs_depro.xlsx")) %>% 
  pivot_longer(everything(), names_to = "source", values_to = "target") %>% 
  mutate(mor = 1) %>% 
  filter(!is.na(target), target != "") %>% 
  distinct(source, target, .keep_all = TRUE) 

# ── Run enrichment ──────────────────────────────────────────────────────────── 

# Methods to run (if additional methods are included then the features with na/inf values need to filtered out inside the loop)

gsva_methods <- c("gsva", "plage", "ssgsea", "zscore")
dcouple_stat <- c("ulm", "mlm", "wsum", "consensus")
minsize_gsva <- 2L 

# Loop over samples (i. 1st loop for count extraction, ii. 2nd loop for gsva methods, iii. 2nd loop for decoupler statistics)

for (sample in names(seurat_list)) {
  message("Processing sample: ", sample)
  obj <- seurat_list[[sample]]
  obj <- UpdateSeuratObject(obj) 
  sample_dir <- file.path(output_dir, sample)
  dir.create(sample_dir, showWarnings = FALSE)
  
  mat <- obj@assays$Spatial$counts %>% as.matrix() %>% `+`(1) %>% log2()
  
  for (m in gsva_methods) { 
    message("  Running GSVA method: ", m)
    res_gsea <- run_gsva(mat = mat, network = idh_geneset, .source = "source", .target = "target", minsize = minsize_gsva, method = m)
    score_df <- res_gsea %>% 
      select(source, condition, score) %>% 
      pivot_wider(id_cols = source, names_from = condition, values_from = score) %>% 
      column_to_rownames("source") %>% 
      t() %>% 
      as.data.frame() 
    
    obj_m <- obj
    for (gs in colnames(score_df)) {
      obj_m@meta.data[[gs]] <- score_df[[gs]]
    }
    
    p <- SpatialFeaturePlot(obj_m, 
                            features = colnames(score_df), 
                            crop = TRUE,
                            pt.size.factor = 2.5,
                            image.alpha = 0) 
    ggsave(file.path(sample_dir, paste0(sample, "_GSVA_", m, ".pdf")), p, width = 12, height = 12) 
    saveRDS(obj_m, file = file.path(sample_dir, paste0(sample, "_GSVA_", m, ".RDS")))
  }
  
  message("  Running decouple (all statistics)") 
  res_decouple <- decouple(mat = mat, network = idh_geneset, .source = "source", .target = "target", minsize = 0)
  
  for (stat in dcouple_stat) {
    message("Processing statistic: ", stat)
    score_df <- res_decouple %>% 
      filter(statistic == stat) %>% 
      select(source, condition, score) %>% 
      pivot_wider(id_cols = source, names_from = condition, values_from = score) %>% 
      column_to_rownames("source") %>% 
      t() %>% 
      as.data.frame() 
    
    obj_d <- obj
    for (gs in colnames(score_df)) {
      obj_d@meta.data[[gs]] <- score_df[[gs]]
    }
    
    p <- SpatialFeaturePlot(obj_d,
                            features = colnames(score_df),
                            crop = TRUE,
                            pt.size.factor = 2.5,
                            image.alpha = 0)
    ggsave(file.path(sample_dir, paste0(sample, "_DECOUPLE_", stat, ".pdf")), p, width = 12, height = 12) 
    saveRDS(obj_d, file = file.path(sample_dir, paste0(sample, "_DECOUPLE_", stat, ".RDS")))
  }
}
