# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(R.utils) 
library(openxlsx) 
library(dplyr) 
library(tidyverse) 
library(ggplot2) 
library(ggpubr) 
library(Seurat) 

# ── Spatial data input ──────────────────────────────────────────────────────── 

# Create spatial seurat objects (Greenwald et al.)

createSpatial_tirosh <- function(sample_name, base_dir) {
  data_path <- file.path(base_dir, sample_name)
  seurat_obj <- Load10X_Spatial(
    data.dir = data_path,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    bin.size = NULL,
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL
  )
  seurat_obj@project.name <- sample_name
  return(seurat_obj)
}

BWH23 <- createSpatial_tirosh("BWH23", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Tirosh samples\\BWH23")
BWH24 <- createSpatial_tirosh("BWH24", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Tirosh samples\\BWH24")
BWH25 <- createSpatial_tirosh("BWH25", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Tirosh samples\\BWH25")
BWH28 <- createSpatial_tirosh("BWH28", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Tirosh samples\\BWH28")
BWH35 <- createSpatial_tirosh("BWH35", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Tirosh samples\\BWH35")
MGH259 <- createSpatial_tirosh("MGH259", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Tirosh samples\\MGH259") 


# Create spatial seurat objects (He et al.) 

wang_lg1 <- "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Wang samples\\LG1\\filtered_feature_bc_matrix"
list.files(wang_lg1) 

createSpatial_wang <- function(sample_name, base_dir, region = "LGG") {
  data_dir <- file.path(base_dir, sample_name, "filtered_feature_bc_matrix")
  expression_matrix <- Read10X(data.dir = data_dir)
  
  seurat_obj <- CreateSeuratObject(counts = expression_matrix, assay = "Spatial", project = sample_name)
  seurat_obj$slice <- 1
  seurat_obj$region <- region
  
  img_path <- file.path(base_dir, sample_name)
  img <- Seurat::Read10X_Image(image.dir = img_path)
  Seurat::DefaultAssay(object = img) <- "Spatial"
  img <- img[colnames(seurat_obj)] 
  seurat_obj[["image"]] <- img
  
  return(seurat_obj)
}

LG1 <- createSpatial_wang("LG1", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Wang samples\\LG1")
LG2 <- createSpatial_wang("LG2", "C:\\01. Hypoxia analysis\\3. Spatial Transcriptomics\\Wang samples\\LG2")


# ── Separately preprocessing ──────────────────────────────────────────────────

PreprocessSpatialSample <- function(seurat_obj, output_file) {
  seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
  saveRDS(seurat_obj, file = output_file)
  return(seurat_obj)
}

BWH23_dub <- PreprocessSpatialSample(BWH23, "BWH23_basic.RDS")
BWH24_dub <- PreprocessSpatialSample(BWH24, "BWH24_basic.RDS")
BWH25_dub <- PreprocessSpatialSample(BWH25, "BWH25_basic.RDS")
BWH28_dub <- PreprocessSpatialSample(BWH28, "BWH28_basic.RDS")
BWH35_dub <- PreprocessSpatialSample(BWH35, "BWH35_basic.RDS")
MGH259_dub <- PreprocessSpatialSample(MGH259, "MGH259_basic.RDS")
LG1_dub <- PreprocessSpatialSample(LG1, "LG1_basic.RDS")
LG2_dub <- PreprocessSpatialSample(LG2, "LG2_basic.RDS")


# ── Merged preprocessing ────────────────────────────────────────────────────── 

# Merge after SCT normalization

merged_seurat <- merge(BWH23_dub, y = list(BWH24_dub, BWH25_dub, BWH28_dub, BWH35_dub, MGH259_dub, LG1_dub, LG2_dub), add.cell.ids = c("BWH23_dub", "BWH24_dub", "BWH25_dub", "BWH28_dub", "BWH35_dub", "MGH259_dub", "LG1_dub", "LG2_dub"), project = "LGG_merge")

# Save the individual clusters  

clust_ind <- SpatialDimPlot(merged_seurat, 
                            pt.size.factor = 3.5, 
                            image.alpha = 0, 
                            ncol = 4, 
                            label = F, 
                            label.size = 0) 
clust_ind
ggsave(filename = "1. Spatial Clusters (Merged Seurat).pdf", plot = clust_ind, width = 10, height = 20, units = c("in")) 


# Set variable features manually

VariableFeatures(merged_seurat) <- SelectIntegrationFeatures(object.list = list(BWH24, BWH25, BWH28, BWH35, MGH259, LG1, LG2), nfeatures = 3000)

merged_seurat <- RunPCA(merged_seurat, assay = "SCT", verbose = FALSE)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "pca", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, verbose = FALSE)
merged_seurat <- RunUMAP(merged_seurat, reduction = "pca", dims = 1:30)


# Create a normalized data layer

merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)

# Save the merged object 

saveRDS(merged_seurat, file = "LGG_merged_10xSeurat.RDS")

# Join layers together

merged_seurat <- JoinLayers(merged_seurat, assay = "Spatial")

# UMI count

spatial_counts <- merged_seurat@assays$Spatial$counts 
spatial_counts <- as.data.frame(spatial_counts)

# ── Spatial region classification ─────────────────────────────────────────────

# Run enrichment analysis in the merged count matrices (ulm)

mat <- spatial_counts %>% as.matrix() 
mat <- log2(mat+1)

genesets <- read.csv("LGG_new_genes.csv")
colnames(genesets) 
genesets <- genesets %>% select("HM_Hypoxia")
genesets <- genesets %>% pivot_longer(c(HM_Hypoxia), names_to = "source", values_to = "target")
genesets <- genesets %>% mutate(mor = 1)
genesets <- genesets %>% filter(target!=" ")
genesets <- genesets %>% filter(target!="") 

# Run all decoupler algorithm  

res_decouple <- decoupleR::decouple(mat = mat, 
                                    network = genesets, 
                                    .source = 'source', 
                                    .target = 'target', 
                                    minsize = 0)

res_decouple %>% dplyr::count(statistic) 
write.xlsx(res_decouple, file = "Result DecoupleR-all (Merged Seurat).xlsx")


# Convert the data into long fromat 

res_decouple <- openxlsx::read.xlsx("Result DecoupleR-all (Merged Seurat).xlsx")

res_long <- res_decouple %>% 
  dplyr::filter(statistic == "ulm") %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% 
  t() %>% 
  as.data.frame() 

res_long <- res_long %>% dplyr::select(HM_Hypoxia) %>% rename(hypoxia_score = HM_Hypoxia)


# Classification of group 

# Calculate median and standard deviation 

med_val <- median(res_long$hypoxia_score, na.rm = TRUE)
std_dev <- sd(res_long$hypoxia_score, na.rm = TRUE)

# Classify samples based on standard deviation thresholds

res_long <- res_long %>% 
  mutate(Hypoxia_class = case_when(hypoxia_score < med_val ~ "Normoxia", 
                                   hypoxia_score >= med_val & hypoxia_score <= (med_val + std_dev) ~ "Mild_hypoxia", 
                                   hypoxia_score > (med_val + std_dev) ~ "True_hypoxia")) 

res_long %>% dplyr::count(Hypoxia_class) 


# Visualization in a density plot  

dens_p1 <- ggplot(res_long, aes(x = hypoxia_score)) +
  geom_density(fill = "orange", alpha = 0.5) + 
  geom_vline(xintercept = med_val, linetype = "solid", color = "black", size = 1.2, alpha = 0.8) +
  geom_vline(xintercept = med_val + std_dev, linetype = "dashed", color = "red", size = 1.2) + 
  theme_minimal()
dens_p1
ggsave(filename = "2. Density Plot Cutoff (Merged Seurat).pdf", plot = dens_p1, width = 6, height = 3, units = c("in"))


# ── Spatial sample sorting ──────────────────────────────────────────────────── 

# Add patient or sample metadata 

rownames(res_long)
res_long <- res_long %>% rownames_to_column(var = "barcodes") 
res_long <- res_long %>% mutate(sample_ID = sub("_.*", "", barcodes))
write.csv(res_long, file = "Result Main Metadata (Merged Seurat).csv", row.names = F)


# Add metadata to seurat object 

merged_seurat$Hypoxia_class <- res_long$Hypoxia_class 
merged_seurat$sample_ID <- res_long$sample_ID 

# Save seurat object 

saveRDS(merged_seurat, file = "LGG_merged_10xSeurat.RDS")


# Proportion of hypoxia classes 

# Read the region distributation metadata 

res_long <- read.xlsx("Result Main Metadata (Merged Seurat).xlsx") 

# Fix the desired order

res_long$Hypoxia_class <- factor(res_long$Hypoxia_class, levels = c("Normoxia", "Mild_hypoxia", "True_hypoxia"))

df_pctg <- res_long %>%
  group_by(sample_ID, Hypoxia_class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_ID) %>%
  mutate(percentage = (count / sum(count)) * 100)

df_pctg <- df_pctg %>%
  group_by(sample_ID) %>%
  mutate(True_hypoxia_pctg = percentage[Hypoxia_class == "True_hypoxia"]) %>%
  ungroup() %>%
  arrange(True_hypoxia_pctg) 


# Stacked bar plot

hp_stack <- ggplot(df_pctg, aes(x = reorder(sample_ID, True_hypoxia_pctg), y = percentage, fill = Hypoxia_class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7")) +
  labs(title = "Proportion of hypoxia classes per sample", x = "Sample ID", y = "Percentage", fill = "Hypoxia class") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
hp_stack 
ggsave(filename = "4. Percentage Plot (Merged Seurat).pdf", plot = hp_stack, width = 4, height = 3, units = c("in")) 

# Surface plots of hypoxic zones 

hp_surf <- SpatialDimPlot(merged_seurat, 
                          group.by = "Hypoxia_class", 
                          pt.size.factor = 3.5, 
                          cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                          image.alpha = 0, 
                          ncol = 4, 
                          label = F, 
                          label.size = 0) 
hp_surf
ggsave(filename = "5. Surface Plot (Merged Seurat-Frontiers).pdf", plot = hp_surf, width = 18, height = 18, units = c("in")) 


DimPlot(merged_seurat, group.by = "Hypoxia_class", reduction = "tnse", label = TRUE) 


# Split seurat object by sample 

split_list <- SplitObject(merged_seurat, split.by = "sample_ID") 

BWH23 <- split_list[["BWH23"]]
BWH24 <- split_list[["BWH24"]]
BWH25 <- split_list[["BWH25"]]
BWH28 <- split_list[["BWH28"]]
BWH35 <- split_list[["BWH35"]]
MGH259 <- split_list[["MGH259"]] 
LG1 <- split_list[["LG1"]]
LG2 <- split_list[["LG2"]]

# Save individual seurat object 

for (sample_name in names(split_list)) {
  saveRDS(split_list[[sample_name]], file = paste0(sample_name, ".RDS"))
}


# ── Better looking plots ──────────────────────────────────────────────────────  

sp1 <- SpatialDimPlot(BWH23, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 5, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp1 
ggsave(filename = "5. Surface Plot (BWH23).pdf", plot = sp1, width = 8, height = 6, units = c("in")) 

sp2 <- SpatialDimPlot(BWH24, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 3, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp2 
ggsave(filename = "5. Surface Plot (BWH24).pdf", plot = sp2, width = 8, height = 6, units = c("in")) 

sp3 <- SpatialDimPlot(BWH25, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 3.75, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp3 
ggsave(filename = "5. Surface Plot (BWH25).pdf", plot = sp3, width = 8, height = 6, units = c("in")) 

sp4 <- SpatialDimPlot(BWH35, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 4.5, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp4 
ggsave(filename = "5. Surface Plot (BWH35).pdf", plot = sp4, width = 8, height = 6, units = c("in")) 

sp5 <- SpatialDimPlot(BWH28, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 3.5, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp5 
ggsave(filename = "5. Surface Plot (BWH28).pdf", plot = sp5, width = 8, height = 6, units = c("in")) 

sp6 <- SpatialDimPlot(MGH259, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 3, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp6 
ggsave(filename = "5. Surface Plot (MGH259).pdf", plot = sp6, width = 8, height = 6, units = c("in")) 

sp7 <- SpatialDimPlot(LG1, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 3.5, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp7 
ggsave(filename = "5. Surface Plot (LG1).pdf", plot = sp7, width = 8, height = 6, units = c("in")) 

sp8 <- SpatialDimPlot(LG2, 
                      group.by = "Hypoxia_class", 
                      pt.size.factor = 3, 
                      cols = c(Normoxia = "#6EC207", Mild_hypoxia = "#FFEB00", True_hypoxia = "#ED3EF7"),
                      image.alpha = 0, 
                      ncol = 4, 
                      label = F, 
                      label.size = 0) 
sp8 
ggsave(filename = "5. Surface Plot (LG2).pdf", plot = sp8, width = 8, height = 6, units = c("in")) 
