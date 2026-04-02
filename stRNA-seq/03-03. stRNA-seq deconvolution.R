# Author : Depro Das, 3DBM Lab, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(dplyr) 
library(tidyverse) 
library(ggplot2) 
library(ggpubr) 
library(Seurat) 
library(DOTr) 

# ── Data input ──────────────────────────────────────────────────────────────── 

# Visium data preparation 

sample_dir <- "C:\\01. Research\\01. mIDH (hypoxia)\\02. Analysis\\03. Analysis (stRNA-seq)\\02. Spatial enrichment"
st_visum <- readRDS(file.path(sample_dir, "MGH259.RDS"))
st_count <- st_visum@assays$Spatial$counts 
st_coord <- GetTissueCoordinates(st_visum)

# Single-cell data preparation 

sc_data <- readRDS("mIDH total (CCA integrated Seurat 4.3.0).RDS")
sc_data <- subset(sc_data, subset = !ann_mps.mye_general %in% c("Unknown", "Glioma")) 
sc_count <- sc_data@assays[["RNA"]]@counts
sc_coord <- sc_data@meta.data[["ann_mps.mye_general"]]
class(sc_coord) 

# ── Run DOT ─────────────────────────────────────────────────────────────────── 

# Create and run DOT object 

dot.srt <- setup.srt(srt_data = st_count, srt_coords = st_coord) 
dot.ref <- setup.ref(ref_data = sc_count, ref_annotations = sc_coord, ref_subcluster_size = 10) 
dot_obj <- create.DOT(dot.srt, dot.ref)
saveRDS(dot_obj, "MGH259_DOTr.RDS") 

# Abundance weight: a larger value more closely matches the abundance of cell types in the spatial data to those in the reference data

# Maximum size of spots (20 is usually sufficiently large for Visium slides) 

dot_loww.res <- run.DOT.lowresolution(dot_obj, ratios_weight = 0, max_spot_size = 20, iterations = 10)
dot_high.res <- run.DOT.highresolution(dot_obj, ratios_weight = 0, iterations = 30)
dim(dot_high.res@weights)

# Plotting celltypes (basic) 

plot_data <- GetTissueCoordinates(st_visum) 
plot_data$celltype <- colnames(dot_high.res@weights)[apply(dot_high.res@weights, 1, which.max)]

p_basic <- ggplot(plot_data, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 1.5) +
  theme_bw() 
p_basic  

# ── Higher annotation ploting ───────────────────────────────────────────────── 

# Celltype map at a higher annotation level 

groups <- sapply(colnames(dot_high.res@weights), function(x) stringr::str_split(x, " ")[[1]][1])
agg_weights <- t(rowsum(t(dot_high.res@weights), groups))
plot_data$group <- colnames(agg_weights)[apply(agg_weights, 1, which.max)]

p_high_ann <- ggplot(plot_data, aes(x = x, y = y, color = group))+
  geom_point(size = 1.25) + 
  scale_color_manual(values = c("#99C945", "#CF3E53", "#F9A729", "#2CA02C", "#86BBD8", "#A26DC2", "#8C564B", "#CC61B0", "#5D69B1", "#66D2CE", "#FF999D", "#229F77")) +
  theme_bw() + 
  theme(panel.background = element_rect(fill = 'white'), 
        panel.grid = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank()) 
p_high_ann 
ggsave(filename = "1. MGH259_surface plot all cells.pdf", plot = p_high_ann, width = 5.5, height = 4)

# Add high resolution classification to spatial data  

st_visum <- AddMetaData(st_visum, metadata = plot_data$group, col.name = "ann_mps.mye_general")

# Data modification for plotting cell fractions  

fraction_idh <- dot_high.res@weights %>% as.data.frame() 
spatial_meta <- st_visum@meta.data
st_merged_df <- cbind(fraction_idh, spatial_meta) 
st_longgg_df <- st_merged_df %>%
  pivot_longer(cols = all_of(c('AC_like', 'OC_like', 'NPC_like', 'Cycling', 'Macrophage', 'Microglia', 'DCs', 'MDSC', 'Oligo', 'Endo', 'TCells', 'Pericytes')), 
               names_to = 'Celltypes', 
               values_to = 'Fraction')

# Plotting cell fractions 

paircomparison <- list(c("True_hypoxia", "Normoxia"), c("True_hypoxia", "Mild_hypoxia"), c("Mild_hypoxia", "Normoxia"))
st_longgg_df$Hypoxia_class <- factor(st_longgg_df$Hypoxia_class, levels = c("Normoxia", "Mild_hypoxia", "True_hypoxia"))

p_frac <- ggplot(st_longgg_df, aes(x = factor(Hypoxia_class), y = Fraction, fill = factor(Celltypes))) + 
  geom_bar(position = "fill", stat = "identity", width = 0.9) + 
  scale_fill_manual(values = c("#CF3E53", "#F9A729", "#2CA02C", "#86BBD8", "#A26DC2", "#8C564B", "#CC61B0", "#5D69B1", "#41644A", "#99C945", "#1B56FD", "#E53888", "#229F77", "black")) +
  theme_bw() + 
  stat_compare_means(comparisons = paircomparison) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
p_frac
ggsave(filename = "2. MGH259_fractions all cells.pdf", plot = p_frac, width = 4, height = 5)
