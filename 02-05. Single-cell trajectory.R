# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(Seurat)
library(Matrix)
library(sparseMatrixStats)
library(slingshot) 
library(ggplot2) 
library(ggplotify) 
library(viridis) 
library(patchwork) 
library(ComplexHeatmap) 
library(circlize) 

# ── Prepare data ──────────────────────────────────────────────────────────────

sc_inew <- readRDS("mIDH glioma (CCA integrated Seurat 4.3.0).rds")

# Dimension metadata (for plotting the labels later)

df_umap <- as.data.frame(sc_inew@reductions$umap@cell.embeddings)
colnames(df_umap) <- c("UMAP_1", "UMAP_2")
df_umap$cluster <- factor(sc_inew$QAD_class) 

# Cluster centroids (for plotting the labels later)

mm <- sparse.model.matrix(~ 0 + factor(sc_inew$QAD_class))
colnames(mm) <- levels(factor(sc_inew$QAD_class))
centroids_2d <- as.matrix(t(t(sc_inew@reductions$umap@cell.embeddings) %*% mm) / Matrix::colSums(mm))

centroids_df <- data.frame(UMAP_1 = centroids_2d[, 1],
                           UMAP_2 = centroids_2d[, 2],
                           cluster = rownames(centroids_2d),
                           stringsAsFactors = FALSE) 

# ── Run trajectory ──────────────────────────────────────────────────────────── 

# Run trajectory (dist.method = "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn")

set.seed(123) 
lineages <- as.SlingshotDataSet(getLineages(data = sc_inew@reductions$umap@cell.embeddings,
                                            clusterLabels = sc_inew$QAD_class,
                                            dist.method = "mnn", 
                                            end.clus = c("D_stage"), 
                                            start.clus = c("Q_stage"))) 

lineages@reducedDim <- sc_inew@reductions$umap@cell.embeddings

# Lineage dataframe 

lineage_dfs <- do.call(rbind, lapply(seq_along(lineages@lineages), function(i) {
  clus_path <- lineages@lineages[[i]]
  clus_path <- clus_path[clus_path %in% centroids_df$cluster]
  if (length(clus_path) < 2) return(NULL) 
  df_finals <- centroids_df[match(clus_path, centroids_df$cluster), ] 
  df_finals$lineage <- paste0("Lineage_", i) 
  df_finals$ordered <- seq_len(nrow(df_finals))
  df_finals
}))

# Plot 

set.seed(123) 
pal <- c(scales::hue_pal()(8), RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2")) 
pal <- rep(sample(pal, length(pal)), 200)

p_lineage <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cluster), size = 0.5, alpha = 0.8) +
  geom_path(data = lineage_dfs, aes(group = lineage), color = "black", linewidth = 1) +
  geom_point(data = lineage_dfs, color = "black", size = 4) +
  geom_text(data = centroids_df, aes(label = cluster), color = "white", fontface = "bold", size = 3) +
  scale_color_manual(values = pal) +
  theme_classic() 
p_lineage 
ggsave("01. Lineages QAD umap.pdf", plot = p_lineage, width = 4, height = 2.8, units = "in")

# ── Define principal curves ─────────────────────────────────────────────────── 

d.curves <- as.SlingshotDataSet(getCurves(data = lineages,
                                          thresh = 1e-1,
                                          stretch = 1e-1,
                                          allow.breaks = F,
                                          approx_points = 100))
d.curves 

# Curve dataframe

curve_df <- do.call(rbind, lapply(seq_along(d.curves@curves), function(i) { 
  crv <- d.curves@curves[[i]]
  data.frame(UMAP_1 = crv$s[, 1], 
             UMAP_2 = crv$s[, 2], 
             curve = paste0("Curve_", i))}))

p_curves <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cluster), size = 0.5, alpha = 0.7) +
  geom_path(data = curve_df, aes(group = curve), color = "black", linewidth = 1.2) +
  scale_color_manual(values = pal) +
  theme_classic()
p_curves
ggsave("02. Curves QAD umap.pdf", plot = p_curves, width = 4.23, height = 3, units = "in")

# ── Pseudotime ────────────────────────────────────────────────────────────────  

# Compute differentiation pseudotime

pseudot.run <- slingPseudotime(d.curves, na = TRUE)

# Average pseudotime (for multiple curves) 

df_umap$pseudotime <- rowMeans(pseudot.run, na.rm = TRUE)
df_umap$pseudotime <- df_umap$pseudotime / max(df_umap$pseudotime, na.rm = TRUE)

p_pseu.time <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = pseudotime), size = 0.5) +
  geom_path(data = curve_df, aes(group = curve), color = "black", linewidth = 1) +
  scale_color_viridis(option = "D", name = "Pseudotime", direction = 1) +
  theme_classic()
p_pseu.time
ggsave("03. Pseudotime QAD average umap.pdf", plot = p_pseu.time, width = 4.25, height = 3, units = "in")

# ── Separately assess each curve ────────────────────────────────────────────── 

colnames(slingPseudotime(d.curves)) 

nums_lineages <- ncol(slingPseudotime(d.curves))
lineage_plots <- list() 

for (i in 1:nums_lineages) {
  pseudotime_curv <- slingPseudotime(d.curves)[, i] 
  pseudotime_plot <- pseudotime_curv
  pseudotime_plot[is.na(pseudotime_plot)] <- NA  # non-lineage cells as null 
  df_umap$pt_plot <- pseudotime_plot
  df_lineg <- df_umap[!is.na(df_umap$pt_plot), ]
  curve_df <- data.frame(UMAP_1 = d.curves@curves[[i]]$s[,1],
                         UMAP_2 = d.curves@curves[[i]]$s[,2]) 
  
  p_ling.i <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = pt_plot), size = 0.5, na.rm = FALSE) +
    scale_color_viridis(option = "D", na.value = "grey80", name = "Pseudotime", direction = 1) +
    geom_density_2d(data = df_lineg, aes(UMAP_1, UMAP_2), color = "black", linewidth = 0.6) + 
    geom_path(data = curve_df, aes(UMAP_1, UMAP_2), linewidth = 1.3, color = "black") +
    theme_classic() +
    ggtitle(colnames(slingPseudotime(d.curves))[i])
  
  lineage_plots[[i]] <- p_ling.i
} 

names(lineage_plots) <- colnames(slingPseudotime(d.curves))

# Lineage plots together 

p.lin <- wrap_plots(lineage_plots, ncol = 1)
p.lin
ggsave("04. Pseudotime QAD separate umap.pdf", plot = p.lin, width = 8, height = 6, units = "in")

# Extract pseudotime

pseudo_df <- as.data.frame(pseudot.run)
pseudo_df <- cbind(pseudo_df, df_umap[, c("UMAP_1", "UMAP_2", "cluster")])
write.csv(pseudo_df, file = "Result pseudotime (QAD).csv")

# ── QAD stages mapping ──────────────────────────────────────────────────────── 

# QAD results 

metadata_cols <- c('Q_stage', 'A_stage', 'D_stage', 'metaprog_4_class', 'QAD_class', 'metaprog_QAD')
stage_feature <- FetchData(object = sc_inew, vars = metadata_cols) 
stage_feature <- cbind(pseudo_df, stage_feature) 

# Progeny results 

progn_results <- read.csv("Result Pathway activity (mlm).csv")
progn_results <- progn_results %>% 
  select(source, condition, score) %>% 
  pivot_wider(names_from = source, values_from = score) %>% 
  column_to_rownames(var = "condition") 

# Data manipulation (common) 

combin_feature <- cbind(stage_feature, progn_results) 
combin_feature$pseudotime <- combin_feature$Lineage1 

# Data manipulation (heatmap) 

qad.ht_feature <- combin_feature %>% 
  filter(!is.na(pseudotime)) %>% filter(metaprog_QAD %in% c('AC_Q', 'AC_A', 'AC_D')) %>% 
  select(pseudotime, metaprog_QAD, metaprog_4_class) %>% 
  filter(pseudotime <= 9) 

# Map pseudotime values to row indices 

pt_final <- as.numeric(qad.ht_feature$pseudotime) 
pt_matix <- as.matrix(pt_final)
colnames(pt_matix) <- "pseudotime"

pt_marks <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
pt_marks <- pt_marks[pt_marks >= min(pt_final, na.rm = TRUE) & pt_marks <= max(pt_final, na.rm = TRUE)]

pt_indix <- sapply(pt_marks, function(x) {
  which.min(abs(pt_final - x))
})

pt_indix <- unique(pt_indix)
pt_marks <- pt_marks[seq_along(pt_indix)]

psu_ht.row_ann <- rowAnnotation(pseudotime = anno_mark(at = pt_indix, 
                                                       labels = as.character(pt_marks), 
                                                       labels_gp = gpar(fontsize = 8), 
                                                       link_width = unit(2, "mm")))

qad_ht.row_ann <- rowAnnotation(metaprog_QAD = qad.ht_feature$metaprog_QAD, 
                                metaprog_4_class = qad.ht_feature$metaprog_4_class, 
                                col = list(metaprog_QAD = c(AC_Q = "#ffe000", AC_A = "#a3ba45", AC_D = "#8dd8ff"), metaprog_4_class = c(AC_like = "#e41a1c", OC_like = "#984ea3", NPC_like = "#e41a1c", Cycling = "#377eb8")), 
                                annotation_width = unit(c(3, 3), "mm"))

set.seed(123) 
p_ht.qad <- Heatmap(pt_matix, 
                    col = colorRamp2(seq(min(pt_final), max(pt_final), length.out = 5), viridis(5)),
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    show_row_names = FALSE,
                    left_annotation  = psu_ht.row_ann,
                    right_annotation = qad_ht.row_ann, 
                    border = TRUE) 
p_ht.qad
p_ht.qad <- as.ggplot(p_ht.qad)
ggsave("05. QAD mapping to pseudotime.pdf", plot = p_ht.qad, width = 3, height = 6, units = "in")

# ── Pathway activity mapping ────────────────────────────────────────────────── 

stage_columns <- c('Q_stage', 'A_stage', 'D_stage')
progn_columns <- c('Androgen','EGFR','Estrogen','Hypoxia','JAK-STAT', 'MAPK','NFkB','PI3K','TGFb','TNFa','Trail','VEGF','WNT','p53')

activ_feature <- combin_feature %>% filter(!is.na(pseudotime)) %>% filter(metaprog_QAD %in% c('AC_Q', 'AC_A'))

activ_feature <- activ_feature %>%
  mutate(across(all_of(c(stage_columns, progn_columns)), scale)) %>%
  mutate(cell_id = row_number()) 

# Stages with one pathway per column

activ_mappped <- lapply(progn_columns, function(pathway) {
  activ_feature %>% 
    select(pseudotime, all_of(stage_columns), all_of(pathway)) %>%
    pivot_longer(cols = c(all_of(stage_columns), all_of(pathway)), names_to = "feature", values_to = "value") %>%
    mutate(facet_group = pathway, feature_type = ifelse(feature %in% stage_columns, "stage", "pathway"))
}) %>% 
  bind_rows()

# Plot mapped activity 

p_act.mappped <- ggplot(activ_mappped, aes(x = pseudotime, y = value, color = feature)) +
  geom_smooth(data = subset(activ_mappped, feature_type == "stage"), method = "loess", se = T, size = 1) +
  geom_smooth(data = subset(activ_mappped, feature_type == "pathway"), method = "loess", se = T, size = 1, color = "#3366ff") +
  facet_wrap(~facet_group, scales = "free_y") +
  theme_bw() + 
  xlim(0, 9) +
  scale_color_manual(values = c("Q_stage" = "#ff00f3", "A_stage" = "#15cc26", "D_stage" = "orange"))
p_act.mappped
ggsave("06. Activity mapping.pdf", plot = p_act.mappped, width = 12, height = 6, units = "in")

