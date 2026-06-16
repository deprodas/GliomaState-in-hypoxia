# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(liana) 
library(tidyverse) 
library(SingleCellExperiment)
library(reticulate)
library(magrittr) 
library(openxlsx)
library(ExperimentHub) 
library(Seurat) 
library(ggplot2) 
library(ggplotify) 
library(patchwork) 
library(ComplexHeatmap) 
library(circlize)
library(viridis)
library(nichenetr) 
library(ggrepel) 
library(cowplot)

# ── Prepare data ──────────────────────────────────────────────────────────────   

sc_data <- readRDS("mIDH total (CCA integrated Seurat 4.3.0).RDS")
sc_data <- subset(sc_data, subset = !ann_mps.mye_general %in% c("Unknown", "Glioma"))

# ── Liana basic (consensus) ─────────────────────────────────────────────────── 

show_methods()
show_resources() 
DefaultAssay(sc_data) <- "RNA"
Idents(sc_data) <- "ann_mps.mye_general"  

sce_obj <- as.SingleCellExperiment(sc_data) 

# Run liana on consensus (default) resource

liana_test <- liana_wrap(sce_obj, resource = c("Consensus"), idents_col = "ann_mps.mye_general", min_cells = 2) 
liana_test %>% dplyr::glimpse() 

# Aggregate these results and plot 

liana_agg <- liana_test %>% liana_aggregate()
liana_fil <- liana_agg %>% dplyr::filter(receptor.complex == "EGFR")
liana_fil %>% dplyr::count(target) 

dot.p_conc.def <- liana_fil %>% 
  liana_dotplot(source_groups = c('Microglia', 'Macrophage', 'DCs', 'MDSC'),
                target_groups = c('AC_like', 'NPC_like', 'OC_like', 'Cycling'), 
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dot.p_conc.def 
ggsave(filename = "1. LR EGFR specific (DefaultDB Liana).pdf", plot = dot.p_conc.def, width = 16, height = 10, units = c("in")) 

# Plot manually 

dot.df_conc.def <- dot.p_conc.def[["data"]] 

dot.p_conc.manual <- ggplot(dot.df_conc.def, aes(x = target, y = interaction, size = specificity, color = magnitude)) +
  geom_point(alpha = 1) +  
  theme_minimal() + 
  scale_colour_gradient2(low = "Darkblue", mid = "white", high = "red", midpoint = 0.75) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~source, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dot.p_conc.manual
ggsave(filename = "1. LR EGFR specific (DefaultDB ggplot).pdf", plot = dot.p_conc.manual, width = 10, height = 6, units = c("in")) 

# Frequency heatmap 

liana_trunc <- liana_agg %>% dplyr::filter(aggregate_rank <= 0.01) 
heat_freq(liana_trunc) 

# ── Liana based on cellchat database ────────────────────────────────────────── 

# Run cellchat (alone) 

ccdb_test <- liana_wrap(sce_obj,
                        # method = 'call_cellchat',
                        resource = c('CellChatDB'), 
                        idents_col = "ann_mps.mye_general", 
                        permutation.params = list(nperms = 100, parallelize = FALSE, workers = 4),
                        expr_prop = 0.05)

# Interactions of interest

ccdb_agg <- ccdb_test %>% liana_aggregate()
ccdb_fil <- ccdb_agg %>% dplyr::filter(receptor.complex == "EGFR")
ccdb_fil %>% dplyr::count(target) 
write.csv(ccdb_agg, file = "Result Liana ann_mps.mye_general (CellChatDB).csv")

# Frequency heatmap 

ccdb_trunc <- ccdb_agg %>% dplyr::filter(aggregate_rank <= 0.01) 
heat_freq(ccdb_trunc) 

# Interactions by source and target cells

dot.p_cc.def <- ccdb_fil %>% 
  liana_dotplot(source_groups = c('Microglia', 'Macrophage', 'DCs', 'MDSC'),
                target_groups = c('AC_like', 'NPC_like', 'OC_like', 'Cycling'), 
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dot.p_cc.def 
ggsave(filename = "2. LR EGFR specific (CellchatDB Liana).pdf", plot = dot.p_cc.def, width = 14, height = 6, units = c("in"))

# Plot manually 

dot.df_cc.def <- dot.p_cc.def[["data"]] 
dot.df_cc.def$new_pvalue <- ifelse(dot.df_cc.def$cellphonedb.pvalue == 0, 1e-5, dot.df_cc.def$cellphonedb.pvalue)

dot.p_cc.manual <- ggplot(dot.df_cc.def, aes(x = target, y = interaction, size = -log10(new_pvalue), color = specificity)) +
  geom_point(alpha = 1) +  
  theme_minimal() + 
  scale_color_viridis(option = "D") + 
  facet_wrap(~source, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dot.p_cc.manual 

bar.p_cc.manual <- ggplot(dot.df_cc.def, aes(x = target, y = specificity, fill = target)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~source, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") 
bar.p_cc.manual

p_cc.comb.1 <- bar.p_cc.manual / dot.p_cc.manual + plot_layout(heights = c(1.5, 3))
p_cc.comb.1
ggsave(filename = "2. LR EGFR specific (CellchatDB ggplot).pdf", plot = p_cc.comb.1, width = 10, height = 3.5, units = c("in"))

# Add cellchat communication probability 

cc.obj_dir <- "~/Desktop/1. mIDH (hypoxia)/2. Analysis/02. Analysis (Single-cell)/07. Cell-cell (Cellchat)" 
cellchat <- readRDS(file.path(cc.obj_dir, "Cellchat ann_mps.mye_general pop-T (CCA integrated Seurat 4.3.0).RDS"))

str(cellchat@net$pval)

prob_array <- cellchat@net$prob
pval_array <- cellchat@net$pval
class(prob_array) 
sources <- dimnames(prob_array)[[1]]
targets <- dimnames(prob_array)[[2]]
intract <- dimnames(prob_array)[[3]]

df_indices <- as.data.frame(which(!is.na(prob_array), arr.ind = TRUE))
index_mtrx <- as.matrix(df_indices[, c("dim1", "dim2", "dim3")])

df_indices$comprob <- prob_array[index_mtrx]
df_indices$pvalues <- pval_array[index_mtrx]
df_indices$sources <- sources[df_indices$dim1]
df_indices$targets <- targets[df_indices$dim2]
df_indices$intract <- intract[df_indices$dim3]

# Filter liana interactions on cellchatdb 

df_specific <- df_indices %>% 
  dplyr::filter(intract %in% c("HBEGF_EGFR", "AREG_EGFR", "TGFA_EGFR") & 
                  sources %in% c('Microglia', 'Macrophage', 'DCs', 'MDSC') & 
                  targets %in% c('AC_like', 'NPC_like', 'OC_like', 'Cycling')) 

df_specific$new_pvalue <- ifelse(df_specific$pval == 0, 1e-5, df_specific$pval)

bar.p_cc.prob <- ggplot(df_specific, aes(x = targets, y = comprob, fill = targets)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~sources, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") 
bar.p_cc.prob

p_cc.comb.2 <- bar.p_cc.prob / dot.p_cc.manual + plot_layout(heights = c(1.5, 3))
p_cc.comb.2
ggsave(filename = "3. LR EGFR specific-2 (CellchatDB ggplot).pdf", plot = p_cc.comb.2, width = 10, height = 4, units = c("in"))

# ── LR expression in single-cell ────────────────────────────────────────────── 

lr_genes <- c("EGFR", "TGFA", "HBEGF", "AREG") 
DefaultAssay(sc_data) <- 'RNA'

# Basic expression 

v0_lr <- VlnPlot(sc_data, 
                 group.by = "ann_mps.mye_general", 
                 fill.by = "ident",
                 sort = "decreasing", 
                 features = lr_genes, 
                 pt.size = 0) 
v0_lr 
ggsave(filename = "4. LR genes violin (EGFR specific).pdf", plot = v0_lr, width = 12, height = 6, units = c("in")) 

r0_lr <- RidgePlot(sc_data, 
                   group.by = "ann_mps.mye_general", 
                   features = lr_genes,
                   sort = "decreasing") 
r0_lr
ggsave(filename = "4. LR genes ridge (EGFR specific).pdf", plot = r0_lr, width = 12, height = 8, units = c("in")) 

# Stack expression 

v1_lr <- VlnPlot(sc_data, 
                 group.by = "ann_mps.mye_general", 
                 fill.by = "ident",
                 features = lr_genes, 
                 sort = "decreasing", 
                 pt.size = 0, 
                 stack = TRUE, 
                 flip = TRUE) 
v1_lr 
ggsave(filename = "4. LR genes stack (EGFR specific).pdf", plot = v1_lr, width = 6, height = 4, units = c("in")) 

# Non-zero expression 

vlnplotCustom <- function(seurat_obj, features, group.by = "ident", ncol = 2) {
  plots <- lapply(features, function(gene) { 
    expr_data <- FetchData(seurat_obj, vars = c(gene, group.by))
    colnames(expr_data)[2] <- "group"
    expr_data <- expr_data[expr_data[[gene]] > 0, ] 
    
    group_order <- expr_data %>%
      group_by(group) %>%
      summarize(med = median(.data[[gene]])) %>%
      arrange(med) %>%
      pull(group)
    
    expr_data$group <- factor(expr_data$group, levels = group_order)
    
    p <- ggplot(expr_data, aes(x = group, y = .data[[gene]], fill = group)) +
      geom_violin(scale = "width", trim = TRUE) +
      theme_bw() +
      ylab(paste(gene, "Expression (non-zero)")) +
      xlab("") +
      ggtitle(gene) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") 
    
    return(p)
  }) 
  wrap_plots(plots, ncol = ncol)
} 

v2_lr <- vlnplotCustom(sc_data, 
                       features = c("EGFR", "HBEGF", "AREG"), 
                       group.by = "ann_mps.mye_general", 
                       ncol = 5) 
v2_lr
ggsave(filename = "5. LR genes violin-nozero ggplot (EGFR specific).pdf", plot = v2_lr, width = 15, height = 10, units = c("in")) 

# ── Liana between two cell categories ───────────────────────────────────────── 

# Subset seurat object to keep the interaction between desired cell types 

sc.sub_data <- subset(sc_data, subset = ann_mps.mye_general == c("Microglia", "AC_like"))
sce.sub_obj <- as.SingleCellExperiment(sc.sub_data, assay = "RNA")

# Run liana on subset seurat 

liana_nich_res <- liana_wrap(sce.sub_obj, 
                             permutation.params = list(nperms = 100, parallelize = FALSE, workers = 4),
                             expr_prop = 0.05, 
                             idents_col = "ann_mps.mye_general", 
                             resource = c('CellChatDB')) %>% liana_aggregate()

# Filter and plot results 

re_tumor_reslt <- liana_nich_res %>%
  subset(source == "Microglia" & target == "AC_like") %>%
  dplyr::rename(ligand = ligand.complex, receptor = receptor.complex)

n <- 20
top_n_re_tumor <- re_tumor_reslt %>%
  arrange(aggregate_rank) %>%
  slice_head(n = n) %>%
  mutate(id = fct_inorder(paste0(ligand, " -> ", receptor)))

p.top.int <- top_n_re_tumor %>%
  ggplot(aes(y = aggregate_rank, x = id)) +
  geom_bar(stat = "identity") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 1))
p.top.int
ggsave(filename = "6. Top interactions (microglia-re).pdf", plot = p.top.int, width = 8, height = 4, units = c("in")) 

# ── Validate liana by nichenet ligand potential ─────────────────────────────── 

# Selected top ligands and filter to those included in nicheNet's ligand-target matrix

ligand_target_matrix <- readRDS("ligand_target_matrix.rds")

ligands <- unique(top_n_re_tumor$ligand)
ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]
ligands 

# Prepare nichenet inputs

express_mat <- sc_sub@assays[["RNA"]]@counts %>% as.data.frame()
sample_info <- sc_sub@meta.data 
sample_info$barcodes <- rownames(sample_info) 

# Before running nicheNet, we also need to define a list of background genes 

background_genes <- express_mat[, rownames(sample_info)[sample_info$ann_mps.mye_general == "AC_like"], drop = FALSE] %>%
  apply(1, function(x) log2(mean(10 * (2^x - 1)) + 1)) %>%
  .[. >= 4] %>%
  names() 

# Gene set of interest (top receptors from cellchat for which the expression is affected due to communication with other cells) 

ccdb_tar.lig <- ccdb_agg %>% dplyr::filter(source %in% "Microglia" & target %in% "AC_like") 
geneset_oi <- ccdb_tar.lig %>% pull(receptor.complex)

# Run nichenet 

nichenet_act <- predict_ligand_activities(geneset = geneset_oi,
                                          background_expressed_genes = background_genes,
                                          ligand_target_matrix = ligand_target_matrix, 
                                          potential_ligands = ligands)

# Plot default heatmap 

vis_liana_nichenet <- top_n_re_tumor %>%
  inner_join(nichenet_act, by = c("ligand" = "test_ligand")) %>%
  arrange(pearson) %>%
  mutate(ligand = fct_inorder(ligand))

nichenet_scores_plot <- vis_liana_nichenet %>%
  group_by(ligand) %>%
  summarize(pearson = mean(pearson)) %>%
  ggplot(aes(y = ligand, x = pearson)) +
  geom_bar(stat = "identity") +
  ggtitle("NicheNet") +
  xlab("Pearson's score") +
  theme_cowplot() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))

liana_receptor_heatmap <- vis_liana_nichenet %>%
  ggplot(aes(y = ligand, x = receptor, fill = aggregate_rank)) +
  geom_tile() +
  theme_cowplot() +
  ggtitle("LIANA") +
  ylab("Ligand") + xlab("Receptor") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour = "gray", linetype = 2),
        legend.position = "left")

ht.custom_niche <- plot_grid(liana_receptor_heatmap, nichenet_scores_plot, align = "h", nrow = 1, rel_widths = c(0.8,0.3))
ht.custom_niche
ggsave(filename = "7. Top L-R custom heatmap (microglia-re).pdf", plot = ht.custom_niche, width = 9.5, height = 6, units = c("in"))


# Plot complex heatmap 

heatmap_mat <- vis_liana_nichenet %>%
  select(ligand, receptor, aggregate_rank) %>%
  pivot_wider(names_from = receptor, values_from = aggregate_rank) %>%
  column_to_rownames("ligand") %>%
  as.matrix() 

heatmap_mat[is.na(heatmap_mat)] <- 0

pearson_vals <- vis_liana_nichenet %>%
  group_by(ligand) %>%
  summarize(pearson = mean(pearson)) %>%
  deframe()

pearson_vals <- pearson_vals[rownames(heatmap_mat)] 

row_ha <- rowAnnotation(pearson = anno_barplot(pearson_vals, border = FALSE, gp = gpar(fill = "steelblue")), 
                        annotation_label = "Pearson")

ht.complex_niche <- Heatmap(heatmap_mat,
                            name = "Aggregate rank", 
                            right_annotation = row_ha,
                            cluster_rows = TRUE,
                            cluster_columns = TRUE,
                            col = colorRamp2(c(0.0, 0.25, 0.5, 0.75, 1.0), c("white", "#caf0f8", "#00b4d8", "#0077b6", "#03045e")),
                            row_names_side = "left", 
                            border = T)
ht.complex_niche <- as.ggplot(ht.complex_niche) 
ht.complex_niche 
ggsave(filename = "7. Top L-R complex heatmap (microglia-re).pdf", plot = ht.complex_niche, width = 8, height = 6, units = c("in"))

# Plot barplots

vis_liana_long <- vis_liana_nichenet %>%
  select(id, pearson, aggregate_rank) %>%
  pivot_longer(cols = c(pearson, aggregate_rank), names_to = "metric", values_to = "value")

p.agr_bar <- ggplot(vis_liana_long, aes(x = value, y = reorder(id, value), fill = metric)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_x") +
  theme_minimal(base_size = 13)
p.agr_bar 
ggsave(filename = "8. L-R barplot (microglia-re).pdf", plot = p.agr_bar, width = 7, height = 5.5, units = c("in"))

