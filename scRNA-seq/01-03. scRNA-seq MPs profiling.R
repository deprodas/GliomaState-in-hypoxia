# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(Seurat) 
library(patchwork) 
library(ComplexHeatmap) 
library(circlize) 
library(ggplotify)
library(corrplot) 
library(dittoSeq) 
library(decoupleR) 
library(presto) 
library(msigdbr) 
library(fgsea) 
library(viridis)
library(tibble) 

# ── Data input ────────────────────────────────────────────────────────────────  

sc_inew <- readRDS("mIDH glioma (CCA integrated Seurat 4.3.0).rds") 
gs.prog <- read.xlsx("mIDH programs_depro.xlsx") 
gs.stdy <- read.xlsx("mIDH programs_itay.xlsx") 
gs.hypx <- read.xlsx("Significant hypoxia geneset.xlsx") 

DimPlot(sc_inew, reduction = "umap", label = T) 


# ── Defining transcriptional states ─────────────────────────────────────────── 

# Jaccard index 

gs.prog_jc <- lapply(gs.prog, function(x) unique(na.omit(x)))
gs.stdy_jc <- lapply(gs.stdy, function(x) unique(na.omit(x)))
gs.combine <- c(gs.prog_jc, gs.stdy_jc)
set_gnames <- names(gs.combine)

j.pairwise <- expand.grid(set1 = set_gnames, set2 = set_gnames, stringsAsFactors = FALSE) # expand to same pairs

j.pairwise$jaccard <- mapply(function(s1, s2) {
  genes_s1 <- gs.combine[[s1]]
  genes_s2 <- gs.combine[[s2]]
  length(intersect(genes_s1, genes_s2)) / length(union(genes_s1, genes_s2))
}, j.pairwise$set1, j.pairwise$set2)

jac_pair_df <- j.pairwise

jac_matx_df <- jac_pair_df %>%
  pivot_wider(names_from = set2, values_from = jaccard, values_fill = 0) %>%
  as.data.frame() %>% 
  column_to_rownames(var = "set1") %>% 
  as.matrix() 

# Plot jaccard heatmap 

set.seed(123) 

ht_jac <- Heatmap(jac_matx_df,
                  col = colorRamp2(c(0, 0.5, 1), c("white", "orange", "red")),
                  cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", jac_matx_df[i, j]), x, y, gp = gpar(fontsize = 10)) },
                  cluster_rows = T, 
                  cluster_columns = T, 
                  row_names_side = "right",
                  column_names_side = "top",
                  heatmap_legend_param = list(title = "Jaccard scores"), 
                  border = T)
ht_jac
ht_jac <- as.ggplot(ht_jac) 
ggsave(filename = "01. Jaccard heatmap.pdf", plot = ht_jac, width = 7.1, height = 6, units = c("in"))


# ── Module scoring and classification ───────────────────────────────────────── 

# Prepare genesets and run module scoring 

gs.comp <- list(AC_like = as.character(gs.prog$MetaProg1), 
                OC_like = as.character(gs.prog$MetaProg2), 
                Cycling = as.character(gs.prog$MetaProg4), 
                NPC_like = as.character(gs.prog$MetaProg3), 
                Hypoxia = as.character(gs.hypx$Hypoxia_all))  

gs.list <- lapply(gs.comp, function(x) list(as.character(x))) 

DefaultAssay(sc_inew) <- "SCT"

for (nm in names(gs.list)) {
  sc_inew <- AddModuleScore(sc_inew, features = gs.list[[nm]], name = nm) 
  oldname <- grep(paste0("^", nm), colnames(sc_inew@meta.data), value = TRUE)
  colnames(sc_inew@meta.data)[colnames(sc_inew@meta.data) == oldname] <- nm
}

# Umap enrichment of programs 

HP_ft <- FeaturePlot(sc_inew, 
                     reduction = "umap", 
                     min.cutoff = 0.03,
                     max.cutoff = 0.12,
                     features = c('Hypoxia')) 
HP_ft 

AS_ft <- FeaturePlot(sc_inew, 
                     reduction = "umap", 
                     min.cutoff = 0.0,
                     max.cutoff = 0.8,
                     features = c('AC_like')) 
AS_ft 

OL_ft <- FeaturePlot(sc_inew, 
                     reduction = "umap", 
                     min.cutoff = 0.0,
                     max.cutoff = 1.0,
                     features = c('OC_like')) 
OL_ft 

NPC_ft <- FeaturePlot(sc_inew, 
                      reduction = "umap", 
                      min.cutoff = 0.1,
                      max.cutoff = 0.6,
                      features = c('NPC_like')) 
NPC_ft 

Cyc_ft <- FeaturePlot(sc_inew, 
                      reduction = "umap", 
                      min.cutoff = 0.0,
                      max.cutoff = 1.0,
                      features = c('Cycling')) 
Cyc_ft

prg_ft <- wrap_plots(HP_ft, AS_ft, OL_ft, NPC_ft, Cyc_ft, ncol = 3) 
prg_ft
ggsave("02. Module scores umap.pdf", plot = prg_ft, width = 14.1, height = 8, units = "in")  


# Our features 

imp_features <- c('Hypoxia', 'AC_like', 'OC_like', 'NPC_like', 'Cycling') 

all_features <- FetchData(object = sc_inew, vars = imp_features) 
all_features <- scale(all_features) %>% as.data.frame()

# Assign cellular / transcriptional states  

prg_features <- all_features %>% dplyr::select(c('AC_like', 'OC_like', 'NPC_like', 'Cycling'))
meta_classes <- colnames(prg_features)[apply(prg_features, 1, which.max)]
meta_classes <- as.data.frame(meta_classes)
meta_classes <- cbind(prg_features, meta_classes)

sc_inew <- AddMetaData(object = sc_inew, metadata = meta_classes$meta_classes, col.name = "metaprog_4_class") 

mp_umap <- DimPlot(sc_inew, 
                   reduction = "umap", 
                   group.by = c("metaprog_4_class"), 
                   cols = c('AC_like' = "#DE1A58", 'OC_like' = "#FF9A00", 'NPC_like' = "#8FA31E", 'Cycling' = "#9112BC"),
                   pt.size = 0.2, 
                   label = T) 
mp_umap 
ggsave("03. Assigned programs umap.pdf", plot = mp_umap, width = 4.1, height = 3, units = "in") 


# ── Plot markers ──────────────────────────────────────────────────────────────  

# Test the module genes for differential statistics  

gs.prog_jc <- lapply(gs.prog, function(x) unique(na.omit(x)))
names(gs.prog_jc) <- paste0("MetaProg", 1:length(gs.prog_jc))

markers_df <- stack(gs.prog_jc) 
colnames(markers_df) <- c("gene", "program")
markers_df <- markers_df[markers_df$gene %in% rownames(sc_inew), ]

# Run differential expression with metamodule genes 

Idents(sc_inew) <- "metaprog_4_class"
DefaultAssay(sc_inew) <- "RNA"

dif_mp.gene <- FindAllMarkers(sc_inew,
                              features = unique(markers_df$gene), 
                              only.pos = TRUE, 
                              logfc.threshold = 0.25,
                              min.pct = 0.1) 

sig_mp.gene <- dif_mp.gene %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(avg_log2FC)) 
sig_mp.gene <- unique(sig_mp.gene$gene) 

# Plot heatmap 

sc_inew$metaprog_4_class <- factor(sc_inew$metaprog_4_class, levels = c('Cycling', 'AC_like', 'NPC_like', 'OC_like'))

DefaultAssay(sc_inew) <- "RNA"
sc_inew <- ScaleData(sc_inew, features = sig_mp.gene)

ht_mp.genes <- DoHeatmap(sc_inew,
                         features = sig_mp.gene,
                         group.by = 'metaprog_4_class',
                         size = 2) + 
  scale_fill_gradientn(colors = c("darkblue", "white", "red"))
ht_mp.genes 
ggsave("04. MP markers heatmap.pdf", plot = ht_mp.genes, width = 8, height = 20, units = "in")  

# Plot dotplot 

picked_mp.genes <- list('AC_like' = c('APOE', 'VIM', 'FOS', 'ID4', 'JUN'), 
                        'OC_like' = c('CNP', 'GPR17', 'SLC1A1', 'ZEB2', 'BIN1'), 
                        'Cycling' = c('SLBP', 'HMGB2', 'NASP', 'TMPO', 'PCNA'),  
                        'NPC_like' = c('SOX11', 'EIF4B', 'CHD7', 'TCF4', 'BAZ2B')) 

sc_inew@active.ident <- factor(sc_inew@active.ident, 
                               levels = c('AC_like', 'OC_like', 'NPC_like', 'Cycling')) 

dp_mp.genes <- DotPlot(sc_inew,
                       features = picked_mp.genes, 
                       cols = c("#FFFF00", "#FF00FF")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dp_mp.genes
ggsave("04. MP markers dotplot.pdf", plot = dp_mp.genes, width = 9, height = 3, units = "in") 


# ── Correlation ───────────────────────────────────────────────────────────────  

cor_features <- all_features %>% dplyr::select(c('Hypoxia', 'AC_like', 'OC_like', 'NPC_like', 'Cycling'))

cor_res <- cor(cor_features, method = "pearson")
corrplot(as.matrix(cor_res), method = "circle", type = "lower", col = colorRampPalette(c("darkblue", "white", "red"))(20))

cor_new <- cor_res %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% c('AC_like', 'OC_like', 'NPC_like', 'Cycling')) %>% 
  select(Hypoxia) %>% 
  rownames_to_column(var = "MPs") 

cor_long <- cor_new %>% pivot_longer(cols= !c("MPs"), names_to = 'Features', values_to = 'correlation')

p.cor_mp <- ggplot(cor_long, aes(x = Features, y = MPs, fill = correlation)) + 
  geom_point(aes(size = correlation), alpha = 1, shape = 21) + 
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p.cor_mp 
ggsave("05. MP corrplot.pdf", plot = p.cor_mp, width = 2.5, height = 3, units = "in")  


# ── Pathway activity inference ──────────────────────────────────────────────── 

sc.g_matrix <- as.matrix(sc_inew@assays$RNA@data) 
progeny_net <- get_progeny(organism = 'human', top = 500) 

# Run mlm

res_progeny <- decoupleR::run_mlm(mat = sc.g_matrix, 
                                  net = progeny_net, 
                                  .source = 'source', 
                                  .target = 'target', 
                                  .mor = 'weight', 
                                  minsize = 5)
res_progeny 
write.csv(res_progeny, file = "Result Pathway activity (mlm).csv", row.names = F) 

# Add as assay and plot pathway features

sc_inew[['pathwaysmlm']] <- res_progeny %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.) 

DefaultAssay(object = sc_inew) <- "pathwaysmlm"

sc_inew <- ScaleData(sc_inew)
sc_inew@assays$pathwaysmlm@data <- sc_inew@assays$pathwaysmlm@scale.data

progeny_umap <- FeaturePlot(sc_inew, 
                            reduction = "umap", 
                            features = c('Androgen', 'EGFR', 'Estrogen', 'Hypoxia', 'JAK-STAT', 'MAPK', 'NFkB', 'PI3K', 'TGFb', 'TNFa', 'Trail', 'WNT', 'p53', 'VEGF'), 
                            ncol = 3, 
                            order = T, 
                            min.cutoff = 0.7,
                            max.cutoff = 2.7,
                            keep.scale = "all") 
progeny_umap
ggsave(progeny_umap, file = "06. Pathway activity umap.pdf", width = 9.45, height = 12.5, units = "in") 

# Data handling and plot heatmap 

Idents(object = sc_inew) <- "metaprog_4_class"

progeny_res_long <- t(as.matrix(sc_inew@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(metaclust = Idents(sc_inew)) %>%
  pivot_longer(cols = -metaclust, names_to = "source", values_to = "score") %>%
  group_by(metaclust, source) %>%
  summarise(mean = mean(score)) 

progeny_res_matx <- progeny_res_long %>%
  pivot_wider(id_cols = 'metaclust', names_from = 'source', values_from = 'mean') %>%
  column_to_rownames('metaclust') %>%
  as.matrix()

# Plot progeny heatmap 

set.seed(123) 
ht_progeny <- Heatmap(progeny_res_matx, 
                      col = colorRampPalette(c("Darkblue", "white","red"))(100), 
                      border = T) 
ht_progeny 
ht_progeny <- as.ggplot(ht_progeny)
ggsave("06. Pathway activity heatmap.pdf", plot = ht_progeny, width = 6.5, height = 2.75, units = "in")  


# ── TFs ───────────────────────────────────────────────────────────────────────  

tfs_net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)

# Run ulm 

res_tfs <- decoupleR::run_ulm(mat = sc.g_matrix, 
                              net = tfs_net, 
                              .source = 'source', 
                              .target = 'target',
                              .mor = 'mor', 
                              minsize = 5) 
res_tfs 
write.csv(res_tfs, file = "Result TF activity (ulm).csv", row.names = F) 

# Add as assay 

sc_inew[['tfsulm']] <- res_tfs %>%
  tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = sc_inew) <- "tfsulm"

sc_inew <- ScaleData(sc_inew) 
sc_inew@assays$tfsulm@data <- sc_inew@assays$tfsulm@scale.data

# Data handling 

n_tfs <- 500
Idents(object = sc_inew) <- "metaprog_4_class"

tfs_res_long <- t(as.matrix(sc_inew@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = Seurat::Idents(sc_inew)) %>%
  tidyr::pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

tfs_res_name <- tfs_res_long %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

tfs_res_matx <- tfs_res_long %>%
  dplyr::filter(source %in% tfs_res_name) %>%
  tidyr::pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix() 

# Check important TFs for plotting selection 

tfs_res_matx %>% t() %>% as_tibble(rownames = "TF") %>% arrange(desc(AC_like)) %>% slice_head(n = 200) %>% pull(TF)

# Plot TFs heatmap (with important TFs labeled) 

imp_tfs_df <- read.xlsx("TFs to plot.xlsx")
tf_ploting <- dput(as.character(imp_tfs_df$Gene)) 

t.tfs_heat <- t(tfs_res_matx) 

set.seed(123) 
mark_at <- which(rownames(t.tfs_heat) %in% tf_ploting)
row_ha1 <- rowAnnotation(foo = anno_mark(at = which(rownames(t.tfs_heat) %in% tf_ploting), 
                                         labels = rownames(t.tfs_heat)[rownames(t.tfs_heat) %in% tf_ploting], 
                                         labels_gp = gpar(fontsize = 10),
                                         side = "right", 
                                         link_width = unit(10, "mm"), 
                                         padding = unit(0.02, "mm"), 
                                         extend = 0.2), 
                         width = unit(5, "cm"))

ht_tfs <- Heatmap(t.tfs_heat, 
                  col = colorRampPalette(c("Darkblue", "white","red"))(100), 
                  na_col = "grey",
                  cluster_rows = T,
                  cluster_columns = T,
                  right_annotation = row_ha1, 
                  row_names_max_width = unit(5,"in"),
                  row_title_gp = gpar(fontsize = 5),
                  column_title_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 12),
                  row_names_gp = gpar(fontsize = 3),
                  show_row_names = F, 
                  show_column_names = T, 
                  border = T)
ht_tfs 
ht_tfs <- as.ggplot(ht_tfs)
ggsave("07. TFs heatmap.pdf", plot = ht_tfs, width = 5, height = 8, units = "in") 


# ── Msigdb enrichment ───────────────────────────────────────────────────────── 

# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html

metaprog_diff.res <- wilcoxauc(sc_inew, 'metaprog_4_class')
write.csv(metaprog_diff.res, file = "Result wilcoxauc metaprog_4_class.csv")

dplyr::count(metaprog_diff.res, group)

# Prepare msigdb gene sets 

msigdb_list <- list(Hallmark = msigdbr(species = "Homo sapiens", category = "H"), 
                    Reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"), 
                    GOBP_now = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), 
                    Wikipath = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS"))

run_geneset <- lapply(msigdb_list, function(df){
  split(df$gene_symbol, df$gs_name)
}) 

# Ranked gene lists for each category 

groups_nm <- unique(metaprog_diff.res$group)

ranks_by_group <- lapply(groups_nm, function(g) {
  metaprog_diff.res %>%
    filter(group == g) %>%
    arrange(desc(auc)) %>%
    select(feature, auc) %>%
    deframe()
})

names(ranks_by_group) <- groups_nm 


# Run fgsea and save all results (min and maxsize are not necessary) 

set.seed(123)
gsea_res_list <- list()

for (gs_cat in names(run_geneset)) {
  for (g in names(ranks_by_group)) {
    gsea_res_list[[length(gsea_res_list) + 1]] <- fgsea(pathways = run_geneset[[gs_cat]], stats = ranks_by_group[[g]], nperm = 10000, minSize = 15, maxSize = 300) %>%
      as_tibble() %>%
      mutate(group = g, geneset_category = gs_cat) 
  }
}

gsea_res_dfs <- bind_rows(gsea_res_list) 
gsea_res_dfs <- gsea_res_dfs %>% mutate(log10_pval = -log10(pval)) 
gsea_res_dfs <- gsea_res_dfs %>% filter(geneset_category %in% c('Hallmark', 'GOBP_now', 'Wikipath'))
gsea_res_dfs <- gsea_res_dfs %>% mutate(leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1)))

write.csv(gsea_res_dfs, file = "Result fgsea (all database, all categories).csv") 


# Plot enrichment heatmaps 

set.seed(123)
groups <- unique(gsea_res_dfs$group)

for (grp in groups) {
  message("plotting: ", grp)
  
  df_top <- gsea_res_dfs %>%
    filter(group == grp) %>%
    group_by(geneset_category) %>%
    slice_max(order_by = NES, n = 5, with_ties = F) %>%
    ungroup()
  
  mat_df <- df_top %>%
    select(pathway, geneset_category, NES) %>%
    pivot_wider(names_from = geneset_category, values_from = NES) %>%
    column_to_rownames("pathway") %>%
    as.matrix() 
  
  row_bar <- rowAnnotation(pval = anno_barplot(df_top$pval, gp = gpar(fill = "orange"), border = FALSE), 
                           annotation_width = unit(2, "cm")) 
  row_ann <- rowAnnotation(geneset = df_top$geneset_category,
                           col = list(geneset = c("Hallmark" = "#e41a1c", "GOBP_now" = "#377eb8", "Reactome" = "#4daf4a", "Wikipath" = "#984ea3")), 
                           border = T) 
  
  ht_gse <- Heatmap(mat_df, 
                    col = viridis::plasma(100),
                    na_col = "white",
                    cluster_rows = F,
                    cluster_columns = F, 
                    border = T,
                    left_annotation = c(row_bar, row_ann), 
                    column_title = paste(grp, "- pathways per geneset")) 
  
  pdf(paste0("08. fgsea ", grp, ".pdf"), width = 7, height = 7)
  draw(ht_gse) 
  dev.off()
}


# ── Plot metaprograms / metadata ────────────────────────────────────────────── 

# Quantification of programs (per cluster / pathology / samples) 

imp_features <- c('Hypoxia', 'AC_like', 'OC_like', 'NPC_like', 'Cycling') 

vl_clust <- VlnPlot(sc_inew, 
                    features = imp_features, 
                    sort = 'decreasing', 
                    group.by = 'seurat_clusters', 
                    pt.size = 0) 
vl_clust
ggsave("09. Programs violin (clusters).pdf", plot = vl_clust, width = 12, height = 6, units = "in")  

vl_sampl <- VlnPlot(sc_inew, 
                    features = imp_features, 
                    sort = 'decreasing', 
                    group.by = 'DataSet', 
                    pt.size = 0) 
vl_sampl
ggsave("09. Programs violin (samples).pdf", plot = vl_sampl, width = 12, height = 8, units = "in") 

vl_types <- VlnPlot(sc_inew, 
                    features = imp_features, 
                    sort = 'decreasing', 
                    group.by = 'Types', 
                    pt.size = 0) 
vl_types
ggsave("09. Programs violin (pathology).pdf", plot = vl_types, width = 6, height = 8, units = "in") 

# Proportion of programs (per cluster / pathology / samples)  

bar_pathl <- dittoBarPlot(sc_inew, "metaprog_4_class", group.by = "Types") 
bar_sampl <- dittoBarPlot(sc_inew, "metaprog_4_class", group.by = "DataSet") 
bar_clust <- dittoBarPlot(sc_inew, "metaprog_4_class", group.by = "seurat_clusters") 

bar_mps <- wrap_plots(bar_pathl, bar_sampl, bar_clust, ncol = 3) 
bar_mps
ggsave("10. Programs barplots.pdf", plot = bar_mps, width = 12, height = 5, units = "in") 


# Save new seurat object 

DefaultAssay(sc_inew) <- "SCT"
saveRDS(sc_inew, file = "mIDH glioma (CCA integrated Seurat 4.3.0).rds") 
