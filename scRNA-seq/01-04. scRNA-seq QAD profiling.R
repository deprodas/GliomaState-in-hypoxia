# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(Seurat) 
library(patchwork) 
library(ComplexHeatmap) 
library(circlize) 
library(ggplotify)
library(ggrepel) 
library(presto)
library(msigdbr) 
library(fgsea) 

# ── Data input ────────────────────────────────────────────────────────────────  

sc_inew <- readRDS("mIDH glioma (CCA integrated Seurat 4.3.0).rds") 
gs.qads <- read.xlsx("GBM_QAD_ptalign_geneset.xlsx") 

DimPlot(sc_inew, reduction = "umap", label = T) 

# ── Module scoring ────────────────────────────────────────────────────────────  

# Prepare genesets and run module scoring 

gs.comp <- list(Q_stage = as.character(gs.qads$Q_stage), 
                A_stage = as.character(gs.qads$A_stage), 
                D_stage = as.character(gs.qads$D_stage)) 

gs.list <- lapply(gs.comp, function(x) list(as.character(x))) 

DefaultAssay(sc_inew) <- "SCT"

for (nm in names(gs.list)) {
  sc_inew <- AddModuleScore(sc_inew, features = gs.list[[nm]], name = nm) 
  oldname <- grep(paste0("^", nm), colnames(sc_inew@meta.data), value = TRUE)
  colnames(sc_inew@meta.data)[colnames(sc_inew@meta.data) == oldname] <- nm
}

# Umap enrichment of programs 

Qt_ft <- FeaturePlot(sc_inew, 
                     reduction = "umap", 
                     min.cutoff = 0.0,
                     max.cutoff = 0.8,
                     features = c('Q_stage')) 
Qt_ft 

At_ft <- FeaturePlot(sc_inew, 
                     reduction = "umap", 
                     min.cutoff = 0.0,
                     max.cutoff = 0.15,
                     features = c('A_stage')) 
At_ft 

Dt_ft <- FeaturePlot(sc_inew, 
                     reduction = "umap", 
                     min.cutoff = 0.0,
                     max.cutoff = 0.2,
                     features = c('D_stage')) 
Dt_ft 

prg_ft <- wrap_plots(Qt_ft, At_ft, Dt_ft, ncol = 3) 
prg_ft
ggsave("01. Module scores umap.pdf", plot = prg_ft, width = 12.25, height = 3.5, units = "in")  


# Quantification of QAD-stages (per cluster / pathology / samples) 

vln_features <- c('Q_stage', 'A_stage', 'D_stage') 

vl_metap <- VlnPlot(sc_inew, 
                    features = vln_features, 
                    sort = 'decreasing', 
                    group.by = 'metaprog_4_class', 
                    pt.size = 0) 
vl_metap
ggsave("02. Programs violin (mps).pdf", plot = vl_metap, width = 12, height = 4, units = "in")  

vl_clust <- VlnPlot(sc_inew, 
                    features = vln_features, 
                    sort = 'decreasing', 
                    group.by = 'seurat_clusters', 
                    pt.size = 0) 
vl_clust
ggsave("02. Programs violin (clusters).pdf", plot = vl_clust, width = 16, height = 4, units = "in")  

vl_sampl <- VlnPlot(sc_inew, 
                    features = vln_features, 
                    sort = 'decreasing', 
                    group.by = 'DataSet', 
                    pt.size = 0) 
vl_sampl
ggsave("02. Programs violin (samples).pdf", plot = vl_sampl, width = 12, height = 4, units = "in") 

vl_types <- VlnPlot(sc_inew, 
                    features = vln_features, 
                    sort = 'decreasing', 
                    group.by = 'Types', 
                    pt.size = 0) 
vl_types
ggsave("02. Programs violin (pathology).pdf", plot = vl_types, width = 5, height = 4, units = "in") 


# ── Enrichment and correlation ──────────────────────────────────────────────── 

# Our features 

imp_features <- c('Q_stage', 'A_stage', 'D_stage', 'AC_like', 'OC_like', 'NPC_like', 'Cycling', 'metaprog_4_class')
all_features <- FetchData(object = sc_inew, vars = imp_features)

# Plot enrichment heatmap  

avg_features <- all_features %>%
  group_by(metaprog_4_class) %>%
  summarise(Q_stage = mean(Q_stage, na.rm = TRUE), 
            A_stage = mean(A_stage, na.rm = TRUE), 
            D_stage = mean(D_stage, na.rm = TRUE)) %>%
  ungroup()

ht.matix_qad <- avg_features %>%
  column_to_rownames("metaprog_4_class") %>%
  as.matrix() %>% 
  t() %>% 
  scale() 

set.seed(123) 
ht.enr_qad <- Heatmap(as.matrix(ht.matix_qad), 
                      col = colorRampPalette(c("Darkblue", "white","red"))(100), 
                      border = T) 
ht.enr_qad
ht.enr_qad <- as.ggplot(ht.enr_qad)
# ggsave("03. Enrichment heatmap.pdf", plot = ht.enr_qad, width = 3.6, height = 2.5, units = "in")  


# Correlation 

num_features <- all_features %>% dplyr::select(c('Q_stage', 'A_stage', 'D_stage', 'AC_like', 'OC_like', 'NPC_like', 'Cycling'))
num_features <- scale(num_features) %>% as.data.frame() 

cor_res <- cor(num_features, method = "pearson")

cor_new <- cor_res %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% c('AC_like', 'OC_like', 'NPC_like', 'Cycling')) %>% 
  select(Q_stage, A_stage, D_stage) 

set.seed(123) 
ht.cor_qad <- Heatmap(as.matrix(cor_new), 
                      col = colorRampPalette(c("Darkblue", "white","red"))(100), 
                      border = T) 
ht.cor_qad
ht.cor_qad <- as.ggplot(ht.cor_qad)
# ggsave("04. Correlation heatmap.pdf", plot = ht.cor_qad, width = 3, height = 2.5, units = "in")  


# ── Classify again (QAD stages) ─────────────────────────────────────────────── 

# Assign only QAD stages  

qad_features <- all_features %>% dplyr::select(c('Q_stage', 'A_stage', 'D_stage')) 
qad_assigned <- colnames(qad_features)[apply(qad_features, 1, which.max)]
qad_assigned <- as.data.frame(qad_assigned)
qad_assigned <- cbind(qad_features, qad_assigned)

sc_inew <- AddMetaData(object = sc_inew, metadata = qad_assigned$qad_assigned, col.name = "QAD_class") 

# QAD cell proportion 

qad_label <- cbind(qad_assigned, all_features) 
propr_qad <- qad_label %>% select(c(metaprog_4_class, qad_assigned)) 

qad_perct <- propr_qad %>% 
  group_by(metaprog_4_class, qad_assigned) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(percent = 100 * n / sum(n)) %>% 
  ungroup() 

p.pie_qad <- ggplot(qad_perct, aes(x = "", y = percent, fill = qad_assigned)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ metaprog_4_class) +
  theme_void() +
  geom_text(aes(label = paste0(round(percent, 1), "%")), position = position_stack(vjust = 0.5), size = 4) 
p.pie_qad
ggsave("05. QAD proportion pie.pdf", plot = p.pie_qad, width = 10, height = 6, units = "in")  


# Assign final classification (AC-like sub-classification/ Q-, A-, D-stages) 

colnames(qad_label) <- make.unique(colnames(qad_label))

qad_label <- qad_label %>% 
  mutate(metaprog_QAD = if_else(metaprog_4_class == "AC_like" & qad_assigned == "Q_stage", "AC_Q", 
                        if_else(metaprog_4_class == "AC_like" & qad_assigned == "A_stage", "AC_A", 
                        if_else(metaprog_4_class == "AC_like" & qad_assigned == "D_stage", "AC_D", metaprog_4_class))))

sc_inew <- AddMetaData(sc_inew, metadata = qad_label$metaprog_QAD, col.name = "metaprog_QAD")


# Plot classifications 

qad_umap <- DimPlot(sc_inew, 
                    reduction = "umap", 
                    group.by = c("QAD_class", "metaprog_QAD"), 
                    cols = c('Q_stage' = "#DE1A58", 'A_stage' = "#134686", 'D_stage' = "#FCC61D", 'AC_Q' = "#FF5555", 'AC_A' = "#8FA31E", 'AC_D' = "#5D2F77", 'NPC_like' = "#3396D3", 'OC_like' = "#FCC61D", 'Cycling' = "#C47BE4"))  
qad_umap
ggsave("06. Assigned QAD umap.pdf", plot = qad_umap, width = 7.9, height = 3, units = "in")  


# ── Previous study enrichment ───────────────────────────────────────────────── 

# Enrichment on total glioma population 

gs_pubgsc.gene <- read.xlsx("Published studies genesets.xlsx") 
available_gene <- rownames(sc_inew)

gs_pubgsc.list <- lapply(gs_pubgsc.gene, function(x) {
  unique(na.omit(as.character(x)))
}) 

gs_pubgsc.keep <- lapply(gs_pubgsc.list, function(gs) {
  intersect(gs, available_gene)
})

DefaultAssay(sc_inew) <- "SCT"

for (nm in names(gs_pubgsc.keep)) {
  sc_inew <- AddModuleScore(object = sc_inew, features = list(gs_pubgsc.keep[[nm]]), name = nm) 
  oldname <- grep(paste0("^", nm), colnames(sc_inew@meta.data), value = TRUE)
  colnames(sc_inew@meta.data)[colnames(sc_inew@meta.data) == oldname] <- nm
}

feature_pubgsc <- names(gs_pubgsc.keep)

# Plot umaps 

pubgsc_umap <- FeaturePlot(sc_inew, 
                           reduction = "umap", 
                           features = feature_pubgsc, 
                           ncol = 4, 
                           order = T, 
                           min.cutoff = 0.01, 
                           keep.scale = "feature")
pubgsc_umap 
ggsave("07. Published signature umap.pdf", plot = pubgsc_umap, width = 13.25, height = 8, units = "in")  


# Enrichment on AC-like sub-classification/ Q-, A-, D-stages 

# Data handling 

qad.mp_namings <- c(feature_pubgsc, 'metaprog_QAD')
qad.mp_feature <- FetchData(object = sc_inew, vars = qad.mp_namings)
qad.mp_feature <- qad.mp_feature %>% filter(metaprog_QAD %in% c('AC_Q', 'AC_A', 'AC_D')) 

num_col_qad.mp <- names(qad.mp_feature)[sapply(qad.mp_feature, is.numeric)]
cat_col_qad.mp <- "metaprog_QAD"

avgscor_qad.mp <- qad.mp_feature %>%
  group_by(!!sym(cat_col_qad.mp)) %>%
  summarise(across(all_of(num_col_qad.mp), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames(var = cat_col_qad.mp) %>% 
  scale()

# Plot heatmap 

pub.study_meta <- read.xlsx("Published studies metadata.xlsx")
col.ann_qad.mp <- HeatmapAnnotation(Studies = dput(as.character(pub.study_meta$Studies)), 
                                    col = list(Studies = c('Whitfield_et_al_2002' = '#FF6F00', 'Zhang_et_al_2014' = 'royalblue', 'Pollen_et_al_2015' = 'limegreen', 'Xie_et_al_2020' = 'gold','Cheung_et_al_2013' = '#40E0D0','Xie_et_al_2022' = 'deeppink1')),
                                    border = TRUE)

ht_qad.mp <- Heatmap(as.matrix(avgscor_qad.mp), 
                     top_annotation = col.ann_qad.mp, 
                     col = colorRampPalette(c("Darkblue", "white","red"))(100), 
                     column_split = pub.study_meta$Studies,
                     border = TRUE) 
ht_qad.mp
ht_qad.mp <- as.ggplot(ht_qad.mp)
ggsave("08. Published signature heatmap.pdf", plot = ht_qad.mp, width = 7.75, height = 4.5, units = "in")


# ── Msigdb enrichment ───────────────────────────────────────────────────────── 

# Enrichment on AC-like sub-classification/ Q-, A-, D-stages 

DefaultAssay(sc_ac.pop) <- "SCT"
sc_ac.pop <- subset(sc_inew, subset = metaprog_4_class == "AC_like") 

difres_ac <- wilcoxauc(sc_ac.pop, 'metaprog_QAD')
write.csv(difres_ac, file = "Result wilcoxauc MP-QAD.csv")

# Prepare msigdb gene sets 

msigdb_list <- list(Hallmark = msigdbr(species = "Homo sapiens", category = "H"), 
                    Reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"), 
                    GOBP_now = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), 
                    Wikipath = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS"))

run_geneset <- lapply(msigdb_list, function(df){
  split(df$gene_symbol, df$gs_name)
}) 

# Ranked gene lists for each category

groups_nm <- unique(difres_ac$group)

ranks_by_group <- lapply(groups_nm, function(g) {
  difres_ac %>%
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
gsea_res_dfs <- gsea_res_dfs %>% mutate(leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1)))

write.csv(gsea_res_dfs, file = "Result fgsea MP-QAD (all database, all categories).csv") 


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
  
  row_ha <- rowAnnotation(geneset = df_top$geneset_category,
                          col = list(geneset = c("Hallmark" = "#e41a1c", "GOBP_now" = "#377eb8", "Reactome" = "#4daf4a", "Wikipath" = "#984ea3")), 
                          border = T) 
  
  ht_gse <- Heatmap(mat_df, 
                    col = viridis::plasma(100),
                    na_col = "white",
                    cluster_rows = F,
                    cluster_columns = F, 
                    border = T,
                    left_annotation = row_ha,
                    column_title = paste(grp, "- pathways per geneset")) 
  
  pdf(paste0("09. fgsea ", grp, ".pdf"), width = 5, height = 7)
  draw(ht_gse) 
  dev.off()
} 

# Save new seurat object 

DefaultAssay(sc_inew) <- "SCT"
saveRDS(sc_inew, file = "mIDH glioma (CCA integrated Seurat 4.3.0).rds") 
