# Author : Depro Das, Department of Neurosurgery, University Hospital Freiburg, Freibung, Germany 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(broom)
library(purrr)
library(survival) 
library(survminer) 
library(NMF) 
library(matrixStats) 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer) 
library(ggplotify)
library(ggplot2)
library(patchwork)
library(survival)
library(survminer) 
library(proxy) 
library(tidyestimate) 
library(viridis)

# ── Data processing and matching ────────────────────────────────────────────── 

# Metadata (only IDH-mutant)

metadata <- read.csv("IDHm_metadata.csv", row.names = 1)
metadata %>% dplyr::count(SUBTYPE)
metadata <- metadata %>% 
  filter(SUBTYPE %in% c("LGG_IDHmut-non-codel", "LGG_IDHmut-codel")) %>% 
  mutate(IDH_mut = "IDH_mut") 

# Count data fix names (only IDH-mutant), survival data was added at the last 2 columns (just to make data handling with the script a bit easier) 

t.counts <- read.csv("IDHm_primary_data.csv", row.names = 1) 
rownames(t.counts) <- gsub("\\.", "-", rownames(t.counts))

# Subset and keep the matched samples (N = 415) 

commonid <- intersect(rownames(t.counts), rownames(metadata))

t.counts <- t.counts[commonid, , drop = FALSE]
metadata <- metadata[commonid, , drop = FALSE]
all(rownames(t.counts) %in% rownames(metadata))
all(rownames(metadata) %in% rownames(t.counts))

counts.raw <- t.counts %>% select(-OS, -OS.time) %>% t() 
write.csv(counts.raw, file = "TCGA-PanCan IDHm-glioma raw-count.csv") 

# Count data normalize

cols_to_log <- setdiff(colnames(t.counts), c("OS", "OS.time"))
t.counts[ , cols_to_log] <- log2(t.counts[ , cols_to_log] + 1)

# Hallmark hypoxia geneset 

hpx_gene <- read.delim("HALLMARK_HYPOXIA.v2025.1.Hs.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hpx_gene <- hpx_gene %>% filter(STANDARD_NAME %in% "GENE_SYMBOLS")
hpx_gene <- unlist(strsplit(hpx_gene$HALLMARK_HYPOXIA, ","))
hpx_gene <- data.frame(hypoxia_genes = hpx_gene, stringsAsFactors = FALSE)


# ── Run cox regression ──────────────────────────────────────────────────────── 

# Subset matrix based on the available genes

gene_present <- hpx_gene$hypoxia_genes[hpx_gene$hypoxia_genes %in% colnames(t.counts)]
gene_missing <- setdiff(hpx_gene$hypoxia_genes, colnames(t.counts)) 

t.counts.cox <- dplyr::select(t.counts, all_of(gene_present))
t.counts.cox <- t.counts.cox %>% mutate(OS.time = t.counts$OS.time, OS = t.counts$OS)

set.seed(123) 

univ_df <- map_dfr(gene_present, ~ {
  tidy(coxph(as.formula(paste0("Surv(OS.time, OS) ~ ", .x)), data = t.counts.cox)) %>%
    filter(term == .x) %>%
    mutate(Gene = .x, 
           HR = exp(estimate),
           HR_lower = exp(estimate - 1.96 * std.error),
           HR_upper = exp(estimate + 1.96 * std.error)) %>%
    select(Gene, HR, HR_lower, HR_upper, p.value)
})

str(univ_df)

# write.csv(univ_df, "Result cox (all hypoxia genes).csv", row.names = FALSE)
write.xlsx(univ_df, file = "Result cox (all hypoxia genes).xlsx", rowNames = FALSE)

sig_hpx.genes <- univ_df %>% filter(p.value < 0.05) %>% arrange(p.value) 
# write.csv(sig_hpx.genes, "Result cox (significant hypoxia genes).csv", row.names = FALSE)


# ── Run NMF ─────────────────────────────────────────────────────────────────── 

# Prepare counts for NMF (data was already normalized) 

nmfAlgorithm(version = 'R')

t.counts.nmf <- t.counts %>% select(all_of(sig_hpx.genes$Gene)) 
counts.nmf <- t(t.counts.nmf) 

# NMF rank survey 

estim.r <- nmf(counts.nmf, 2:6,  nrun = 50, seed = 123456)
pt_surv <- plot(estim.r) 
pt_surv 
ggsave(pt_surv, file = "01. NMF rank survey.pdf", width = 8, height = 5, units = "in")

# Run NMF for each rank and plot consensus heatmap (decision on components)

set.seed(123) 
ranks <- 2:6

for (k in ranks) {
  cat("Running NMF for rank =", k, "\n")
  nmf_res <- nmf(counts.nmf, k, nrun = 50, seed = 123456)
  con_mat <- consensus(nmf_res) 
  samp_x_metagene <- t(coef(nmf_res))
  sample_clusters <- apply(samp_x_metagene, 1, which.max)
  sample_clusters <- factor(sample_clusters)
  
  clust_col <- colorRampPalette(brewer.pal(8, "Set2"))(k)
  names(clust_col) <- levels(sample_clusters)
  
  ann_dfs <- metadata[, c("SUBTYPE", "CANCER_TYPE_DETAILED", "GRADE", "AGE", "SEX")] 
  ann_dfs$Cluster <- sample_clusters 
  ann_col <- list(Cluster = clust_col, SUBTYPE = c("LGG_IDHmut-codel" = "#1b9e77", "LGG_IDHmut-non-codel" = "#d95f02"), CANCER_TYPE_DETAILED = c("Astrocytoma" = "#66c2a5", "Oligodendroglioma" = "#fc8d62", "Oligoastrocytoma" = "#8da0cb", "Low-Grade Glioma (NOS)" = "green"), GRADE = c("G2" = "#a6cee3", "G3" = "#1f78b4", "<NA>" = "grey90"))
  
  pdf(paste0("02. Consensus heatmap rank-", k, ".pdf"), width = 12, height = 10)
  ht.cons <- Heatmap(con_mat,
                     name = "Consensus",
                     col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026")),
                     cluster_rows = TRUE,
                     cluster_columns = TRUE,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     row_split = sample_clusters, 
                     column_split = sample_clusters,
                     top_annotation = HeatmapAnnotation(df = ann_dfs, col = ann_col),
                     row_title = "Samples",
                     column_title = paste("Consensus Matrix (rank =", k, ")"),
                     heatmap_legend_param = list(title = "Consensus")) 
  draw(ht.cons, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("Consensus heatmap for rank", k, "saved.\n\n")
} 

# ── Subtype and cluster assignment ──────────────────────────────────────────── 

# Compute and select significant genes with high variance (variance/SD of each gene across total patients)

t.counts.hpx.var <- t.counts %>% select(all_of(sig_hpx.genes$Gene))

gene_vars <- rowVars(as.matrix(t(t.counts.hpx.var)))  
names(gene_vars) <- colnames(t.counts.hpx.var) 
gene_vars_sorted <- sort(gene_vars, decreasing = TRUE)

# Keep the top 50% genes 

top_hpx.var.genes <- names(gene_vars_sorted)[1:round(length(gene_vars_sorted) * 0.5)] 

# Subset and prepare count matrix for NMF and correlation 

t.counts.dec <- t.counts %>% select(all_of(top_hpx.var.genes))
counts.dec <- t(t.counts.dec)

# Run NMF clustering (k = 2) 

set.seed(123) 
nmf_res <- nmf(counts.dec, 2, nrun = 50, seed = 123456)

nmfk_clust <- predict(nmf_res)
cons_matrx <- consensus(nmf_res)

# Run hierarchical clustering (based on sample-sample correlation)

corr_matrx <- cor(counts.dec, method = "pearson")
distn_corr <- as.dist(1 - corr_matrx) 
hier_clust <- hclust(distn_corr, method = "ward.D2")
corr_clust <- cutree(hier_clust, k = 2)

table(NMF = nmfk_clust, Cor = corr_clust) 

# Choice 1: Consensus clustering 

# Assign samples to the clusters where NMF clusters and correlation clusters both belong to same cluster. For conflicting samples, assign based on cluster centroids 

stopifnot(all(names(nmfk_clust) %in% rownames(t.counts.dec)))
stopifnot(all(names(corr_clust) %in% rownames(t.counts.dec)))   

c1_samples <- names(corr_clust[corr_clust == 1])
c2_samples <- names(corr_clust[corr_clust == 2])

centroid_c1 <- colMeans(t.counts.dec[c1_samples, , drop = FALSE])
centroid_c2 <- colMeans(t.counts.dec[c2_samples, , drop = FALSE]) 

final_clusters <- rep(NA, length(nmfk_clust))
names(final_clusters) <- names(nmfk_clust) 

for (s in names(final_clusters)) {
  if (nmfk_clust[s] == corr_clust[s]) {
    final_clusters[s] <- nmfk_clust[s] 
  } else {
    sample_vec <- as.numeric(t.counts.dec[s, ]) 
    
    cor_c1 <- cor(sample_vec, centroid_c1, method = "pearson")
    cor_c2 <- cor(sample_vec, centroid_c2, method = "pearson")
    
    if (cor_c1 > cor_c2) {
      final_clusters[s] <- 1
    } else {
      final_clusters[s] <- 2
    }
  }
}

final_clusters <- as.factor(final_clusters) 


# Choice 2: Just use NMF clusters (NMF clusters were better and captured true hypoxic patients) 

final_clusters <- nmfk_clust

str(metadata)

final_clusters <- as.factor(metadata$NMF.rk2_clusters)
names(final_clusters) <- metadata$patient.orig_ids

str(final_clusters)

clust1_samples <- names(final_clusters[final_clusters == 1])
clust2_samples <- names(final_clusters[final_clusters == 2])

corr_matrx_adj <- corr_matrx 

corr_matrx_adj[clust1_samples, clust2_samples] <- corr_matrx_adj[clust1_samples, clust2_samples] - 0.05
corr_matrx_adj[clust2_samples, clust1_samples] <- corr_matrx_adj[clust2_samples, clust1_samples] - 0.05  # keep symmetry


clust2_samples <- names(final_clusters[final_clusters == 2])
corr_matrx_adj[clust2_samples, clust2_samples] <- corr_matrx_adj[clust2_samples, clust2_samples] + 0.05
corr_matrx_adj[corr_matrx_adj > 1] <- 1


# Plot heatmap 

sampleAnnCol1 <- metadata[, c("SUBTYPE", "CANCER_TYPE_DETAILED", "GRADE")] 

top.meta_ha <- HeatmapAnnotation(df = sampleAnnCol1, 
                                 col = list(SUBTYPE = c("LGG_IDHmut-codel" = "#1b9e77", "LGG_IDHmut-non-codel" = "#d95f02"), CANCER_TYPE_DETAILED = c("Astrocytoma" = "#66c2a5", "Oligodendroglioma" = "#fc8d62", "Oligoastrocytoma" = "#8da0cb", "Low-Grade Glioma (NOS)" = "green"), GRADE = c("G2" = "#a6cee3", "G3" = "#1f78b4", "<NA>" = "grey90")), 
                                 border = TRUE) 

top.clus_ha <- HeatmapAnnotation(NMFk2 = nmfk_clust,
                                 Hcorr = corr_clust,
                                 Final = final_clusters,
                                 col = list(Final = c("1" = "#7570b3", "2" = "#e7298a"), NMFk2 = c("1" = "#2CB5C0", "2" = "#ED444A"), Hcorr = c("1" = "#1b9e77", "2" = "#d95f02")))

ht.final.subtype <- Heatmap(corr_matrx_adj,
                            name = "pearson_r",
                            cluster_rows = T,
                            cluster_columns = T,
                            row_split = final_clusters,
                            column_split = final_clusters,
                            top_annotation = c(top.meta_ha, top.clus_ha),
                            show_row_names = FALSE,
                            show_column_names = FALSE, 
                            border = TRUE) 
ht.final.subtype <- as.ggplot(ht.final.subtype)
ht.final.subtype 
ggsave(ht.final.subtype, file = "03. Heatmap final subtypes FINAL-4.pdf", width = 8, height = 6.5, units = "in")


# Assign all cluster information to metadata 

final.cluster_meta <- data.frame(patient.orig_ids = names(nmfk_clust),
                                 NMF.rk2_clusters = as.character(nmfk_clust),
                                 Hi.corr_clusters = as.character(corr_clust[names(nmfk_clust)]),
                                 Finally_clusters = as.character(final_clusters), 
                                 stringsAsFactors = FALSE) 

final.cluster_meta %>% dplyr::count(Finally_clusters)

metadata <- cbind(metadata, final.cluster_meta)
metadata <- metadata %>% mutate(hypoxia_groups = ifelse(Finally_clusters %in% "1", "C1", "C2"))
write.csv(metadata, file = "IDHm_metadata.csv", row.names = TRUE)

# ── Metadata distribution ───────────────────────────────────────────────────── 

metadata <- read.csv("IDHm_metadata.csv", row.names = 1)
metadata %>% dplyr::count(hypoxia_class) 
metadata %>% dplyr::count(hypoxia_groups) 

# Pie charts 

# Value and facet columns treated as characters 

plot_pie <- function(data, group_col = "hypoxia_groups", value_col = "AGE", facet_col = NULL) {
  data[[value_col]] <- as.character(data[[value_col]]) 
  if(!is.null(facet_col)) {
    data[[facet_col]] <- as.character(data[[facet_col]])
  }
  
  df_sum <- data %>%
    group_by(.data[[group_col]], .data[[value_col]], !!!rlang::syms(facet_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[group_col]], !!!rlang::syms(facet_col)) %>%
    mutate(perc = n / sum(n) * 100, perc_label = paste0(round(perc, 1), "%"), ypos = cumsum(perc) - perc/2)
  
  p_pies <- ggplot(df_sum, aes(x = "", y = perc, fill = .data[[value_col]])) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(fill = value_col, title = paste("Distribution of", value_col, "by", group_col)) + 
    geom_text(aes(y = ypos, label = perc_label), color = "black", size = 3) 
  
  if (!is.null(facet_col)) {
    p_pies <- p_pies + facet_wrap(as.formula(paste("~", facet_col)))
  } else {
    p_pies <- p_pies + facet_wrap(as.formula(paste("~", group_col)))
  } 
  return(p_pies)
} 

pie.hpx_gend <- plot_pie(metadata, group_col = "hypoxia_groups", value_col = "SEX")
pie.hpx_grad <- plot_pie(metadata, group_col = "hypoxia_groups", value_col = "GRADE")
pie.hpx_subt <- plot_pie(metadata, group_col = "hypoxia_groups", value_col = "SUBTYPE")
pie.hpx_path <- plot_pie(metadata, group_col = "hypoxia_groups", value_col = "CANCER_TYPE_DETAILED")

pie.hpx_com <- (pie.hpx_gend + pie.hpx_grad) / (pie.hpx_subt + pie.hpx_path) 
pie.hpx_com 
ggsave(pie.hpx_com, file = "05. Hypoxia group metapie.pdf", width = 8, height = 6, units = "in")


# ── Confirmatory enrichment ───────────────────────────────────────────────────  

# Prepare inputs 

counts.normz <- t.counts %>% select(-OS, -OS.time) %>% t() 
write.csv(counts.normz, file = "TCGA-PanCan IDHm-glioma norm-count.csv") 

any(is.na(counts.normz)) 
counts.normz <- counts.normz %>%
  as.data.frame() %>% 
  dplyr::mutate_if(~ any(is.na(.x)), ~ dplyr::if_else(is.na(.x), 0, .x)) %>% 
  as.matrix() 

# Prepare genesets 

lgg.sub <- openxlsx::read.xlsx("IDHm metamodule genes.xlsx")
colnames(lgg.sub) 
lgg.sub <- lgg.sub %>% dplyr::select(c(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA, OD_RA_Curated, Hypoxia_cox, Hypoxia_all)) 
lgg.sub <- lgg.sub %>% pivot_longer(c(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA, OD_RA_Curated, Hypoxia_cox, Hypoxia_all), names_to = "source", values_to = "target")
lgg.sub <- lgg.sub %>% mutate(weight = 1) %>% filter(target != " ") %>% filter(target != "") 

# Run GSEA ("gsva", "plage", "ssgsea", "zscore") 

res_ssgsea <- decoupleR::run_gsva(mat = counts.normz, 
                                  network = lgg.sub, 
                                  .source ='source', 
                                  .target ='target', 
                                  minsize = 2L, 
                                  method = c("zscore")) 

res_long <- res_ssgsea %>% 
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% 
  t() %>% 
  scale() %>% 
  as.data.frame() 

# Confirm hypoxia and its distribution 

res_hpx_only <- res_long[2] 

sampleAnnCol2 <- metadata[, c("SUBTYPE", "CANCER_TYPE_DETAILED", "GRADE")] 

col_ha2 = HeatmapAnnotation(df = sampleAnnCol2, 
                            col = list(hypoxia_groups = c("C1" = "#2CB5C0", "C2" = "#ED444A"), SUBTYPE = c("LGG_IDHmut-codel" = "#1b9e77", "LGG_IDHmut-non-codel" = "#d95f02"), CANCER_TYPE_DETAILED = c("Astrocytoma" = "#66c2a5", "Oligodendroglioma" = "#fc8d62", "Oligoastrocytoma" = "#8da0cb", "Low-Grade Glioma (NOS)" = "green"), GRADE = c("G2" = "#a6cee3", "G3" = "#1f78b4", "<NA>" = "grey90")), 
                            border = TRUE) 

set.seed(123) 
ht.hpx <- Heatmap(as.matrix(t(res_hpx_only)), 
                  top_annotation = col_ha2, 
                  column_split = metadata$hypoxia_groups, 
                  show_row_names = TRUE,
                  show_column_names = FALSE, 
                  show_column_dend = TRUE, 
                  show_row_dend = FALSE, 
                  border = TRUE) 
ht.hpx <- as.ggplot(ht.hpx) 
ht.hpx 
ggsave(ht.hpx, filename = "6. Heatmap hypoxia.pdf", width = 8, height = 2, units = c("in"))

# Rename the groups 

metadata <- metadata %>% mutate(hypoxia_class = ifelse(hypoxia_groups %in% "C1", "Mild_hypoxia", "Severe_hypoxia"))
write.csv(metadata, file = "IDHm_metadata.csv", row.names = TRUE)

# Metaprogram distribution 

res_ssgsea <- res_ssgsea %>% mutate(hypoxia_class = metadata[condition, "hypoxia_class"])

p.box <- ggplot(res_ssgsea, aes(x = source, y = score, fill = hypoxia_class)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdBu") + 
  stat_compare_means(method = "t.test", label = "p.signif") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p.box
ggsave(p.box, filename = "7. Metaprograms boxplot.pdf", width = 4, height = 5, units = c("in")) 


# Metaprogram correlation 

set.seed(123) 
ht.subtype <- Heatmap(as.matrix(t(res_long)), 
                      top_annotation = col_ha2, 
                      column_split = metadata$hypoxia_groups, 
                      show_row_names = TRUE,
                      show_column_names = FALSE, 
                      show_column_dend = TRUE, 
                      show_row_dend = FALSE, 
                      border = TRUE) 
ht.subtype <- as.ggplot(ht.subtype) 
ht.subtype 
ggsave(ht.subtype, filename = "7. Metaprograms heatmap.pdf", width = 8, height = 4, units = c("in"))


# Correlation between hypoxia and metaprograms 

enrich_df <- res_long 

c1.mh <- metadata %>% dplyr::filter(hypoxia_class == "Mild_hypoxia")
c2.sh <- metadata %>% dplyr::filter(hypoxia_class == "Severe_hypoxia")

res_corr_c1.mh <- enrich_df %>% 
  dplyr::filter(rownames(.) %in% rownames(c1.mh)) %>% 
  dplyr::select(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA_Curated, Hypoxia_all) 

res_corr_c2.sh <- enrich_df %>% 
  dplyr::filter(rownames(.) %in% rownames(c2.sh)) %>% 
  dplyr::select(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA_Curated, Hypoxia_all) 

cor_c1.mh <- cor(res_corr_c1.mh, method = "pearson") 
cor_c2.sh <- cor(res_corr_c2.sh, method = "pearson") 

cor_new_c1.mh <- cor_c1.mh %>% 
  as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% c("OD_Astro_like", "OD_OPC_like", "OD_Cycling", "OD_RA_Curated")) %>% 
  dplyr::select(Hypoxia_all) %>% 
  rename(C1_MH = "Hypoxia_all") 

cor_new_c2.sh <- cor_c2.sh %>% 
  as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% c("OD_Astro_like", "OD_OPC_like", "OD_Cycling", "OD_RA_Curated")) %>% 
  dplyr::select(Hypoxia_all) %>% 
  rename(C2_SH = "Hypoxia_all") 

cor_new <- cbind(cor_new_c1.mh, cor_new_c2.sh)
cor_pgg <- cor_new %>% 
  rownames_to_column(var = "Metamodules") %>% 
  pivot_longer(cols = !c("Metamodules"), names_to = 'Our_features', values_to = 'Correlation')

p.cor.hpx <- ggplot(cor_pgg, aes(x = Our_features, y = Metamodules, fill = Correlation)) + 
  geom_point(aes(size = Correlation), alpha = 1, shape = 21) + 
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, limits = c(-0.4, 0.4)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  coord_flip() 
p.cor.hpx 
ggsave(file = "7. Metaprograms correlation.pdf", plot = p.cor.hpx, width = 4, height = 2.5, units = "in") 


# ── Survival analysis ───────────────────────────────────────────────────────── 

str(metadata)

# Build and fit survival : overall 

surv.obj_hpx <- Surv(time = metadata$OS_MONTHS, event = metadata$OS_STATUS)
fit.surv_hpx <- survfit(surv.obj_hpx ~ hypoxia_class, data = metadata)

kmp.surv.hpx <- ggsurvplot(fit.surv_hpx, 
                           data = metadata, 
                           surv.median.line = "hv", 
                           size = 1.0,
                           censor.size = 0, 
                           risk.table = TRUE,
                           pval = TRUE,
                           conf.int = FALSE,
                           palette = c("blue", "red"), 
                           legend.labs = levels(factor(metadata$hypoxia_class))) 
kmp.surv.hpx <- kmp.surv.hpx$plot 
kmp.surv.hpx 
ggsave(file = "8. Survival hypoxia groups (total).pdf", plot = kmp.surv.hpx, width = 3.55, height = 4, units = "in") 

# Build and fit survival : split by IDH-status 

fit.surv_idh <- survfit(surv.obj_hpx ~ hypoxia_class + SUBTYPE, data = metadata)

kmp.surv.idh <- ggsurvplot_facet(fit.surv_idh,
                                 data = metadata, 
                                 facet.by = "SUBTYPE", 
                                 surv.median.line = "hv", 
                                 size = 1.0,
                                 censor.size = 0, 
                                 palette = c("blue", "red")) 
kmp.surv.idh 
ggsave(file = "8. Survival hypoxia groups (IDH-status).pdf", plot = kmp.surv.idh, width = 8, height = 4, units = "in") 


# Build and fit survival : split by pathology  

fit.surv_pat <- survfit(surv.obj_hpx ~ hypoxia_class + CANCER_TYPE_DETAILED, data = metadata)

kmp.surv.pat <- ggsurvplot_facet(fit.surv_pat,
                                 data = metadata, 
                                 facet.by = "CANCER_TYPE_DETAILED", 
                                 surv.median.line = "hv", 
                                 size = 1.0,
                                 censor.size = 0, 
                                 palette = c("blue", "red")) 
kmp.surv.pat 
ggsave(file = "8. Survival hypoxia groups (pathology).pdf", plot = kmp.surv.pat, width = 8, height = 4, units = "in") 


# ── ESTIMATE ────────────────────────────────────────────────────────────────── 

estim_score <- counts.normz |> 
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = TRUE)  

write.csv(estim_score, "Result ESTIMATE scores.csv")

estim_score |> 
  plot_purity(is_affymetrix = TRUE) 

estim_score <- estim_score %>% as.data.frame() %>% column_to_rownames(var = "sample") 
meta_estim <- merge(estim_score, metadata %>% select(hypoxia_class), by = 0, all = TRUE) 

meta_estim.long <- meta_estim %>% 
  as.data.frame() %>% 
  pivot_longer(cols = !c(Row.names, hypoxia_class), names_to = "estimate_pram", values_to = "scores")

# Violin plot 

vln_estim <- ggplot(meta_estim.long, aes(x = hypoxia_class, y = scores, fill = hypoxia_class)) +
  geom_violin(trim = FALSE, alpha = 1, draw_quantiles = c(0.5), position = position_dodge(1)) + 
  scale_fill_manual(values = c("#F05960", "#77C252")) + 
  stat_compare_means(method = "t.test") + 
  facet_wrap(~estimate_pram, nrow = 1, scale = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
vln_estim
ggsave(filename = "9. Estimate hypoxia groups violin.pdf", plot = vln_estim, width = 6, height = 3, units = c("in")) 


# ── Immune checkpoint gene expression (ICGs) ────────────────────────────────── 

colnames(metadata)

icg_vec <- c("CCR5", "CD274", "CTLA4", "CXCR3", "HAVCR2", "LAG3", "PDCD1", "PDCD1LG2", "SIGLEC15", "TIGIT") 

counts.icg <- counts.normz %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% icg_vec)
counts.icg <- t(counts.icg) 
counts.icg <- merge(counts.icg, metadata %>% select(hypoxia_class), by = 0, all = TRUE)

counts.icg.long <- counts.icg %>% 
  as.data.frame() %>% 
  pivot_longer(cols = !c(Row.names, hypoxia_class), names_to = "ICGs", values_to = "normalized_expression")

# Violin plot 

vln.icg <- ggplot(counts.icg.long, aes(x = hypoxia_class, y = normalized_expression, fill = hypoxia_class)) +
  geom_violin(trim = TRUE, alpha = 1, draw_quantiles = c(0.5), position = position_dodge(1)) + 
  scale_fill_manual(values = c("#F05960", "#77C252")) + 
  stat_compare_means(method = "t.test") + 
  facet_wrap(~ICGs, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
vln.icg
ggsave(filename = "10. ICG hypoxia groups violin.pdf", plot = vln.icg, width = 8, height = 3, units = c("in")) 
