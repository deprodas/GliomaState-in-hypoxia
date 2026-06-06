# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

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

# Treatment information 

treatment <- read.csv("data_timeline_treatment.csv")

treatment <- treatment %>% 
  group_by(PATIENT_ID) %>%
  summarise(TREATMENT_TYPE = paste(unique(TREATMENT_TYPE), collapse = "_"))

treatment <- treatment %>% filter(treatment$PATIENT_ID %in% rownames(metadata))
treatment %>% dplyr::count(TREATMENT_TYPE) %>% print(n = 25)

small_trt <- treatment %>%
  dplyr::count(TREATMENT_TYPE) %>% 
  filter(n < 10) %>% pull(TREATMENT_TYPE)

treatment <- treatment %>%
  mutate(TREATMENT_new = ifelse(TREATMENT_TYPE %in% small_trt, "Other treatments", TREATMENT_TYPE))

treatment %>%
  dplyr::count(TREATMENT_new)

metadata <- metadata %>%
  rownames_to_column(var = 'PATIENT_ID') %>% 
  left_join(treatment, by = "PATIENT_ID") %>%
  mutate(TREATMENT_new = ifelse(is.na(TREATMENT_new), "not_available", TREATMENT_new)) %>% 
  column_to_rownames(var = 'PATIENT_ID') 

write.csv(metadata, file = "IDHm_metadata.csv")

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
write.csv(univ_df, "Result cox (all hypoxia genes).csv", row.names = FALSE)

sig_hpx.genes <- univ_df %>% filter(p.value < 0.05) %>% arrange(p.value) 
write.csv(sig_hpx.genes, "Result cox (significant hypoxia genes).csv", row.names = FALSE)

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
metadata$SUBTYPE_GRADE <- interaction(metadata$SUBTYPE, metadata$GRADE)
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
    mutate(perc = n / sum(n) * 100, perc_label = paste0("N=", n, "\n(", round(perc, 1), "%)"), ypos = cumsum(perc) - perc/2)
  
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
pie.hpx_tret <- plot_pie(metadata, group_col = "hypoxia_groups", value_col = "TREATMENT_new")
pie.hpx_path <- plot_pie(metadata, group_col = "hypoxia_groups", value_col = "CANCER_TYPE_DETAILED")

pie.hpx_com <- (pie.hpx_gend + pie.hpx_grad) / (pie.hpx_subt + pie.hpx_path) / pie.hpx_tret
pie.hpx_com 
ggsave(pie.hpx_com, file = "05. Hypoxia group metapie.pdf", width = 14, height = 10, units = "in")

# Focusing on other treatment groups

other_treat <- metadata %>% dplyr::filter(TREATMENT_new %in% c("Other treatments"))

pie.oth_trt <- plot_pie(other_treat, group_col = "hypoxia_groups", value_col = "TREATMENT_TYPE")
pie.oth_trt
ggsave(pie.oth_trt, file = "05. Hypoxia group Other treatment pie.pdf", width = 14, height = 10, units = "in")

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
lgg.sub <- lgg.sub %>% dplyr::select(c(MetaProg1_AC, MetaProg2_OPC, MetaProg3_NPC, MetaProg4_Cycling, Hypoxia_cox)) 
lgg.sub <- lgg.sub %>% pivot_longer(c(MetaProg1_AC, MetaProg2_OPC, MetaProg3_NPC, MetaProg4_Cycling, Hypoxia_cox), names_to = "source", values_to = "target")
lgg.sub <- lgg.sub %>% mutate(weight = 1) %>% filter(target != " ") %>% filter(target != "") 

duplicated_genes <- lgg.sub$target[duplicated(lgg.sub$target)]
duplicated_genes

lgg.sub <- lgg.sub[!duplicated(lgg.sub[, c("source", "target")]), ]

# Run GSEA ("gsva", "plage", "ssgsea", "zscore") 

res_ssgsea <- decoupleR::run_gsva(mat = counts.normz, 
                                  network = lgg.sub, 
                                  .source ='source', 
                                  .target ='target', 
                                  minsize = 5L, 
                                  method = c("ssgsea")) 

res_long <- res_ssgsea %>% 
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% 
  t() %>% 
  scale() %>% 
  as.data.frame() 

# Confirm hypoxia and its distribution 

res_hpx_only <- res_long['Hypoxia_cox']

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
# ggsave(ht.hpx, filename = "6. Heatmap hypoxia.pdf", width = 8, height = 2, units = c("in"))

# Rename the groups 

metadata <- metadata %>% mutate(hypoxia_class = ifelse(hypoxia_groups %in% "C1", "Mild_hypoxia", "Severe_hypoxia"))

# ── Association with MPs ──────────────────────────────────────────────────────   

# Distribution 

res_ssgsea <- res_ssgsea %>% mutate(SUBTYPE_GRADE = metadata[condition, "SUBTYPE_GRADE"])

p.box <- ggplot(res_ssgsea, aes(x = source, y = score, fill = SUBTYPE_GRADE)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdBu") + 
  stat_compare_means(method = "t.test", label = "p.signif") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p.box
ggsave(p.box, filename = "7. Metaprograms boxplot.pdf", width = 4, height = 5, units = c("in"))

# Enrichment 

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

# Enrichment - average (4 groups)

enrich_df <- cbind(res_long, metadata) 

enrich_mp <- enrich_df %>%
  select(SAMPLE_ID, SUBTYPE_GRADE, hypoxia_class, MetaProg1_AC, MetaProg2_OPC, MetaProg3_NPC, MetaProg4_Cycling) %>%
  pivot_longer(cols = starts_with("MetaProg"), names_to = "Metaprogram", values_to = "Enrichment")

mps_summm <- enrich_mp %>%
  group_by(SUBTYPE_GRADE, hypoxia_class, Metaprogram) %>%
  summarise(Enrichment = mean(Enrichment, na.rm = TRUE), .groups = "drop")

mps_summm <- na.omit(mps_summm)

p.enr.dot <- ggplot(mps_summm, aes(x = hypoxia_class, y = Metaprogram, fill = Enrichment)) +
  geom_point(aes(size = abs(Enrichment)), shape = 21, color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~ SUBTYPE_GRADE, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p.enr.dot
ggsave(p.enr.dot, filename = "7. Metaprograms dotplot.pdf", width = 6, height = 3, units = c("in"))

# Correlation between hypoxia and metaprograms (4 groups) 

enrich_df <- cbind(res_long, metadata) 

cor_subdf <- enrich_df %>%
  pivot_longer(cols = starts_with("MetaProg"), names_to = "Metaprogram", values_to = "Metaprogram_value") %>%
  group_by(SUBTYPE_GRADE, hypoxia_class, Metaprogram) %>%
  summarise(Correlation = cor(Metaprogram_value, Hypoxia_cox, method = "pearson"), .groups = "drop")

cor_subdf <- na.omit(cor_subdf)

p.cor.dot <- ggplot(cor_subdf, aes(x = hypoxia_class, y = Metaprogram, fill = Correlation)) +
  geom_point(aes(size = abs(Correlation)), shape = 21, color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, limits = c(-0.67, 0.67)) +
  facet_wrap(~ SUBTYPE_GRADE, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p.cor.dot
ggsave(file = "7. Metaprograms correlation (4 categories).pdf", plot = p.cor.dot, width = 6, height = 3, units = "in") 

# ── Survival and cox analysis ───────────────────────────────────────────────── 

metasurv <- cbind(res_hpx_only, metadata)
metasurv <- metasurv[complete.cases(metasurv[, c('OS_MONTHS', 'OS_STATUS', 'hypoxia_class', 'SUBTYPE_GRADE', 'TREATMENT_new')]), ]

# Univariate survival / overall

surv.obj_hpx <- Surv(time = metasurv$OS_MONTHS, event = metasurv$OS_STATUS)
fit.surv_hpx <- survfit(surv.obj_hpx ~ hypoxia_class, data = metasurv)

kmp.surv.hpx <- ggsurvplot(fit.surv_hpx, 
                           data = metasurv, 
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

# Split by pathology 

fit.surv_pat <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ hypoxia_class, data = metasurv)
kmp.surv.pat <- ggsurvplot_facet(fit.surv_pat,
                                 data = metadata, 
                                 facet.by = "SUBTYPE", 
                                 surv.median.line = "hv", 
                                 pval = TRUE,
                                 size = 1.0,
                                 censor.size = 0) 
kmp.surv.pat 
ggsave(file = "8. Survival hypoxia groups (pathology).pdf", plot = kmp.surv.pat, width = 8, height = 4, units = "in") 

# Split by grade 

fit.surv_grd <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ hypoxia_class, data = metasurv)
kmp.surv.grd <- ggsurvplot_facet(fit.surv_grd,
                                 data = metadata, 
                                 facet.by = "GRADE", 
                                 surv.median.line = "hv", 
                                 pval = TRUE,
                                 size = 1.0,
                                 censor.size = 0) 
kmp.surv.grd 
ggsave(file = "8. Survival hypoxia groups (grade).pdf", plot = kmp.surv.grd, width = 8, height = 4, units = "in") 

# Split by hypoxia groups for pathology and grade (4 groups)

fit.surv_hpg <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ SUBTYPE_GRADE, data = metasurv)
kmp.surv.hpg <- ggsurvplot_facet(fit.surv_hpg,
                                 data = metasurv,
                                 facet.by = "hypoxia_class",
                                 surv.median.line = "hv",
                                 pval = TRUE,
                                 size = 1,
                                 censor.size = 0) 
kmp.surv.hpg 
ggsave(file = "8. Survival hypoxia groups (hypoxia in pathology & grade).pdf", plot = kmp.surv.hpg, width = 8, height = 4, units = "in") 

# Split by pathology and grade (4 groups) for hypoxia groups 

fit.surv_ptg <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ hypoxia_class, data = metasurv)
kmp.surv.ptg <- ggsurvplot_facet(fit.surv_ptg,
                                 data = metadata, 
                                 facet.by = "SUBTYPE_GRADE", 
                                 surv.median.line = "hv", 
                                 pval = TRUE, 
                                 size = 1.0,
                                 censor.size = 0, 
                                 nrow = 1) 
kmp.surv.ptg 
ggsave(file = "8. Survival hypoxia groups (pathology & grade in hypoxia).pdf", plot = kmp.surv.ptg, width = 12, height = 4, units = "in") 

# Multivariable cox regression 

surv.obj_hpx <- Surv(time = metasurv$OS_MONTHS, event = metasurv$OS_STATUS)

cox_nrm <- coxph(surv.obj_hpx ~ hypoxia_class + GRADE + SUBTYPE + TREATMENT_new, data = metasurv)
for_nrm <- ggforest(cox_nrm, data = metasurv) 
for_nrm 
ggsave(file = "8. Cox regression-norm.pdf", plot = for_int, width = 6, height = 2.5, units = "in") 

# ── ESTIMATE ────────────────────────────────────────────────────────────────── 

estim_score <- counts.normz %>% 
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = TRUE)  

write.csv(estim_score, "Result ESTIMATE scores.csv")

estim_score %>% 
  plot_purity(is_affymetrix = TRUE) 

estim_score <- estim_score %>% as.data.frame() %>% column_to_rownames(var = "sample") 
meta_estimm <- cbind(estim_score, metadata %>% select(SUBTYPE_GRADE, hypoxia_class)) 
estimm_long <- meta_estimm %>% 
  as.data.frame() %>% 
  pivot_longer(cols = !c(hypoxia_class, SUBTYPE_GRADE), names_to = "estimate_pram", values_to = "scores") %>% 
  na.omit()

# Box plot 

box_estim <- ggplot(estimm_long, aes(x = hypoxia_class, y = scores, fill = SUBTYPE_GRADE)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) + 
  facet_wrap(~estimate_pram, nrow = 1, scale = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
box_estim
ggsave(filename = "9. Estimate hypoxia groups bar (4 class).pdf", plot = box_estim, width = 8, height = 3, units = c("in")) 

# ── Immune checkpoint gene expression (ICGs) ────────────────────────────────── 

colnames(metadata)

icg_vec <- c('CCR5', 'CD274', 'CTLA4', 'CXCR3', 'HAVCR2', 'LAG3', 'PDCD1', 'PDCD1LG2', 'SIGLEC15', 'TIGIT') 

counts.icg <- counts.normz %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% icg_vec)
counts.icg <- t(counts.icg) 
counts.icg <- cbind(counts.icg, metadata %>% select(hypoxia_class, SUBTYPE_GRADE))

counts.icg.long <- counts.icg %>% 
  as.data.frame() %>% 
  pivot_longer(cols = !c(hypoxia_class, SUBTYPE_GRADE), names_to = "ICGs", values_to = "normalized_expression") %>% 
  na.omit()

# Box plot 

box_icg <- ggplot(counts.icg.long, aes(x = hypoxia_class, y = normalized_expression, fill = SUBTYPE_GRADE)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) + 
  facet_wrap(~ICGs, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
box_icg
ggsave(filename = "10. ICG hypoxia groups box (4 class).pdf", plot = box_icg, width = 10, height = 4, units = c("in")) 
