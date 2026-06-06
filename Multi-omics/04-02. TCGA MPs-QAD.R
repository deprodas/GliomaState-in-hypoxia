# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(ggplot2) 
library(ggplotify)
library(patchwork) 
library(limma) 
library(edgeR) 
library(ggrepel) 
library(fgsea) 
library(msigdbr) 
library(clusterProfiler) 
library(GseaVis) 
library(tibble) 
library(ComplexHeatmap) 
library(circlize) 
library(viridis) 
library(decoupleR)
library(glmnet)

# ── Prepare data ────────────────────────────────────────────────────────────── 

# Metadata (only IDH-mutant) 

metadata <- read.csv("IDHm_metadata.csv", row.names = 1) 
metadata$SUBTYPE_GRADE <- interaction(metadata$SUBTYPE, metadata$GRADE)
metadata %>% dplyr::count(SUBTYPE_GRADE) 
metadata <- metadata[!is.na(metadata$SUBTYPE_GRADE), ] 

# Count data (only IDH-mutant) 

htseq_counts <- read.csv("TCGA-PanCan IDHm-glioma raw-count.csv", row.names = 1) 
colnames(htseq_counts) <- gsub('\\.', '-', colnames(htseq_counts)) 
htseq_counts <- htseq_counts[, rownames(metadata)] 

# Convert NA values to 0 in count matrix

htseq_counts[is.na(htseq_counts)] <- 0
all(is.numeric(htseq_counts))
all(rownames(metadata) %in% colnames(htseq_counts)) 
all(colnames(htseq_counts) %in% rownames(metadata)) 

# Checking input format 

class(metadata) 
class(htseq_counts)

# ── Run limma ───────────────────────────────────────────────────────────────── 

set.seed(123) 

groups_cat <- unique(metadata$SUBTYPE_GRADE)
cat.lima_res <- list()

for (g in groups_cat) {
  cat("running limma for:", g, "\n")
  meta_sub <- metadata %>% filter(SUBTYPE_GRADE == g)
  meta_sub <- meta_sub %>% filter(hypoxia_class %in% c("Severe_hypoxia", "Mild_hypoxia")) 
  if (length(unique(meta_sub$hypoxia_class)) < 2) {
    cat("skipping", g, "- only one hypoxia group present\n")
    next
  }
  counts_sub <- htseq_counts[, rownames(meta_sub)]
  d0_int <- DGEList(counts_sub)
  d0_int <- calcNormFactors(d0_int, method = "TMM")
  keeeep <- rowSums(cpm(d0_int) > 1) >= 2
  d_noww <- d0_int[keeeep, ]
  design <- model.matrix(~ 0 + hypoxia_class, data = meta_sub)
  colnames(design) <- gsub("hypoxia_class", "", colnames(design))
  log2cpm <- cpm(d_noww, log = TRUE, prior.count = 1)
  fit.1 <- lmFit(log2cpm, design)
  contr <- makeContrasts(Severe_hypoxia - Mild_hypoxia, levels = colnames(coef(fit.1)))
  tmp.1 <- contrasts.fit(fit.1, contr)
  tmp.1 <- eBayes(tmp.1) 
  top.table <- topTable(tmp.1, sort.by = "P", n = Inf)
  top.table$SUBTYPE_GRADE <- g 
  cat.lima_res[[g]] <- top.table
  write.csv(top.table, paste0("Result limma SH_vs_MH_", g, ".csv"))
}

# Combine all 

comb_res <- bind_rows(cat.lima_res)

# ── Plot limma ──────────────────────────────────────────────────────────────── 

limmares <- comb_res %>%
  na.omit() %>%
  rownames_to_column("Gene") %>%
  mutate(Gene_clean = sub("\\..*", "", Gene), nlog10 = -log10(P.Value)) 

# Set condition 

limmares$expression = ifelse(limmares$P.Value < 0.05 & abs(limmares$logFC) >= 0.5, 
                      ifelse(limmares$logFC > 0.5,'Up-regulated','Down-regulated'), 'Stable') 
limmares %>% group_by() %>% dplyr::count(expression)
write.csv(limmares, "Result limma SH_vs_MH_4_SUBTYPE_GRADE.csv", row.names = F)

# Top genes (up-regulated and down-regulated)

topgene_num <- 10 
topgene_lab <- limmares %>%
  group_by(SUBTYPE_GRADE) %>%
  group_modify(~ bind_rows(.x %>% filter(expression == "Up-regulated") %>% arrange(P.Value, desc(abs(logFC))) %>% head(topgene_num), 
                           .x %>% filter(expression == "Down-regulated") %>% arrange(P.Value, desc(abs(logFC))) %>% head(topgene_num))) %>% 
  ungroup() 

# Volcano plot 

p.vol_lim1 <- ggplot(limmares, aes(x = logFC, y = nlog10, fill = AveExpr)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
  geom_point(shape = 21, size = 1.2, alpha = 0.9, stroke = 0.1) +
  scale_fill_viridis_c(option = "D") +
  facet_wrap(~ SUBTYPE_GRADE, scales = "free", nrow = 1) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_bw(base_size = 12) +
  labs(x = "log2 FC", y = "-log10(P-value)") +
  geom_text_repel(data = topgene_lab, aes(label = Gene_clean), size = 2.5, fontface = "italic", max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2, segment.size = 0.2) 
p.vol_lim1
ggsave(filename = "01. Volcano limma-1.pdf", plot = p.vol_lim1, width = 10, height = 3, units = c("in")) 

# ── Read msigdb specific genes ──────────────────────────────────────────────── 

homo_genesets <- msigdbr(species = "Homo sapiens") 

see_genesets <- homo_genesets %>%
  dplyr::distinct(gs_collection, gs_subcollection) %>%
  dplyr::arrange(gs_collection, gs_subcollection) 

get_genesets <- function(species, category, subcategory = NULL, class_name) {
  msigdbr(species = species, category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, gene_symbol) %>% 
    dplyr::mutate(class = class_name)
}

hall_dbs <- get_genesets("Homo sapiens", category = "H", subcategory = NULL, class_name = "Hallmark")
gobp_dbs <- get_genesets("Homo sapiens", category = "C5", subcategory = "GO:BP", class_name = "GOBP")
kegg_dbs <- get_genesets("Homo sapiens", category = "C2", subcategory = "CP:KEGG_MEDICUS", class_name = "KEGG")
wiki_dbs <- get_genesets("Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS", class_name = "Wiki")
rect_dbs <- get_genesets("Homo sapiens", category = "C2", subcategory = "CP:REACTOME", class_name = "Reactome")
cell_dbs <- get_genesets("Homo sapiens", category = "C8", subcategory = NULL, class_name = "Celltype")

gsea_dbs <- list(Hallmark = hall_dbs, GOBP = gobp_dbs, KEGG = kegg_dbs, Wiki = wiki_dbs, Reactome = rect_dbs, Celltype = cell_dbs)

# ── Run GSEA ────────────────────────────────────────────────────────────────── 

# GSEA function  

run_gsea <- function(gene_list, term2gene, db_name, subtype_name) {
  res_ob <- GSEA(gene_list, TERM2GENE = term2gene, verbose = FALSE)
  if (is.null(res_ob)) return(NULL)
  res_df <- as_tibble(res_ob@result)
  res_df$database <- db_name
  res_df$SUBTYPE_GRADE <- subtype_name
  write.csv(res_df, paste0("Result GSEA_", subtype_name, "_", db_name, ".csv"), row.names = FALSE)
  return(res_df)
}

all.gsea_res <- list()

for (subtype in unique(comb_res$SUBTYPE_GRADE)) {
  cat("Running GSEA for:", subtype, "\n")
  limma_subbb <- comb_res %>% filter(SUBTYPE_GRADE == subtype) %>% rownames_to_column("SYMBOL") %>% dplyr::select(SYMBOL, logFC) %>% na.omit()
  gene_ranked <- limma_subbb$logFC
  names(gene_ranked) <- limma_subbb$SYMBOL
  gene_ranked <- sort(gene_ranked, decreasing = TRUE)
  for (db in names(gsea_dbs)) {
    gsea_res <- run_gsea(gene_ranked, gsea_dbs[[db]], db, subtype)
    if (!is.null(gsea_res)) {
      all.gsea_res[[paste(subtype, db, sep = "_")]] <- gsea_res
    }
  }
}

# Combine GSEA

gsea_all <- bind_rows(all.gsea_res)
write.csv(gsea_all, "Result GSEA SH_vs_MH_4_SUBTYPE_GRADE.csv", row.names = F)

# ── Plot GSEA ───────────────────────────────────────────────────────────────── 

# Bubble plot 

bubb_paths <- gsea_all %>%
  filter(p.adjust < 0.05) %>% 
  mutate(neglog10p = -log10(pvalue)) %>% 
  group_by(SUBTYPE_GRADE, database) %>%
  slice_max(order_by = abs(NES), n = 10) %>%
  ungroup() 

bubb_plots <- list()

for (db in unique(bubb_paths$database)) {
  df_sub <- bubb_paths %>% filter(database == db)
  p_bubb <- ggplot(df_sub, aes(x = Description, y = SUBTYPE_GRADE, size = neglog10p, color = enrichmentScore)) +
    geom_point() +
    scale_color_viridis(option = "D") +
    scale_size(range = c(2, 8)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  bubb_plots[[db]] <- p_bubb
}

p.bubb_all <- wrap_plots(bubb_plots, ncol = 3)
p.bubb_all
ggsave(filename = "02. GSEA all bubble.pdf", plot = p.bubb_all, width = 40, height = 25, units = c("in"), limitsize = FALSE) 

# Heatmap 

heat_paths <- gsea_all %>%
  filter(database %in% c('Hallmark', 'Wiki', 'Reactome')) %>%
  filter(pvalue < 0.05) %>%
  group_by(SUBTYPE_GRADE) %>%
  slice_max(order_by = abs(NES), n = 30) %>%
  ungroup() 

heat_paths <- gsea_all %>% filter(Description %in% gene_sets)

heat_data <- heat_paths %>%
  select(database, Description, SUBTYPE_GRADE, NES) %>%
  pivot_wider(id_cols = c(database, Description), names_from = SUBTYPE_GRADE, values_from = NES, values_fill = 0)

heat_matx <- heat_data %>%
  select(-database) %>%
  column_to_rownames('Description') %>%
  t() %>% scale() %>% t() %>% 
  as.matrix() 

annotation_row <- data.frame(Database = heat_data$database)
rownames(annotation_row) <- heat_data$Description

set.seed(123)
ht_fil <- ComplexHeatmap::pheatmap(heat_matx, 
                                   color = colorRampPalette(c("darkblue", "white", "red"))(100), 
                                   row_split = annotation_row$Database, 
                                   annotation_row = annotation_row)
ht_fil
ht_fil <- as.ggplot(ht_fil)
ggsave(filename = "02. GSEA heatmap.pdf", plot = ht_fil, width = 5.75, height = 8.5, units = c("in"), limitsize = FALSE) 

# ── Limma (SH-MH) ─────────────────────────────────────────────────────────────

# Limma 

set.seed(123) 

groups_hpx <- unique(metadata$CANCER_TYPE)
hpx.limma_res <- list()

for (g in groups_hpx) {
  cat("running limma for:", g, "\n")
  meta_sub <- metadata %>% filter(hypoxia_class %in% c("Severe_hypoxia", "Mild_hypoxia")) 
  if (length(unique(meta_sub$hypoxia_class)) < 2) {
    cat("skipping", g, "- only one hypoxia group present\n")
    next
  }
  counts_sub <- htseq_counts[, rownames(meta_sub)]
  d0_int <- DGEList(counts_sub)
  d0_int <- calcNormFactors(d0_int, method = "TMM")
  keeeep <- rowSums(cpm(d0_int) > 1) >= 2
  d_noww <- d0_int[keeeep, ]
  design <- model.matrix(~ 0 + hypoxia_class, data = meta_sub)
  colnames(design) <- gsub("hypoxia_class", "", colnames(design))
  log2cpm <- cpm(d_noww, log = TRUE, prior.count = 1)
  fit.1 <- lmFit(log2cpm, design)
  contr <- makeContrasts(Severe_hypoxia - Mild_hypoxia, levels = colnames(coef(fit.1)))
  tmp.1 <- contrasts.fit(fit.1, contr)
  tmp.1 <- eBayes(tmp.1) 
  top.table <- topTable(tmp.1, sort.by = "P", n = Inf)
  top.table$SUBTYPE_GRADE <- g 
  hpx.limma_res[[g]] <- top.table
  write.csv(top.table, paste0("Result limma SH_vs_MH_", g, ".csv"))
}

hpx.limma_data <- bind_rows(hpx.limma_res)

# ── GSEA using fgsea (SH-MH) ──────────────────────────────────────────────────

midh_genes <- read.xlsx("mIDH programs_depro.xlsx")
midh_fixed <- lapply(names(midh_genes), function(x) {
  na.omit(midh_genes[[x]])
})
names(midh_fixed) <- names(midh_genes)

# Prepare ranked gene list 

hpx.gene_ranks <- hpx.limma_data$logFC
names(hpx.gene_ranks) <- rownames(hpx.limma_data)
hpx.gene_ranks <- tapply(hpx.gene_ranks, names(hpx.gene_ranks), mean)
hpx.gene_ranks <- sort(hpx.gene_ranks, decreasing = TRUE)

# Run fgsea 

fgsea_hpx.ress <- fgsea(pathways = midh_fixed, stats = hpx.gene_ranks, minSize = 5, maxSize = 500, nperm = 10000) 
fgsea_hpx.ress <- fgsea_hpx.ress[order(fgsea_hpx.ress$NES, decreasing = TRUE), ]

# Plot 

n_plot <- 20
fgsea_hpx.path <- fgsea_hpx.ress$pathway
rank_positions <- data.frame(gene = names(hpx.gene_ranks), position = 1:length(hpx.gene_ranks))

for(pw in fgsea_hpx.path){
  lead_genes <- fgsea_hpx.ress %>%
    filter(pathway == pw) %>%
    pull(leadingEdge) %>%
    .[[1]] 
  best_genes <- data.frame(gene = lead_genes) %>%
    filter(gene %in% names(hpx.gene_ranks)) %>%
    mutate(rank_weighted = abs(hpx.gene_ranks[gene])) %>%
    arrange(desc(rank_weighted)) %>%
    slice(1:min(n_plot, n())) %>%
    left_join(rank_positions, by = "gene") 
  
  nes <- round(fgsea_hpx.ress$NES[fgsea_hpx.ress$pathway == pw], 2)
  pvl <- signif(fgsea_hpx.ress$pval[fgsea_hpx.ress$pathway == pw], 2)
  fdr <- signif(fgsea_hpx.ress$padj[fgsea_hpx.ress$pathway == pw], 2) 
  
  p_gseaplot <- plotEnrichment(midh_fixed[[pw]], hpx.gene_ranks) +
    theme_classic() +
    labs(title = paste0(pw), subtitle = paste0("NES = ", nes, " | p = ", pvl, " | FDR = ", fdr)) +
    geom_vline(data = best_genes, aes(xintercept = position), linetype = "dashed", color = "red", alpha = 0.5) +
    geom_text(data = best_genes, aes(x = position, y = 0, label = gene), angle = 90, vjust = -0.3, size = 3, color = "red")
  ggsave(filename = paste0("03. GSEA_", pw, ".pdf"), plot = p_gseaplot, width = 6, height = 4)
}

# Save leading-edge genes

n_save <- 50
hpx.lead_genes <- list()

for(pw in fgsea_hpx.path){
  lead_genes <- fgsea_hpx.ress %>%
    filter(pathway == pw) %>%
    pull(leadingEdge) %>%
    .[[1]]
  best_genes <- data.frame(gene = lead_genes) %>%
    filter(gene %in% names(hpx.gene_ranks)) %>%
    mutate(rank_weighted = abs(hpx.gene_ranks[gene])) %>%
    arrange(desc(rank_weighted)) %>%
    slice(1:min(n_save, n())) %>%
    pull(gene)
  
  if(length(best_genes) < n_save){
    best_genes <- c(best_genes, rep(NA, n_save - length(best_genes)))
  }
  hpx.lead_genes[[pw]] <- best_genes
}

hpx.leadedge_df <- as.data.frame(hpx.lead_genes, stringsAsFactors = FALSE)
write.xlsx(hpx.leadedge_df, "Result Leading-edge genes (MPs).xlsx")

# ── MP, QAD validation ssgsea ─────────────────────────────────────────────────

colnames(midh_genes) 
midh_sgsea <- midh_genes %>% pivot_longer(c(everything()), names_to = "source", values_to = "target")
midh_sgsea <- midh_sgsea %>% mutate(weight = 1) %>% filter(target != " ") %>% filter(target != "") 

duplicated_genes <- midh_sgsea$target[duplicated(midh_sgsea$target)]
duplicated_genes

midh_sgsea <- midh_sgsea[!duplicated(midh_sgsea[, c("source", "target")]), ]

# Run ssGSEA ("gsva", "plage", "ssgsea", "zscore") 

res_ssgsea <- decoupleR::run_gsva(mat = htseq_counts, 
                                  network = midh_sgsea, 
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

enrich_df <- cbind(res_long, metadata) 

enrich_mp <- enrich_df %>%
  select(SAMPLE_ID, SUBTYPE_GRADE, hypoxia_class, MetaProg1_AC, MetaProg2_OPC, MetaProg3_NPC, MetaProg4_Cycling, Q_GBM, A_GBM, D_GBM, EGFR_50, EGFR_all) %>%
  pivot_longer(cols = c(MetaProg1_AC, MetaProg2_OPC, MetaProg3_NPC, MetaProg4_Cycling, Q_GBM, A_GBM, D_GBM, EGFR_50, EGFR_all), names_to = "Metaprogram", values_to = "Enrichment")

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
ggsave(p.enr.dot, filename = "04. MP-QAD dotplot.pdf", width = 5, height = 4, units = c("in"))

