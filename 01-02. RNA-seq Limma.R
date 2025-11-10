# Author : Depro Das, Department of Neurosurgery, University Hospital Freiburg, Freibung, Germany 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(ggplot2) 
library(ggplotify)
library(patchwork) 
library(limma) 
library(edgeR) 
library(fgsea) 
library(msigdbr) 
library(clusterProfiler) 
library(tibble) 
library(ComplexHeatmap) 
library(circlize) 
library(viridis) 

# ── Prepare data ────────────────────────────────────────────────────────────── 

# Count data (only IDH-mutant) 

htseq_counts <- read.csv("TCGA-PanCan IDHm-glioma raw-count.csv", row.names = 1) 
colnames(htseq_counts) <- gsub('\\.', '-', colnames(htseq_counts)) 

# Convert NA values to 0 in count matrix

htseq_counts[is.na(htseq_counts)] <- 0
all(is.numeric(htseq_counts))

# Metadata (only IDH-mutant) 

metadata <- read.csv("IDHm_metadata.csv", row.names = 1) 
colnames(metadata)
metadata %>% dplyr::count(hypoxia_class) 

all(rownames(metadata) %in% colnames(htseq_counts)) 
all(colnames(htseq_counts) %in% rownames(metadata)) 


# ── Run limma ───────────────────────────────────────────────────────────────── 

set.seed(123) 
head(htseq_counts)

# Create DGEList object 

d0 <- DGEList(htseq_counts)

# Pre-processing (calculate normalization factors) 

d0 <- calcNormFactors(d0, method = "TMM")

# Filter out low-expressed genes 

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 

# Create a design matrix 

design <- model.matrix(~ 0 + hypoxia_class, data = metadata)
colnames(design) <- gsub("hypoxia_class", "", colnames(design))

# Log normalize the count matrix 

log2cpm <- cpm(d0, log = TRUE, prior.count = 1)

# Fit the expression matrix to a linear model

fit <- lmFit(log2cpm, design) 

# Compute contrast 

contr <- makeContrasts(Severe_hypoxia - Mild_hypoxia, levels = colnames(coef(fit)))
contr 

# Estimate contrast for each gene 

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

write.csv(top.table, "Results limma (SH-MH).csv") 


# ── Plot limma ──────────────────────────────────────────────────────────────── 

limmadata <- top.table %>% na.omit()
limmadata$nlog10 <- -log10(limmadata$P.Value) 

# Set condition 

limmadata$expression = ifelse(limmadata$P.Value < 0.05 & abs(limmadata$logFC) >= 0.5, 
                       ifelse(limmadata$logFC > 0.5,'Up-regulated','Down-regulated'), 'Stable') 
limmadata %>% dplyr::count(expression)


# Top genes (up-regulated and down-regulated)

top <- 10 
top_limma.genes <- bind_rows(limmadata %>% 
                               filter(expression == 'Up-regulated') %>% 
                               arrange(P.Value, desc(abs(logFC))) %>% 
                               head(top), 
                             limmadata %>% 
                               filter(expression == 'Down-regulated') %>% 
                               arrange(P.Value, desc(abs(logFC))) %>% 
                               head(top)) 
# Volcano plot 

p.vol_lim <- ggplot(data = limmadata, aes(x = logFC, y = nlog10, fill = expression)) + 
  geom_vline(xintercept = c(-0.5 , 0.5), lty = 2, col = "black", lwd = 0.5) +
  geom_hline(yintercept = 0.05, lty = 2, col = "black", lwd = 0.5) + 
  geom_point(size = 1, alpha = 1, pch = 21, stroke = 0.0001) + 
  scale_fill_manual(values = c("#FD841F", "#D6D5A8", "#38E54D")) + 
  xlim(-3, 3) + 
  theme_bw() + 
  xlab("logFC") + 
  ylab("-log10(P-value)") + 
  ggrepel::geom_text_repel(aes(label = rownames(top_limma.genes)), data = top_limma.genes, size = 2, force = 10, fontface = "italic", max.overlaps = Inf) 
p.vol_lim 
ggsave(filename = "01. Volcano limma-1.pdf", plot = p.vol_lim, width = 5, height = 4, units = c("in")) 


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

hall_gsea <- get_genesets("Homo sapiens", category = "H", subcategory = NULL, class_name = "Hallmark")
gobp_gsea <- get_genesets("Homo sapiens", category = "C5", subcategory = "GO:BP", class_name = "GOBP")
kegg_gsea <- get_genesets("Homo sapiens", category = "C2", subcategory = "CP:KEGG_MEDICUS", class_name = "KEGG")
wiki_gsea <- get_genesets("Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS", class_name = "Wiki")
rect_gsea <- get_genesets("Homo sapiens", category = "C2", subcategory = "CP:REACTOME", class_name = "Reactome")


# ── Run GSEA ────────────────────────────────────────────────────────────────── 

# Create pre-ranked gene from our limma result  

str(limmadata)

limma_fil.df <- limmadata %>% rownames_to_column(var = "SYMBOL") %>% dplyr::select(SYMBOL, logFC) 
limma_ranked <- limma_fil.df$logFC 
names(limma_ranked) <- as.character(limma_fil.df$SYMBOL)
limma_ranked <- sort(limma_ranked, decreasing = TRUE)  

# Run enrichment 

run_gsea <- function(gene_list, term2gene, db_name) {
  res_ob <- GSEA(gene_list, TERM2GENE = term2gene, verbose = FALSE)
  res_df <- as_tibble(res_ob@result)
  res_df$database <- db_name
  write.csv(res_df, paste0("Result GSEA (", db_name, ").csv"), row.names = FALSE)
  return(res_df)
}

gsea_results <- list(Hallmark = hall_gsea, GOBP = gobp_gsea, KEGG = kegg_gsea, Wiki = wiki_gsea, Reactome = rect_gsea)

gsea_dfs <- lapply(names(gsea_results), function(db) {
  run_gsea(limma_ranked, gsea_results[[db]], db)
})

gsea_all <- bind_rows(gsea_dfs)
write.csv(gsea_all, "Results GSEA-all (SH-MH).csv") 


# ── Plot GSEA ───────────────────────────────────────────────────────────────── 

# Prepare data and select pathways 

gsea_all %>% count(database) 
gsea_fil <- gsea_all %>% filter(database %in% c("Hallmark", "Reactome", "KEGG"))
gsea_fil$nlog10 <- -log10(gsea_fil$pvalue) 

path_top <- gsea_fil %>% group_by(database) %>% arrange(database, desc(NES)) %>% slice_head(n = 7) %>% ungroup() 
path_bot <- gsea_fil %>% group_by(database) %>% arrange(database, NES) %>% slice_head(n = 7) %>% ungroup() 

path_top.long <- path_top %>% 
  group_by(ID, database) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = database, values_from = NES, values_fill = NA) %>%
  column_to_rownames("ID") %>%
  as.data.frame() 

path_bot.long <- path_bot %>% 
  group_by(ID, database) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = database, values_from = NES, values_fill = NA) %>%
  column_to_rownames("ID") %>%
  as.data.frame() 

# Common features  

db_names <- unique(gsea_fil$database) 
db_color <- setNames(RColorBrewer::brewer.pal(length(db_names), "Set1"), db_names)

# Up-regulated 

ann_meta.up <- path_top %>% column_to_rownames(var = "ID") %>% select(nlog10, database) 
row_ann1.up <- rowAnnotation(nlog10 = anno_barplot(ann_meta.up$nlog10, gp = gpar(fill = "yellow2"), bar_width = 0.8, border = TRUE)) 
row_ann2.up <- rowAnnotation(Database = ann_meta.up$database, col = list(Database = db_color), border = TRUE) 

ht.up.sh <- Heatmap(as.matrix(path_top.long), 
                    col = viridisLite::plasma(10),
                    right_annotation = c(row_ann2.up, row_ann1.up),
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    na_col = "white", 
                    border = TRUE, 
                    column_title_gp = gpar(fontsize = 10, fontface = "bold")) 
ht.up.sh <- as.ggplot(ht.up.sh)
ht.up.sh 
ggsave(filename = "02. GSEA up.pdf", plot = ht.up.sh, width = 5.25, height = 6, units = c("in")) 

# Down-regulated 

ann_meta.dn <- path_bot %>% column_to_rownames(var = "ID") %>% select(nlog10, database) 
row_ann1.dn <- rowAnnotation(nlog10 = anno_barplot(ann_meta.dn$nlog10, gp = gpar(fill = "yellow2"), bar_width = 0.8, border = TRUE)) 
row_ann2.dn <- rowAnnotation(Database = ann_meta.dn$database, col = list(Database = db_color), border = TRUE) 

ht.dn.sh <- Heatmap(as.matrix(path_bot.long), 
                    col = viridisLite::plasma(10),
                    right_annotation = c(row_ann2.dn, row_ann1.dn),
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    na_col = "white", 
                    border = TRUE, 
                    column_title_gp = gpar(fontsize = 10, fontface = "bold")) 
ht.dn.sh <- as.ggplot(ht.dn.sh)
ht.dn.sh
ggsave(filename = "02. GSEA down.pdf", plot = ht.dn.sh, width = 5.25, height = 6, units = c("in")) 
