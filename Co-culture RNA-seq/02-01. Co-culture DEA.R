# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(readxl) 
library(tidyverse) 
library(dplyr) 
library(tibble)
library(tidyr) 
library(openxlsx) 
library(ggrepel) 
library(ggplot2)
library(ggplotify) 
library(ggpubr) 
library(ComplexHeatmap)
library(org.Hs.eg.db) 
library(AnnotationDbi) 
library(limma) 
library(edgeR) 
library(viridis) 
library(fgsea) 
library(msigdbr) 
library(clusterProfiler) 
library(circlize) 

# ── Data load and manipulation ──────────────────────────────────────────────── 

path_savefile <- "C:\\01. Hypoxia analysis\\1. Bulk RNA-seq\\Limma & GSEA (SC)\\1. HP_mic-NX_mic\\"
new_data <- read.csv("TMM_normalized_log2CPM.csv", row.names = 1) 
new_data <- new_data %>% dplyr::select(!c(Undetermined, N33, H33))

metadata <- readxl::read_excel("Sample_Identifier.xlsx") 
metadata <- metadata %>% column_to_rownames(var = "Code")
metadata$Group_Condition <- paste(metadata$Group, metadata$Condition, sep = "_")

all(rownames(metadata) %in% colnames(new_data)) 
all(colnames(new_data) %in% rownames(metadata)) 

# ── Limma ─────────────────────────────────────────────────────────────────────  

# Pre-processing and filtering

d0 <- DGEList(new_data)
d0 <- calcNormFactors(d0, method = "TMM")
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 

# Design matrix 

design <- model.matrix(~ 0 + Group_Condition, data = metadata)
colnames(design) <- gsub("Group_Condition", "", colnames(design))

# Fit and compute contrast 

log2cpm <- cpm(d, log = TRUE, prior.count = 1)
fit <- lmFit(log2cpm, design) 
contr <- makeContrasts(LGG_MICROGLIA_Hypoxia - LGG_MICROGLIA_Normoxia, levels = colnames(coef(fit)))
contr

# Estimate contrast for each gene 

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

# ── Volcano Plot ──────────────────────────────────────────────────────────────  

limmadata <- top.table
limmadata <- limmadata %>% na.omit()
limmadata$nlog10 <- -log10(limmadata$P.Value) 
limmadata$expression = ifelse(limmadata$P.Value < 0.05 & abs(limmadata$logFC) >= 1, 
                       ifelse(limmadata$logFC > 1 ,'Up-regulated','Down-regulated'), 'Stable') 

limmadata %>% dplyr::count(expression)
write.csv(limmadata, file = file.path(path_savefile, "Result Limma (HP_mic-NX_mic).csv"), row.names = TRUE) 

# Top Genes (up-regulated and down-regulated)

top <- 10
top_limma.genes <- bind_rows(limmadata %>% 
                               filter(expression == 'Up-regulated') %>% 
                               arrange(P.Value, desc(abs(logFC))) %>% 
                               head(top), 
                             limmadata %>% 
                               filter(expression == 'Down-regulated') %>% 
                               arrange(P.Value, desc(abs(logFC))) %>% 
                               head(top)) 
# Volcano - manual code 

vol_limm <- ggplot(data = limmadata, aes(x = logFC, y = nlog10, col = AveExpr)) + 
  geom_vline(xintercept = c(-1 , 1), lty = 2, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5) + 
  geom_point(size = 2, alpha = 0.5, stroke = 0.0) + 
  scale_color_viridis() + 
  theme_bw() + 
  ggrepel::geom_text_repel(aes(label = SYMBOL), data = top_limma.genes, size = 2, force = 10, fontface = "italic", max.overlaps = Inf) 
vol_limm 
ggsave(filename = file.path(path_savefile, "1. Volcano Plot (HP_mic-NX_mic).pdf"), plot = vol_limm, width = 5, height = 4) 

# ── Read MSigDB specific genesets ───────────────────────────────────────────── 

homo_genesets <- msigdbr(species = "Homo sapiens") 
view_genesets <- homo_genesets %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat) 

homo.hall_gsea <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol) 
homo.gobp_gsea <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol) 
homo.kegg_gsea <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(gs_name, gene_symbol) 
homo.wiki_gsea <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, gene_symbol) 
homo.rect_gsea <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% dplyr::select(gs_name, gene_symbol) 

homo.hall_gsea <- homo.hall_gsea %>% mutate(class = "Hallmark")  
homo.gobp_gsea <- homo.gobp_gsea %>% mutate(class = "GOBP")  
homo.kegg_gsea <- homo.kegg_gsea %>% mutate(class = "KEGG") 
homo.wiki_gsea <- homo.wiki_gsea %>% mutate(class = "Wiki") 
homo.rect_gsea <- homo.rect_gsea %>% mutate(class = "Reactome") 

homo.hall_genelist <- split(homo.hall_gsea$gene_symbol, homo.hall_gsea$gs_name) 
homo.gobp_genelist <- split(homo.gobp_gsea$gene_symbol, homo.gobp_gsea$gs_name) 
homo.kegg_genelist <- split(homo.kegg_gsea$gene_symbol, homo.kegg_gsea$gs_name) 
homo.wiki_genelist <- split(homo.wiki_gsea$gene_symbol, homo.wiki_gsea$gs_name) 
homo.rect_genelist <- split(homo.rect_gsea$gene_symbol, homo.rect_gsea$gs_name) 

# ── Run GSEA ──────────────────────────────────────────────────────────────────  

# Limma results 

limmadata <- read.csv(paste0(path_savefile, "Result Limma (HP_mic-NX_mic).csv")) 
limmadata <- limmadata %>% dplyr::rename(SYMBOL = "X")
limma_man <- limmadata %>% na.omit() 

# Pre-ranking  

limma_sub <- limma_man %>% dplyr::select(SYMBOL, logFC) 
limma_rnk <- limma_sub$logFC 
names(limma_rnk) <- as.character(limma_sub$SYMBOL)
limma_rnk <- sort(limma_rnk, decreasing = TRUE) 

# Run GSEA - Hallmark 

hall_GSEA_res <- GSEA(limma_rnk, TERM2GENE = homo.hall_gsea, verbose = FALSE)
hall_GSEA_dfs <- as_tibble(hall_GSEA_res@result) 
write.csv(hall_GSEA_dfs, file = file.path(path_savefile, "Result GSEA Hallmark (HP_mic-NX_mic).csv"), row.names = FALSE)  

# Run GSEA - GOBP 

gobp_GSEA_res <- GSEA(limma_rnk, TERM2GENE = homo.gobp_gsea, verbose = FALSE)
gobp_GSEA_dfs <- as_tibble(gobp_GSEA_res@result) 
write.csv(gobp_GSEA_dfs, file = file.path(path_savefile, "Result GSEA GOBP (HP_mic-NX_mic).csv"), row.names = FALSE)  

# Run GSEA - KEGG 

kegg_GSEA_res <- GSEA(limma_rnk, TERM2GENE = homo.kegg_gsea, verbose = FALSE)
kegg_GSEA_dfs <- as_tibble(kegg_GSEA_res@result) 
write.csv(kegg_GSEA_dfs, file = file.path(path_savefile, "Result GSEA KEGG (HP_mic-NX_mic).csv"), row.names = FALSE)  

# Run GSEA - Wikipathways 

wiki_GSEA_res <- GSEA(limma_rnk, TERM2GENE = homo.wiki_gsea, verbose = FALSE)
wiki_GSEA_dfs <- as_tibble(wiki_GSEA_res@result) 
write.csv(wiki_GSEA_dfs, file = file.path(path_savefile, "Result GSEA Wiki (HP_mic-NX_mic).csv"), row.names = FALSE)  

# Run GSEA - Reactome 

rect_GSEA_res <- GSEA(limma_rnk, TERM2GENE = homo.rect_gsea, verbose = FALSE)
rect_GSEA_dfs <- as_tibble(rect_GSEA_res@result) 
write.csv(rect_GSEA_dfs, file = file.path(path_savefile, "Result GSEA Reactome (HP_mic-NX_mic).csv"), row.names = FALSE)  

# All in one plot 

file_list <- list.files(path = path_savefile, pattern = "Hallmark.*\\.csv$", full.names = TRUE) 
read_and_label <- function(file) {
  sample_name <- str_trim(str_extract(file, "(?<=\\().*(?=\\))")) 
  df <- read.csv(file)
  df$condition <- sample_name 
  return(df)
}

# Read all files and add sample name column and merge 

df_listt <- lapply(file_list, read_and_label)
H_merged <- Reduce(function(x, y) full_join(x, y, by = intersect(names(x), names(y))), df_listt)
H_merged$minus_log10_pvalue <- -log10(H_merged$pvalue)

H_merged %>% count(condition)

H_astro <- H_merged %>% dplyr::filter(condition %in% c('HP_ast-NX_ast', 'HP_ast-HP_tum'))
H_micro <- H_merged %>% dplyr::filter(condition %in% c('HP_mic-NX_mic', 'HP_mic-HP_tum'))
H_basic <- H_merged %>% dplyr::filter(condition %in% c('HP_tum-NX_tum'))

# All condition (only significant pathways) 

bub_all <- ggplot(H_merged, aes(x = condition, y = ID, color = NES, size = minus_log10_pvalue)) +
  geom_point(alpha = 0.7) + 
  scale_color_viridis(option = "C") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bub_all
ggsave(filename = "2. All Hallmark pathways.pdf", plot = bub_all, width = 7, height = 12, units = c("in")) 
