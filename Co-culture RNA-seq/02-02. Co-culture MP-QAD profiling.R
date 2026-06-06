# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(readxl) 
library(tidyverse) 
library(ggplot2) 
library(ggpubr) 
library(ComplexHeatmap)
library(ggplotify) 
library(decoupleR)
library(dplyr) 
library(tibble) 
library(ggrepel) 
library(ComplexHeatmap) 
library(circlize) 

# ── Data load and manipulation ──────────────────────────────────────────────── 

countdata <- read.delim('TMM_normalized_log2CPM.txt') 
countdata <- countdata %>% dplyr::select(!c(N33, H33))

metaadata <- readxl::read_excel('Sample_Identifier.xlsx') 
metaadata <- metaadata %>% dplyr::filter(!Code %in% c('N33', 'H33'))

# ── Basic gene expression ───────────────────────────────────────────────────── 

gene_data <- countdata %>% 
  dplyr::filter(rownames(.) %in% c('EEF2', 'EGFR', 'EEF1A1')) %>% 
  t() %>% 
  as.data.frame() 

# Normalization and data manipulation 

gene_data <- log2(gene_data + 1)  
meta_subs <- metaadata %>% column_to_rownames(var = "Code")
gene_data <- cbind(gene_data, meta_subs)
gene_data <- gene_data %>% rownames_to_column(var = "Samples")
colnames(gene_data)
gene_data <- gene_data %>%
  pivot_longer(!c(Samples, Condition, Group, Culture_type), names_to = "selected_genes", values_to = "expression")

# Box plot

pt_g.expr <- ggplot(gene_data, aes(x = Group, y = expression, fill = Condition)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#F72C5B", "#8FD14F")) + 
  theme_bw() + 
  stat_compare_means(method = "t.test") + 
  facet_wrap(~selected_genes) 
pt_g.expr 
ggsave(filename = "1. Gene expression (boxplot-1).pdf", plot = pt_g.expr, width = 15, height = 5, units = c("in"))

# ── DocupleR analysis ─────────────────────────────────────────────────────────

# Prepare data

any(is.na(countdata)) 
counts <- countdata %>%
  as.data.frame() %>% 
  dplyr::mutate_if(~ any(is.na(.x)), ~ dplyr::if_else(is.na(.x), 0, .x)) %>% 
  as.matrix() 

design <- metaadata 

# Run enrichment 

net_prog <- decoupleR::get_progeny(organism = 'human', top = 500)
net_gsea <- net_prog %>% dplyr::select(!c(p_value))
res_prog <- decoupleR::run_gsva(mat = counts, 
                                network = net_gsea, 
                                .source ='source', 
                                .target ='target', 
                                minsize = 2L, 
                                method = c("zscore")) 

# Data manipulation

res_matx <- res_prog %>% 
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

res_matx <- scale(res_matx)

prog_tot <- cbind(res_matx, metaadata) 
prog_tot <- prog_tot %>%
  dplyr::select(!c(Code, Culture_type)) %>%
  dplyr::group_by(Group, Condition) %>%
  dplyr::summarise(across(everything(), median, na.rm = TRUE)) 

prog_tot$Group_Condition <- paste(prog_tot$Group, prog_tot$Condition, sep = "_")
prog_tot <- prog_tot %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Group_Condition") 

# Heatmap 

colnames(prog_tot)
prog_tot <- prog_tot %>% dplyr::select(!c(Group, Condition))
prog_dfs <- prog_tot %>% t() %>% scale()

set.seed(123) 
ht_progn <- Heatmap(prog_dfs, 
                    col = colorRamp2(c(-1.7, 0, 1.7), c("Darkblue", "white","red"))) 
ht_progn
ht_progn <- as.ggplot(ht_progn)
ggsave(filename = "2. Progeny pathways (group-wise).pdf", plot = ht_progn, width = 4, height = 7.5, units = c("in"))

# ── DocupleR custom ───────────────────────────────────────────────────────────

# Custom gene set 

genesets <- read_excel("mIDH programs_depro.xlsx")
colnames(genesets)
genesets <- genesets %>% pivot_longer(c(MetaProg1_AC, MetaProg2_OPC, MetaProg3_NPC, MetaProg4_Cycling), names_to = "source", values_to = "target")
genesets <- genesets %>% mutate(weight = 1) %>% filter(target!=" ") %>% filter(target!="")

# Run enrichment 

res_cust <- decoupleR::run_gsva(mat = counts, 
                                network = genesets, 
                                .source ='source', 
                                .target ='target', 
                                minsize = 2L, 
                                method = c("gsva")) 

# Data manipulation 

res_wide <- res_cust %>% 
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% 
  t() 

cust_tot <- cbind(res_wide, metaadata) 
cust_tot <- cust_tot %>%
  dplyr::select(!c(Code, Culture_type)) %>% 
  dplyr::group_by(Group, Condition) %>%
  dplyr::summarise(across(everything(), median, na.rm = TRUE)) 

cust_tot$Group_Condition <- paste(cust_tot$Group, cust_tot$Condition, sep = "_")
cust_tot <- cust_tot %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Group_Condition") 

colnames(cust_tot)
cust_tot <- cust_tot %>% dplyr::select(!c(Group, Condition)) 

# Heatmap 

cust_dfs <- cust_tot %>% scale() 

set.seed(123) 
ht_custm <- Heatmap(cust_dfs, 
                    col = colorRamp2(c(-1.5, 0, 1.5), c("Darkblue", "white","red"))) 
ht_custm
ht_custm <- as.ggplot(ht_custm) 
ggsave(filename = "2. MPs (group-wise).pdf", plot = ht_custm, width = 5, height = 3.75, units = c("in"))
