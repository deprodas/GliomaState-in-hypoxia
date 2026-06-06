# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(ggplot2) 
library(ggplotify) 
library(ggpubr) 
library(patchwork) 
library(clusterProfiler) 
library(decoupleR) 
library(ComplexHeatmap) 
library(circlize) 

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
htseq_counts[is.na(htseq_counts)] <- 0

# ── Run ssgsea ────────────────────────────────────────────────────────────────  

# Prepare geneset

miller_genesets <- read.xlsx("Geneset Miller_Nature_2023.xlsx")
microg_genesets <- miller_genesets %>%
  slice(-1) %>%
  select(matches("microglia", ignore.case = TRUE)) 

microg_sgsea <- microg_genesets %>% pivot_longer(c(everything()), names_to = "source", values_to = "target")
microg_sgsea <- microg_sgsea %>% mutate(weight = 1) %>% filter(target != " ") %>% filter(target != "") 

duplicated_genes <- microg_sgsea$target[duplicated(microg_sgsea$target)]
duplicated_genes

microg_sgsea <- microg_sgsea[!duplicated(microg_sgsea[, c("source", "target")]), ]

# Run ssgsea 

res_mic.gsea <- decoupleR::run_gsva(mat = htseq_counts, 
                                    network = microg_sgsea, 
                                    .source ='source', 
                                    .target ='target', 
                                    minsize = 5L, 
                                    method = c("ssgsea")) 
res_mic.long <- res_mic.gsea %>% 
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% 
  t() %>% 
  scale() %>% 
  as.data.frame() 

write.csv(res_mic.long, file = "Result ssGSEA (microglia miller_nature).csv")

# ── Correlation ───────────────────────────────────────────────────────────────  

# Combine previous results 

meta_micro <- metadata %>% 
  mutate(sub_grade_hpx = paste(SUBTYPE_GRADE, hypoxia_class, sep = "_")) %>% 
  select(sub_grade_hpx)

res_mp.qad <- read.csv("Result ssGSEA (MPs-QAD-Progeny).csv", row.names = 1)
res_mp.qad <- res_mp.qad %>% 
  select(A_GBM, EGFR_50) %>% 
  dplyr::filter(rownames(.) %in% rownames(res_mic.long))

comb_micro <- cbind(res_mic.long, res_mp.qad, meta_micro) 

# Correlation 

filt_micro <- comb_micro %>% select(!EGFR_50)
long_micro <- filt_micro %>% pivot_longer(cols = -c(A_GBM, sub_grade_hpx), names_to = "pathway", values_to = "score") 
cor_mic_df <- long_micro %>%
  group_by(sub_grade_hpx, pathway) %>%
  summarise(cor_val = cor(A_GBM, score, method = "pearson", use = "complete.obs"), .groups = "drop")

p_mic_qads <- ggplot(cor_mic_df, aes(x = reorder(pathway, cor_val), y = cor_val)) +
  geom_point(size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~sub_grade_hpx, nrow = 1) +
  theme_bw() 
p_mic_qads 
ggsave(filename = "01. Correlation A-stage (microglia miller_nature).pdf", plot = p_mic_qads, width = 14, height = 10, units = c("in")) 

# ── Co-localization correlation ─────────────────────────────────────────────── 

# Microglia scores

microglia_col <- cor_mic_df %>%
  dplyr::filter(cor_val > 0, sub_grade_hpx %in% c("LGG_IDHmut-codel.G3_Severe_hypoxia", "LGG_IDHmut-non-codel.G3_Severe_hypoxia")) %>%
  dplyr::group_by(pathway) %>%
  dplyr::filter(n_distinct(sub_grade_hpx) == 2) %>% 
  dplyr::ungroup() %>%
  dplyr::pull(pathway) %>%
  unique()

poss.corr_mic <- res_mic.long %>%
  mutate(microglia_avg = rowMeans(select(., all_of(microglia_col)), na.rm = TRUE)) %>% 
  select(microglia_avg)

# Add metadata

res_mp.new_df <- res_mp.qad %>% 
  select(EGFR_50) %>% 
  dplyr::filter(rownames(.) %in% rownames(res_mic.long))

meta_pathways <- metadata %>% 
  mutate(sub_grade_hpx = paste(SUBTYPE_GRADE, hypoxia_class, sep = "_")) %>% 
  select(sub_grade_hpx)

plot_pathways <- cbind(poss.corr_mic, res_mp.new_df, meta_pathways)

# Plot correlation 

p_mic_progeny <- ggplot(plot_pathways, aes(x = EGFR_50, y = microglia_avg)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~sub_grade_hpx, scales = "free", nrow = 1) +
  stat_cor(data = subset(plot_pathways), method = "pearson", label.x.npc = "left", label.y.npc = 0.98, size = 3.5) +
  theme_bw() 
p_mic_progeny
ggsave("02. Correlation EGFR signaling (microglia miller_nature).pdf", p_mic_progeny, width = 20, height = 2.75)

# ── Specific marker genes ─────────────────────────────────────────────────────  

# Miller et al. (2025) and Chhatbar et al. (2026)

mic_marker <- c("ISG15", "IRF7", "IFITM3", "MX1", "CXCL10", "CXCL12", "CCL4", "CXCR4", "EGR1", "JUN")

infl_count <- htseq_counts %>%
  filter(rownames(.) %in% mic_marker) %>%
  mutate(across(everything(), ~log2(.x + 1)))

infl_longg <- infl_count %>%
  rownames_to_column(var = "genes") %>% 
  pivot_longer(cols = -genes, names_to = "samples", values_to = "log_expr")

meta_micro %>% dplyr::count(sub_grade_hpx)
meta_micro <- meta_micro %>% rownames_to_column(var = "samples")

infl_joind <- infl_longg %>% inner_join(meta_micro, by = "samples")
infl_avger <- infl_joind %>%
  group_by(genes, sub_grade_hpx) %>%
  summarise(mean_expr = mean(log_expr, na.rm = TRUE), .groups = "drop")

infl_final <- infl_avger %>%
  pivot_wider(names_from = sub_grade_hpx, values_from = mean_expr) %>%
  column_to_rownames("genes") %>%
  as.matrix()

infl_scled <- t(scale(t(infl_final)))
infl_scled[is.na(infl_scled)] <- 0

set.seed(123) 
ht_inf_mic <- ComplexHeatmap::Heatmap(infl_scled, 
                                      col = colorRampPalette(c("darkblue", "white", "red"))(100),
                                      cluster_rows = TRUE,
                                      cluster_columns = FALSE,
                                      show_row_names = TRUE,
                                      show_column_names = TRUE, 
                                      border = T)
ht_inf_mic
ht_inf_mic <- as.ggplot(ht_inf_mic)
ggsave(filename = "03. Markers (microglia miller_chhatbar).pdf", plot = ht_inf_mic, width = 5, height = 7, units = c("in")) 
