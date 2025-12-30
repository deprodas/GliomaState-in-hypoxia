# Author : Depro Das, Department of Neurosurgery, University Hospital Freiburg, Freiburg, Germany 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(ggplot2)
library(ggrepel)
library(tidyverse)
library(openxlsx) 
library(ggpubr) 

# ── Data input and prepare ────────────────────────────────────────────────────

# Metadata (only IDH-mutant) 

meta_df <- read.csv("IDHm_metadata.csv", row.names = 1)
mh_samp <- meta_df %>% filter(hypoxia_class == "Mild_hypoxia") %>% rownames() %>% as.character()
sh_samp <- meta_df %>% filter(hypoxia_class == "Severe_hypoxia") %>% rownames() %>% as.character()

# mRNA data (only IDH-mutant) 

mrna_df <- read.csv("TCGA-PanCan IDHm-glioma raw-count.csv", row.names = 1) 
colnames(mrna_df) <- gsub('\\.', '-', colnames(mrna_df)) 

mrna_log <- mrna_df %>%
  mutate(across(everything(), ~log2(as.numeric(.) + 1)))

# Protein data (only IDH-mutant) 

prot_df <- read.xlsx("data_rppa_zscores.xlsx") 

prot_df <- prot_df %>% mutate(across(-Composite.Element.REF, ~ as.numeric(gsub(" ", "", .))))
colnames(prot_df)[-1] <- gsub("-\\d{2}$", "", colnames(prot_df)[-1])

# Keep mean expression for duplicated protein names 

prot_df <- prot_df %>%
  group_by(Composite.Element.REF) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>% 
  column_to_rownames(var = "Composite.Element.REF")

# Subset IDH-mutants 

prot_df <- prot_df[, intersect(rownames(meta_df), colnames(prot_df))]

# ── Average across genes ────────────────────────────────────────────────────

mrna_mild <- mrna_log[, colnames(mrna_log) %in% mh_samp] %>%
  mutate(gene_names = rownames(.)) %>%
  rowwise() %>%
  mutate(log2_mean_mRNA = mean(c_across(-gene_names), na.rm = TRUE)) %>%
  ungroup() %>%
  select(gene_names, log2_mean_mRNA) %>%
  mutate(hypoxia_class = "Mild_hypoxia")

mrna_sevr <- mrna_log[, colnames(mrna_log) %in% sh_samp] %>%
  mutate(gene_names = rownames(.)) %>%
  rowwise() %>%
  mutate(log2_mean_mRNA = mean(c_across(-gene_names), na.rm = TRUE)) %>%
  ungroup() %>%
  select(gene_names, log2_mean_mRNA) %>%
  mutate(hypoxia_class = "Severe_hypoxia")

# Determine the minimum protein value to shift (avoid negatives)

min_protn <- min(prot_df, na.rm = TRUE)

prot_mild <- prot_df[, colnames(prot_df) %in% mh_samp] %>%
  mutate(prot_names = rownames(.)) %>%
  rowwise() %>%
  mutate(log2_mean_protein = mean(log2(c_across(-prot_names) - min_protn + 1), na.rm = TRUE)) %>%
  ungroup() %>%
  select(prot_names, log2_mean_protein) %>%
  mutate(hypoxia_class = "Mild_hypoxia")

prot_sevr <- prot_df[, colnames(prot_df) %in% sh_samp] %>%
  mutate(prot_names = rownames(.)) %>%
  rowwise() %>%
  mutate(log2_mean_protein = mean(log2(c_across(-prot_names) - min_protn + 1), na.rm = TRUE)) %>%
  ungroup() %>%
  select(prot_names, log2_mean_protein) %>%
  mutate(hypoxia_class = "Severe_hypoxia")

# Merge mRNA and protein for each group

merge_mild <- inner_join(mrna_mild, prot_mild, by = c("gene_names" = "prot_names"))
merge_sevr <- inner_join(mrna_sevr, prot_sevr, by = c("gene_names" = "prot_names"))

merge_df.p <- bind_rows(merge_mild, merge_sevr) %>%
  mutate(highlight = ifelse(gene_names %in% c("EEF2", "HIF3A", "KIF14", "TNR"), "Yes", "No"))

# Plot scatter plot 

p1.dot <- ggplot(merge_df.p, aes(x = log2_mean_mRNA, y = log2_mean_protein)) +
  geom_point(aes(color = highlight), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  geom_text_repel(data = subset(merge_df.p, highlight == "Yes"), aes(label = gene_names), size = 3) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x.npc = "left", label.y.npc = 0.9, size = 5) +
  scale_color_manual(values = c("Yes" = "red", "No" = "grey70")) +
  facet_wrap(~hypoxia_class.x) +
  labs(x = "log2 mean mRNA", y = "log2 mean protein") +
  theme_bw()
p1.dot 
ggsave(filename = "01. Protein scatter.pdf", plot = p1.dot, width = 8, height = 4, units = c("in")) 

# ── Average across samples ────────────────────────────────────────────────────

filt_mrna <- mrna_df["EEF2", , drop = FALSE] %>% t() %>% as.data.frame()
colnames(filt_mrna) <- "mRNA"
filt_mrna$Sample <- rownames(filt_mrna)
filt_mrna <- filt_mrna %>% mutate(mRNA = log2(mRNA + 1))

filt_prot <- prot_df["EEF2", , drop = FALSE] %>% t() %>% as.data.frame()
colnames(filt_prot) <- "Protein"
filt_prot$Sample <- rownames(filt_prot)
min_value <- min(filt_prot$Protein, na.rm = TRUE)
filt_prot <- filt_prot %>% mutate(Protein = log2(Protein + abs(min_value) + 1))

filt_df <- merge(filt_mrna, filt_prot, by = "Sample")
filt_df <- filt_df %>% mutate(Group = case_when(Sample %in% sh_samp ~ "sh", Sample %in% mh_samp ~ "mh", TRUE ~ "Other"))
filt_df <- filt_df %>% filter(Group != "Other")

p2.dot <- ggplot(filt_df, aes(x = mRNA, y = Protein, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x.npc = "left", label.y.npc = 0.9, size = 5) +
  facet_wrap(~Group) +
  ylim(12, 17) +
  theme_bw() + 
  theme(text = element_text(size = 14))
p2.dot 
ggsave(filename = "02. Protein scatter gene specific.pdf", plot = p2.dot, width = 8, height = 4, units = c("in")) 

# ── Protein enrichment ──────────────────────────────────────────────────────── 

min_protn <- min(prot_df, na.rm = TRUE)

prot_log2t <- prot_df %>% mutate(across(everything(), ~ log2(.x - min_protn + 1)))
prot_clean <- prot_log2t 
prot_clean[is.na(prot_clean)] <- 0 

# Custom gene set 

genesets <- read.xlsx("IDHm_new_genes.xlsx") 

genesets <- genesets %>% 
  dplyr::select(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA) %>% 
  pivot_longer(c(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA), names_to = "source", values_to = "target") %>% 
  mutate(weight = 1) %>%
  filter(target!=" " & target!="") 

# Run decoupler 

res_ssgsea <- decoupleR::run_gsva(mat = prot_clean, 
                                  network = genesets, 
                                  .source ='source', 
                                  .target ='target', 
                                  minsize = 2L, 
                                  method = c("ssgsea")) 

meta_filtr <- meta_df %>% select(hypoxia_class) %>% rownames_to_column(var = "SAMPLE_ID")
res_ssgsea <- res_ssgsea %>% left_join(meta_filtr, by = c("condition" = "SAMPLE_ID")) 

p.box <- ggplot(res_ssgsea, aes(x = hypoxia_class, y = score, fill = hypoxia_class)) +
  geom_boxplot() + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) + 
  stat_compare_means(method = "t.test", label = "p.format") + 
  facet_wrap(~source) 
p.box 
ggsave(filename = "02. Protein RE enrichment.pdf", plot = p.box, width = 2, height = 4, units = c("in")) 
