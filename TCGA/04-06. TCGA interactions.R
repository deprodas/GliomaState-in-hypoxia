# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(BulkSignalR) 
library(igraph) 
library(STexampleData) 
library(doParallel) 
library(tidyverse) 
library(openxlsx) 
library(liana) 
library(ggplot2)
library(ggsankey) 
library(stringr)

# ── Set up geneset & database (constant) ────────────────────────────────────── 

# Custom pathways 

net_prog <- read.xlsx("mIDH programs_depro.xlsx") 
net_prog <- net_prog %>%
  pivot_longer(cols = everything(), names_to = "Reactome name", values_to = "Gene name") %>%
  filter(!is.na(`Gene name`)) %>%
  distinct()

pathwy_id <- net_prog %>%
  distinct(`Reactome name`) %>%
  mutate(`Reactome ID` = paste0("Prog_", row_number()))

net_prog <- net_prog %>% left_join(pathwy_id, by = "Reactome name")
net_prog <- net_prog %>% mutate(across(everything(), as.character))
sapply(net_prog, is.list) 
net_prog <- as.data.frame(net_prog)

resetPathways(dataframe = net_prog, resourceName = "Reactome")

# Custom LR database (Liana)

lrdf_liana <- liana::select_resource("Consensus")
lrdb_liana <- lrdf_liana$Consensus %>%
  select(ligand = source_genesymbol, receptor = target_genesymbol) %>%
  separate_rows(receptor, sep = "_") %>%
  distinct() %>%
  as.data.frame()

resetLRdb(lrdb_liana, switch = TRUE)

# ── Prepare inputs (variable) ───────────────────────────────────────────────── 

# Metadata (only IDH-mutant) 

metadata <- read.csv("IDHm_metadata.csv", row.names = 1) 
metadata$SUBTYPE_GRADE <- interaction(metadata$SUBTYPE, metadata$GRADE)
metadata <- metadata[!is.na(metadata$SUBTYPE_GRADE), ] 
metadata <- metadata %>% mutate(sub_grade_hpx = paste(SUBTYPE_GRADE, hypoxia_class, sep = "_"))
metadata %>% dplyr::count(sub_grade_hpx)

# Count data (only IDH-mutant) 

htseq_counts <- read.csv("TCGA-PanCan IDHm-glioma raw-count.csv", row.names = 1) 
colnames(htseq_counts) <- gsub('\\.', '-', colnames(htseq_counts)) 
htseq_counts <- htseq_counts[, rownames(metadata)] 

# Convert NA values to 0 in count matrix

htseq_counts[is.na(htseq_counts)] <- 0
all(is.numeric(htseq_counts))
all(rownames(metadata) %in% colnames(htseq_counts)) 
all(colnames(htseq_counts) %in% rownames(metadata)) 

# ── Run analysis ────────────────────────────────────────────────────────────── 

groups <- unique(metadata$sub_grade_hpx)

lr_result_list <- list()

for (g in groups) {
  message("Running: ", g)
  meta_sub <- metadata %>% filter(sub_grade_hpx == g)
  counts_sub <- htseq_counts[, rownames(meta_sub), drop = FALSE]
  counts_sub[is.na(counts_sub)] <- 0 
  
  closeAllConnections()
  cl_mk <- makeCluster(1)
  registerDoParallel(cl_mk)
  bsrdm <- BSRDataModel(counts = counts_sub)
  bsrdm <- learnParameters(bsrdm, quick = TRUE, max.pw.size = 1000)
  bsrif <- BSRInference(bsrdm, min.cor = 0.00, reference = "REACTOME") 
  lr_df <- LRinter(bsrif)
  lr_df$sub_grade_hpx <- g 
  lr_result_list[[g]] <- lr_df
  stopCluster(cl_mk)
}

lr_all.res <- bind_rows(lr_result_list)
write.xlsx(lr_all.res, file = "Result LR (mIDH programs_depro).xlsx")

# ── Plot ──────────────────────────────────────────────────────────────────────  

lr_all.res <- read.xlsx("Result LR (mIDH programs_depro).xlsx") 
lr_fil.res <- lr_all.res %>%  
  filter(L %in% c('HBEGF') & R %in% c('EGFR') & pw.name %in% c('EGFR_all')) %>%  
  filter(str_detect(sub_grade_hpx, "Severe_hypoxia"))

lr_sank_df <- lr_fil.res %>% make_long(sub_grade_hpx, L, R, value = rank.corr) 

pt_lr_sank <- ggplot(lr_sank_df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, value = value)) +
  geom_sankey(flow.alpha = 0.6, node.color = "black") +
  geom_sankey_label(aes(label = node), size = 3, color = "black", fill = "white") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_blank()) 
pt_lr_sank 
ggsave(pt_lr_sank, filename = "01. Sankey selective.pdf", width = 6, height = 3, units = c("in"))
