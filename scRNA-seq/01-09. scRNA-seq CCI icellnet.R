# Author : Depro Das, 3DBM Lab, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(R.utils)
library(dplyr) 
library(tidyverse) 
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(ggplotify)
library(openxlsx) 
library(icellnet)
library(gridExtra)
library(Seurat) 
library(ComplexHeatmap) 
library(circlize) 

# ── Data input ────────────────────────────────────────────────────────────────  

sc_data <- readRDS("mIDH total (CCA integrated Seurat 4.3.0).RDS")
sc_data <- subset(sc_data, subset = !ann_mps.mye_general %in% c("Unknown", "Glioma"))
DimPlot(sc_data, reduction = 'umap', group.by = 'seurat_clusters', label = T)

# Microglia-3227, Macrophage-586, MDSC-109, DCs-55 

metadata <- sc_data@meta.data 
metadata %>% count(ann_mps.mye_general)

# ── Icellnet database ───────────────────────────────────────────────────────── 

# Set up database 

db <- as.data.frame(read.csv(curl::curl(url = "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
db.name.couple <- name.lr.couple(db, type = "Family")
db.name.subfam <- name.lr.couple(db, type = "Subfamily")
head(db.name.couple)

write.xlsx(db, file = "ICELLNETdb.xlsx")

# Taking into account the total nb of cells in each cluster

filter.percents <- 0 
Idents(sc_data) <- sc_data$ann_mps.mye_general 
average.cleaned <- sc.data.cleaning(object = sc_data, 
                                    db = db, 
                                    filter.perc = filter.percents, 
                                    save_file = F, 
                                    path = "~/Desktop/Research Projects (3DBM)/1. mIDH (hypoxia)/2. Analysis/02. Analysis (Single-cell)/09. Cell-cell (Icellnet)/", 
                                    force.file = F)
saveRDS(average.cleaned, file = "Icellnet ann_mps.mye_general (CCA integrated Seurat 4.3.0).RDS")

# Select senders and receivers   

average.cleaned <- readRDS("Icellnet ann_mps.mye_general (CCA integrated Seurat 4.3.0).RDS")

# Select senders and receivers   

data.icell <- as.data.frame(gene.scaling(as.data.frame(average.cleaned), n = 1, db = db))
colnames(data.icell)
sendercell <- as.data.frame(data.icell[ , c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'Oligo', 'Endo', 'TCells', 'Pericytes', 'BCells')], row.names = rownames(data.icell))
selections <- c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'Oligo', 'Endo', 'TCells', 'Pericytes', 'BCells')

# Compute communication scores 

score.computation.1 <- icellnet.score(direction = "in", 
                                      PC.data = sendercell, 
                                      CC.data = as.data.frame(data.icell[ , c('AC_like')], row.names = rownames(data.icell)),  
                                      PC = selections, 
                                      CC.type = "RNAseq", 
                                      PC.type = "RNAseq",  
                                      db = db)
scores_df <- as.data.frame(score.computation.1[[1]])
lr_result <- score.computation.1[[2]]
write.xlsx(lr_result, file = "LR scores-raw (immune to glioma).xlsx")

# Plot direction and network 

network.plot_in <- network.create(icn.score = scores_df[1], 
                                  scale = c(round(min(scores_df)-1), round(max(scores_df))+1), 
                                  direction = "in") 
network.plot_in
ggsave(network.plot_in, filename = "1. Signaling direction.pdf", width = 7, height = 5, units = c("in")) 

# Filter out 'others' or 'NA' values 

db.name.couple %>% as.data.frame() %>% dplyr::count(Family)
db.name.subfam %>% as.data.frame() %>% dplyr::count(Subfamily)

db.new.couple <- db.name.couple %>% na.omit() %>% as.data.frame() %>% column_to_rownames(var = "Pair")
db.new.subfam <- db.name.subfam %>% na.omit() %>% as.data.frame() %>% column_to_rownames(var = "Pair")

db.new.couple %>% as.data.frame() %>% dplyr::count(Family)
db.new.subfam %>% as.data.frame() %>% dplyr::count(Subfamily)

# ── Family contribution ─────────────────────────────────────────────────────── 

# Visualization of contribution of family of molecules to communication scores 

# Define the y axis range of the barplot 

ymax = round(max(scores_df))+1 

class(db.name.subfam)
class(lr_result)

# Exclude unspecific / NA / others interactions 

lr_subfam <- lr_result %>% as.data.frame() %>% dplyr::filter(rownames(lr_result) %in% rownames(db.new.subfam)) %>% dplyr::select(!c("BCells")) 
lr_subfam <- na.omit(lr_subfam)
db_subfam <- db.name.subfam %>% na.omit() %>% as.data.frame() 

ht.family <- LR.family.score(lr = lr_subfam, 
                             db.couple = db_subfam, 
                             plot = "heatmap", 
                             title = "Family contribution") 
ht.family
ht.family <- ggplotify::as.ggplot(ht.family)
ggsave(ht.family, filename = "2. Icellnet family contribution.pdf", width = 6, height = 12, units = c("in")) 

# ── Overall L/R interactions ────────────────────────────────────────────────── 

# Most contributing and different interactions 

gg_con <- LR.heatmap(lr = lr_result, thresh = 0, topn = 20, sort.by = "sum", title = "Most contributing interactions")
gg_var <- LR.heatmap(lr = lr_result, thresh = 0, topn = 20, sort.by = "var", title = "Most different interactions")
gg_int <- gg_con + gg_var
gg_int
ggsave(gg_int, filename = "3. Most contributing and different interactions (ggplot heatmap).pdf", width = 9.75, height = 6, units = c("in"))

# Plot differently  

pr_con <- LR.selection(lr = lr_result, thresh = 0 , topn = 20 , sort.by = "sum")
pr_var <- LR.selection(lr = lr_result, thresh = 0 , topn = 20 , sort.by = "var")

col_fu <- colorRamp2(c(0, 100), c("white", "red")) 

ht_con <- Heatmap(as.matrix(pr_con), 
                  cluster_rows = T, 
                  cluster_columns = T, 
                  clustering_method_rows = "ward.D", 
                  name = "Score", 
                  clustering_distance_rows = "euclidean", 
                  show_column_dend = T, 
                  show_row_dend = T, 
                  show_column_names = T, 
                  show_row_names = T,
                  column_title = "Most contributing interactions", 
                  col = col_fu) 
ht_con 
ht_con <- ggplotify::as.ggplot(ht_con)

ht_var <- Heatmap(as.matrix(pr_var), 
                  cluster_rows = T, 
                  cluster_columns = T, 
                  clustering_method_rows = "ward.D", 
                  name = "Score", 
                  clustering_distance_rows = "euclidean", 
                  show_column_dend = T, 
                  show_row_dend = T, 
                  show_column_names = T, 
                  show_row_names = T,
                  column_title = "Most different interactions", 
                  col = col_fu) 
ht_var
ht_var <- ggplotify::as.ggplot(ht_var)

ht_int <- ht_con + ht_var
ht_int 
ggsave(ht_int, filename = "3. Most contributing and different interactions (complex heatmap).pdf", width = 12, height = 7.5, units = c("in"))

# ── Specific family interactions ────────────────────────────────────────────── 

lr_pair <- lr_result %>% as.data.frame() %>% dplyr::filter(rownames(lr_result) %in% rownames(db.new.subfam)) 
lr_pair <- cbind(db.new.subfam, lr_pair)
lr_pair <- na.omit(lr_pair)
lr_pair <- lr_pair %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Signaling") %>% 
  filter(Subfamily %in% c("EGFR family") & if_any(everything(), ~ grepl("/ EGFR", ., ignore.case = TRUE)))

lr_heat <- lr_pair
lr_heat <- lr_heat %>% column_to_rownames(var = "Signaling") 

set.seed(123)
ht_lrsp <- Heatmap(lr_heat[,2:7], 
                   col = colorRamp2(c(0, 60), c("white", "red")),
                   column_title = "Most contributing interactions in EGFR family", 
                   border = TRUE) 
ht_lrsp
ht_lrsp <- as.ggplot(ht_lrsp)
ggsave(ht_lrsp, filename = "4. Icellnet pairs (egfr family).pdf", width = 4.5, height = 4.4, units = c("in"))

# ── Specific L/R interactions ───────────────────────────────────────────────── 

lr_custom_heatmap <- function(data.icell, db, int_name, output_prefix = "Communication probability") {
  lr_pt <- LR.viz(data = data.icell, db = db, int = int_name, plot = TRUE)
  lr_df <- lr_pt[["data"]] %>% pivot_wider(id_cols = c("Ligand"), names_from = c("Receptor"), values_from = "value") %>% column_to_rownames(var = "Ligand") 
  lr_df[] <- lapply(lr_df, function(x) as.numeric(unlist(x))) 
  lr_df <- lr_df %>%
    dplyr::select(c('AC_like', 'Cycling', 'NPC_like', 'OC_like')) %>%
    dplyr::filter(rownames(.) %in% c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'Oligo', 'Endo', 'TCells', 'Pericytes', 'BCells'))
  
  comm_h <- Heatmap(as.matrix(lr_df),
                    column_title = int_name,
                    name = "Communication probability",
                    col = colorRamp2(c(60, 70, 80, 90, 100), c("white", "#caf0f8", "#00b4d8", "#0077b6", "#03045e")),
                    border = TRUE) 
  comm_h_g <- as.ggplot(comm_h) 
  
  filename <- paste0(output_prefix, " ", gsub(" / ", "-", int_name), ".pdf")
  ggsave(filename = filename, plot = comm_h_g, width = 5.75, height = 4.5, units = "in") 
  return(list(heatmap = comm_h, plot = comm_h_g, data = lr_df, filename = filename))
}

lr_custom_heatmap(data.icell = data.icell, 
                  db = db, 
                  int_name = "HBEGF / EGFR", 
                  output_prefix = "6. Communication probability") 

lr_custom_heatmap(data.icell = data.icell, 
                  db = db, 
                  int_name = "AREG / EGFR", 
                  output_prefix = "6. Communication probability") 

lr_custom_heatmap(data.icell = data.icell, 
                  db = db, 
                  int_name = "TGFA / EGFR", 
                  output_prefix = "6. Communication probability") 

