# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(tibble)
library(ggplotify) 
library(Seurat) 
library(CellChat) 
library(patchwork) 
library(viridis) 
library(ComplexHeatmap) 
library(circlize) 
library(reshape2) 
library(reticulate) 
library(NMF) 
library(ggsankey)

# ── Data input and create cellchat ──────────────────────────────────────────── 

sc_data <- readRDS("mIDH total (CCA integrated Seurat 4.3.0).RDS")
sc_inew <- readRDS("mIDH glioma (CCA integrated Seurat 4.3.0).rds") 

# Data manipulation (transfer features from tumor classes to total TME populations)

tme_features <- c('ann_mye_general', 'ann_mye_details') 
tme_metadata <- FetchData(object = sc_data, vars = tme_features) 

tum_features <- c('metaprog_4_class', 'QAD_class', 'metaprog_QAD') 
tum_metadata <- FetchData(object = sc_inew, vars = tum_features) 

tme_metadata <- tme_metadata %>%
  rownames_to_column("barcode") %>%
  left_join(tum_metadata %>% rownames_to_column("barcode") %>% select(barcode, ann_mps.mye_general = metaprog_4_class, ann_qad.mye_general = QAD_class, ann_mqa.mye_general = metaprog_QAD), by = "barcode") %>%
  mutate(ann_mps.mye_general = coalesce(ann_mps.mye_general, ann_mye_general), 
         ann_qad.mye_general = coalesce(ann_qad.mye_general, ann_mye_general), 
         ann_mqa.mye_general = coalesce(ann_mqa.mye_general, ann_mye_general)) %>%
  column_to_rownames("barcode")

tme_metadata <- tme_metadata %>%
  rownames_to_column("barcode") %>%
  left_join(tum_metadata %>% rownames_to_column("barcode") %>% select(barcode, ann_mps.mye_details = metaprog_4_class, ann_qad.mye_details = QAD_class, ann_mqa.mye_details = metaprog_QAD), by = "barcode") %>%
  mutate(ann_mps.mye_details = coalesce(ann_mps.mye_details, ann_mye_details), 
         ann_qad.mye_details = coalesce(ann_qad.mye_details, ann_mye_details), 
         ann_mqa.mye_details = coalesce(ann_mqa.mye_details, ann_mye_details)) %>%
  column_to_rownames("barcode")

tme_metadata %>% dplyr::count(ann_mps.mye_general)

# Add metadata 

sc_data <- AddMetaData(object = sc_data, metadata = tme_metadata[, c("ann_mps.mye_general", 
                                                                     "ann_qad.mye_general", 
                                                                     "ann_mqa.mye_general", 
                                                                     "ann_mps.mye_details", 
                                                                     "ann_qad.mye_details", 
                                                                     "ann_mqa.mye_details")])
# Save total seurat object 

DefaultAssay(sc_data) <- "RNA" 
saveRDS(sc_data, file = "mIDH total (CCA integrated Seurat 4.3.0).RDS") 

# Remove unknown cells 

sc_data@meta.data %>% dplyr::count(ann_mps.mye_general)
sc_data <- subset(sc_data, subset = !ann_mps.mye_general %in% c("Unknown", "Glioma"))

# Create cellchat object  

data.in <- GetAssayData(sc_data, assay = "RNA", layer = "data") 
sc_data <- SetIdent(sc_data, value = sc_data@meta.data$ann_mps.mye_general) 
i.label <- Idents(sc_data) 
meta.df <- data.frame(group = i.label, row.names = names(i.label)) 

cellchat <- createCellChat(object = data.in, meta = meta.df, group.by = "group")

# ── Set up object and database ────────────────────────────────────────────────  

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# Use a subset of cellchatdb for analysis

CellChatDB.sec <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
CellChatDB.use <- CellChatDB.sec 

df_pathway <- as.data.frame(CellChatDB.use[["interaction"]])
df_pathway %>% dplyr::count(pathway_name)

cellchat@DB <- CellChatDB.use 

# ── Processing and analysis ───────────────────────────────────────────────────   

cellchat <- subsetData(cellchat) 

future::plan("multisession", workers = 4) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human) 

# Compute the communication probability 

options(future.globals.maxSize = 8000 * 1024^2) 
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = T) 
cellchat <- filterCommunication(cellchat, min.cells = 10) 

# Compute the communication probability on signaling pathway level 

cellchat <- computeCommunProbPathway(cellchat) 
cellchat <- aggregateNet(cellchat) 
cellchat <- netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat, file = "Cellchat ann_mps.mye_general pop-T (CCA integrated Seurat 4.3.0).RDS") 

# ── Predict general principles ──────────────────────────────────────────────── 

cellchat <- readRDS("Cellchat ann_mps.mye_general pop-T (CCA integrated Seurat 4.3.0).RDS")

# Total interaction circular plot 

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)

pdf("1. Total Interactions pop-T (Circle).pdf", width = 10, height = 8)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# Total interaction heatmap 

heat_con.i <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues") 
heat_wgt.i <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues") 
heat_tot.i <- as.ggplot(heat_con.i) + as.ggplot(heat_wgt.i) 
heat_tot.i
ggsave(filename = "2. Total Interactions pop-T (Heatmap).pdf", plot = heat_tot.i, width = 9.25, height = 4.25, units = c("in")) 

# Specific cell-cell interaction / custom heatmap (count and weight) 

customHeatmap <- function(cellchat, matrix_name = "weight", heatmap_colors = NULL) {
  count.tot <- cellchat@net[[matrix_name]]
  count.tot <- count.tot %>% as.data.frame() %>% dplyr::select(c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'TCells', 'Oligo', 'Endo', 'Pericytes', 'BCells'))
  count.tot <- count.tot %>% dplyr::filter(rownames(.) %in% c('AC_like', 'NPC_like', 'OC_like', 'Cycling'))
  
  row_sums <- rowSums(count.tot, na.rm = TRUE)
  col_sums <- colSums(count.tot, na.rm = TRUE)
  row_anno <- rowAnnotation(Sum = anno_barplot(row_sums, gp = gpar(fill = "grey"), width = unit(2, "cm")))
  col_anno <- HeatmapAnnotation(Sum = anno_barplot(col_sums, gp = gpar(fill = "grey"), height = unit(2, "cm"))) 
  
  if (is.null(heatmap_colors)) {
    heatmap_colors <- colorRamp2(c(min(count.tot, na.rm = TRUE), mean(count.tot, na.rm = TRUE), max(count.tot, na.rm = TRUE)), c("white", "lightblue", "blue"))
  }
  
  wg_h <- Heatmap(as.matrix(count.tot),
                  name = matrix_name, 
                  col = heatmap_colors,
                  top_annotation = col_anno,
                  right_annotation = row_anno, 
                  border = TRUE) 
  wg_h.g <- as.ggplot(wg_h)
  return(wg_h.g)
}

cust_ct <- customHeatmap(cellchat, 
                         matrix_name = "count", 
                         heatmap_colors = colorRamp2(c(0, 2.5, 5, 7.5, 10), c("white", "#caf0f8", "#00b4d8", "#0077b6", "#03045e"))) 

cust_wg <- customHeatmap(cellchat, 
                         matrix_name = "weight",
                         heatmap_colors = colorRamp2(c(0.0, 0.00375, 0.0075, 0.01125, 0.015), c("white", "#caf0f8", "#00b4d8", "#0077b6", "#03045e")))

cust_h <- cust_ct + cust_wg 
cust_h 
ggsave(filename = "3. Specific interactions pop-T (Heatmap).pdf", plot = cust_h, width = 12, height = 4, units = c("in"))


# Visualize dominant senders (sources) and receivers (targets) in a 2D space

scat.tot <- netAnalysis_signalingRole_scatter(cellchat)
scat.spc <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("EGF"))
scat.all <- scat.tot + scat.spc 
scat.all
ggsave(filename = "4. Dominant senders and receivers pop-T.pdf", plot = scat.all, width = 7, height = 2.75, units = c("in")) 

# Identify signals contributing the most to outgoing or incoming signaling of certain cell groups

heat_out.r <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.heatmap = "Blues", width = 6, height = 14)
heat_inn.r <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.heatmap = "Blues", width = 6, height = 14)
heat_tot.r <- as.ggplot(heat_out.r) + as.ggplot(heat_inn.r) 
heat_tot.r
ggsave(filename = "5. Outgoing and incoming signaling pop-T (all).pdf", plot = heat_tot.r, width = 10, height = 10, units = c("in")) 

# ── Pathway-specific roles ────────────────────────────────────────────────────  

# Network centrality scores 

all.pathways <- cellchat@netP$pathways

for (pw in all.pathways) {
  pdf(paste0("6. Pathway specific roles pop-T ", pw, ".pdf"), width = 6, height = 4) 
  
  netAnalysis_signalingRole_network(cellchat, 
                                    signaling = pw, 
                                    color.heatmap = "Blues",
                                    cluster.cols = TRUE,
                                    width = 8,
                                    height = 2.25,
                                    font.size = 10) 
  dev.off()
}

# Communication probability at pathway level 

com.path <- cellchat@netP[["pathways"]]
com.prob <- cellchat@netP[["prob"]]
df_commp <- melt(com.prob, varnames = c("source", "target", "pathway_index"), value.name = "probability")
df_commp$pathway <- com.path[df_commp$pathway_index] 

sources.use <- c('Microglia', 'Macrophage', 'DCs', 'MDSC')
targets.use <- c('AC_like', 'NPC_like', 'OC_like', 'Cycling')

for (pw in unique(df_commp$pathway)) {
  df_wide <- df_commp %>%
    dplyr::filter(pathway == pw, source %in% sources.use, target %in% targets.use) %>%
    pivot_wider(id_cols = source, names_from = target, values_from = probability) %>%
    column_to_rownames("source") 
  
  df_wide[] <- lapply(df_wide, function(x) as.numeric(unlist(x)))
  
  if (all(is.na(df_wide)) || sum(df_wide, na.rm = TRUE) == 0) {
    message(paste("Skipping pathway:", pw))
    next
  }
  
  df_wide <- t(df_wide) 
  row_sums <- rowSums(df_wide, na.rm = TRUE)
  col_sums <- colSums(df_wide, na.rm = TRUE)
  row_anno <- rowAnnotation(Sum = anno_barplot(row_sums, gp = gpar(fill = "grey"), width = unit(2, "cm")))
  col_anno <- HeatmapAnnotation(Sum = anno_barplot(col_sums, gp = gpar(fill = "grey"), height = unit(2, "cm"))) 
  mat_vals <- as.matrix(df_wide)
  mat_vals <- mat_vals[!is.na(mat_vals)]
  break_ht <- quantile(mat_vals, probs = c(0, 0.25, 0.5, 0.75, 1)) 
  comm_h <- Heatmap(as.matrix(df_wide),
                    name = "Communication probability",
                    col = colorRamp2(break_ht, c("white", "#caf0f8", "#00b4d8", "#0077b6", "#03045e")),
                    top_annotation = col_anno,
                    right_annotation = row_anno,
                    border = TRUE,
                    column_title = pw) 
  comm_h.g <- as.ggplot(comm_h)
  ggsave(filename = paste0("7. Communication probability pop-T ", pw, ".pdf"), plot = comm_h.g, width = 5.75, height = 4.5, units = "in")
} 

# ── Ligand-receptor pair (more focused on one pathway) ──────────────────────── 

pathways.show <- "EGF" 

# Communication probability at ligand-receptor level 

lr_p.1 <- netVisual_bubble(cellchat, 
                           signaling = pathways.show, 
                           sources.use = c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'Oligo', 'Endo', 'Pericytes'),
                           targets.use = c('AC_like', 'NPC_like', 'OC_like', 'Cycling'), 
                           remove.isolate = FALSE) 
lr_p.1
ggsave(filename = paste0("8. Communication probability lr pop-T ", pathways.show, ".pdf"), plot = lr_p.1, width = 12, height = 2.4, units = c("in"))


# Manual plot (Communication probability) 

str(cellchat@net$pval)

prob_array <- cellchat@net$prob
pval_array <- cellchat@net$pval
class(prob_array) 

sources <- dimnames(prob_array)[[1]]
targets <- dimnames(prob_array)[[2]]
intract <- dimnames(prob_array)[[3]]

# Get array indices and convert that to an explicit index matrix

df_indices <- as.data.frame(which(!is.na(prob_array), arr.ind = TRUE))
index_mtrx <- as.matrix(df_indices[, c("dim1", "dim2", "dim3")])

# Add values to the index matrix as much as i want 

df_indices$comprob <- prob_array[index_mtrx]
df_indices$pvalues <- pval_array[index_mtrx]
df_indices$sources <- sources[df_indices$dim1]
df_indices$targets <- targets[df_indices$dim2]
df_indices$intract <- intract[df_indices$dim3]

# Filter specific interactions 

df_specific <- df_indices %>% 
  dplyr::filter(intract %in% c("HBEGF_EGFR", "AREG_EGFR", "TGFA_EGFR") & 
                  sources %in% c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'Oligo', 'Endo', 'Pericytes') & 
                  targets %in% c('AC_like', 'NPC_like', 'OC_like', 'Cycling')) 

df_specific$new_pvalue <- ifelse(df_specific$pvalues == 0, 1e-5, df_specific$pvalues)

lr_p.2 <- ggplot(df_specific, aes(x = targets, y = intract, color = comprob, size = -log10(new_pvalue))) +
  geom_point(alpha = 1) + 
  scale_colour_gradient2(low = "#4dac26", mid = "white", high = "#d01c8b", midpoint = 0.000135) +
  facet_wrap(~sources, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
lr_p.2

lr_bar.1 <- ggplot(df_specific, aes(x = targets, y = comprob, fill = targets)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~sources, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") 
lr_bar.1

lr_comb <- lr_bar.1 / lr_p.2 + plot_layout(heights = c(1.5, 3))
lr_comb
ggsave(filename = paste0("9. Communication probability lr manual pop-T ", pathways.show, ".pdf"), plot = lr_comb, width = 14, height = 3.5, units = c("in")) 


# Signaling contribution of pairs 

lr_cont <- netAnalysis_contribution(cellchat, signaling = "EGF") 
lr_cont
ggsave(filename = paste0("10. Contribution lr pop-T ", pathways.show, ".pdf"), plot = lr_cont, width = 4, height = 3, units = c("in")) 


# Plot gene expression 

lr_exp <- plotGeneExpression(cellchat, signaling = "EGF", enriched.only = TRUE, type = "violin")
lr_exp
ggsave(filename = paste0("11. Expression lr pop-T ", pathways.show, ".pdf"), plot = lr_exp, width = 4.5, height = 3.5, units = c("in")) 


# ── Global communication patterns ───────────────────────────────────────────── 

# Infer the number of patterns

selectK(cellchat, pattern = "outgoing") 
selectK(cellchat, pattern = "incoming") 

nPatterns_ot = 2
nPatterns_in = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns_ot)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns_in)

# River plots

g.patt_ot <- netAnalysis_river(cellchat, pattern = "outgoing")
g.patt_in <- netAnalysis_river(cellchat, pattern = "incoming")
g.patt_rp <- g.patt_ot + g.patt_in
g.patt_rp
ggsave(filename = "14. Communication patterns all (river plot).pdf", plot = g.patt_rp, width = 14, height = 6, units = c("in")) 

# Dot plots 

d.patt_ot <- netAnalysis_dot(cellchat, pattern = "outgoing")
d.patt_in <- netAnalysis_dot(cellchat, pattern = "incoming")
d.patt_rp <- d.patt_ot + d.patt_in
d.patt_rp
ggsave(filename = "15. Communication patterns all (dot plot).pdf", plot = d.patt_rp, width = 14, height = 3.5, units = c("in")) 

# ── Uniform plot based on outgoing and incoming signaling ───────────────────── 

# Custom alluvial plot 

netP.ot <- cellchat@netP[["pattern"]][["outgoing"]] 
data.ot <- netP.ot[["data"]]
cell.ot <- netP.ot[["pattern"]][["cell"]] 
sigl.ot <- netP.ot[["pattern"]][["signaling"]]

netP.in <- cellchat@netP[["pattern"]][["incoming"]] 
data.in <- netP.in[["data"]]
cell.in <- netP.in[["pattern"]][["cell"]] 
sigl.in <- netP.in[["pattern"]][["signaling"]]

data.new.ot <- as.data.frame(data.ot) %>% rownames_to_column(var = "celltype") %>% pivot_longer(cols = -celltype, names_to = "pathways", values_to = "outgoing_score")
data.new.in <- as.data.frame(data.in) %>% rownames_to_column(var = "celltype") %>% pivot_longer(cols = -celltype, names_to = "pathways", values_to = "incoming_score")

df.ot_send <- data.new.ot %>% 
  dplyr::filter(celltype %in% c('Microglia', 'Macrophage', 'DCs', 'MDSC', 'Oligo', 'Endo', 'Pericytes')) %>% 
  dplyr::rename(senders = celltype) 

df.in_recv <- data.new.in %>% 
  dplyr::filter(celltype %in% c('AC_like', 'NPC_like', 'OC_like', 'Cycling')) %>% 
  dplyr::rename(receivers = celltype) 

# Meaningful and shared interactions 

df.ot_send <- df.ot_send %>% filter(outgoing_score > quantile(df.ot_send$outgoing_score, 0.5, na.rm = TRUE))
df.in_recv <- df.in_recv %>% filter(incoming_score > quantile(df.in_recv$incoming_score, 0.5, na.rm = TRUE))

shared_path <- intersect(df.ot_send$pathways, df.in_recv$pathways)

df.ot_send <- df.ot_send %>% filter(pathways %in% shared_path)
df.in_recv <- df.in_recv %>% filter(pathways %in% shared_path)

# Common format 

df_sig.true <- df.ot_send %>% inner_join(df.in_recv, by = "pathways") 

# Try several scoring 

df_multiplicative <- df_sig.true %>% 
  mutate(score = outgoing_score * incoming_score)

df_geometric.mean <- df_sig.true %>% 
  mutate(score = sqrt(outgoing_score * incoming_score)) 

df_maximum <- df_sig.true %>% 
  mutate(score = pmax(outgoing_score, incoming_score))

df_normal.mean <- df_sig.true %>% 
  mutate(score = (outgoing_score + incoming_score)/2) 

df_harmonic.mean <- df_sig.true %>% 
  mutate(score = 2 * outgoing_score * incoming_score / (outgoing_score + incoming_score))

# Plot meaningful interactions (sankey plot)

df_sig.common <- df_harmonic.mean

df_sig.sankey <- df_sig.common %>% 
  filter(score > median(df_sig.common$score, na.rm = TRUE)) %>% 
  select(c(senders, pathways, receivers, score)) %>% 
  make_long(senders, pathways, receivers, value = score) 

p.lr_snk <- ggplot(df_sig.sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, value = value, fill = factor(node))) +
  geom_sankey(flow.alpha = 0.7) +
  geom_sankey_label(aes(label = node), size = 3, color = 1, fill = "white") +
  theme_sankey(base_size = 12) + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text  = element_blank()) 
p.lr_snk 
ggsave(filename = "16. Meaningful interactions sankey.pdf", plot = p.lr_snk, width = 5, height = 4, units = c("in")) 
