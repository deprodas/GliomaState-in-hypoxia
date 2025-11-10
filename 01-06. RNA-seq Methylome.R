# Author : Depro Das, Department of Neurosurgery, University Hospital Freiburg, Freiburg, Germany 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(ggplot2)
library(ggpubr)
library(limma)
library(edgeR) 
library(DMRcate)
library(Gviz) 
library(RColorBrewer) 

# ── Methylation data preparation ────────────────────────────────────────────── 

# https://github.com/hamidghaedi/Methylation_Analysis 

# https://f1000research.com/articles/6-2055 

met <- read.table(file = "TCGA-LGG.methylation450.tsv.gz", sep = "\t", header = TRUE) 
met <- met %>% column_to_rownames(var = "Composite.Element.REF")
colnames(met) <- gsub("\\.", "-", gsub("(TCGA)\\.(..)\\.(.{4})\\..*", "\\1-\\2-\\3", colnames(met)))

# Remove probes with NA 

probe.sum <- rowSums(is.na(met)) 
table(probe.sum == 0) 
probe.null <- probe.sum[probe.sum == 0]
met <- met[row.names(met) %in% names(probe.null), ] 

# Remove probes that match chromosomes X and Y 

probe.keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(probe.keep)
met <- met[probe.keep, ] 

# Removing probes that have been demonstrated to map to multiple places in the genome (https://www.tandfonline.com/doi/full/10.4161/epi.23470) 

crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]
met <- met[ -which(row.names(met) %in% crs.reac), ]

# Beta value 

bval <- met

# ── Classify and compare groups ─────────────────────────────────────────────── 

# Metadata (only IDH-mutant) 

meta_df <- read.csv("IDHm_metadata.csv", row.names = 1)
barcode <- meta_df %>% rownames() %>% as.character()
colnames(meta_df)

# Removing samples from meth matrixes 

bval <- bval[, colnames(bval) %in% barcode] 
dim(bval) 
meta_df <- meta_df %>% filter(rownames(.) %in% colnames(bval))

# Making sure about samples in clinical and matrixes and their order

table(colnames(bval) %in% row.names(meta_df))

all(row.names(meta_df) %in% colnames(bval))
all(colnames(bval) %in% row.names(meta_df))

# Plot violin (separate groups)

meta_vlue <- meta_df %>% select(hypoxia_class) %>% rownames_to_column(var = "samples")

bval_long <- bval %>%
  as.data.frame() %>% 
  rownames_to_column("ProbeID") %>%
  pivot_longer(cols = -ProbeID, names_to = "samples", values_to = "Beta") %>%
  left_join(meta_vlue, by = c("samples" = "samples")) 

p.vln_hpx <- ggplot(bval_long, aes(x = hypoxia_class, y = Beta, fill = hypoxia_class)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  theme_minimal() + 
  scale_fill_manual(values = c("Mild_hypoxia" = "skyblue", "Severe_hypoxia" = "salmon"))
p.vln_hpx 
ggsave(filename = "01. Beta value plot-1.pdf", plot = p.vln_hpx, width = 6, height = 4, units = c("in")) 

# ── Differential methylation analysis ───────────────────────────────────────── 

# Converting beta values to m-values

mval <- t(apply(bval, 1, function(x) log2(x/(1-x))))
dim(mval)

all(colnames(mval) %in% colnames(bval))
all(colnames(bval) %in% colnames(mval))

saveRDS(mval, file = "mval.RDS", compress = FALSE)
saveRDS(bval, file = "bval.RDS", compress = FALSE) 

# Grouping variable

meta_df$hypoxia_class <- as.factor(meta_df$hypoxia_class)
meta_df$hypoxia_class <- relevel(meta_df$hypoxia_class, ref = "Mild_hypoxia")

design <- model.matrix(~ hypoxia_class, data = meta_df)

fit1 <- lmFit(mval, design)
fit2 <- eBayes(fit1)

# Extracting significantly methylated probes

deff.meth = topTable(fit2, coef = ncol(design), sort.by = "p",number = nrow(mval), adjust.method = "BY")
deff.meth$ProbeID <- rownames(deff.meth)

# Add probe annotation 

ann450k <- read.delim("HM450.hg38.manifest.gencode.v36.probeMap", stringsAsFactors = FALSE) 
deff.meth <- merge(deff.meth, ann450k, by.x = "ProbeID", by.y = "X.id", all.x = TRUE)

# Set condition 

deff.meth <- deff.meth %>% na.omit()
deff.meth$nlog10 <- -log10(deff.meth$P.Value) 

deff.meth$methylation = ifelse(deff.meth$logFC >= 0.4, 'Hypermethylated', 
                        ifelse(deff.meth$logFC <= -0.4, 'Hypomethylated', 'Unchanged')) 
deff.meth %>% dplyr::count(methylation) 
write.csv(deff.meth, "Results DMG (SH-MH).csv", row.names = FALSE) 

# Top methylated genes 

top <- 10 
top_myth.genes <- bind_rows(deff.meth %>% 
                              filter(methylation == 'Hypermethylated') %>% 
                              arrange(P.Value, desc(abs(logFC))) %>% 
                              head(top), 
                            deff.meth %>% 
                              filter(methylation == 'Hypomethylated') %>% 
                              arrange(P.Value, desc(abs(logFC))) %>% 
                              head(top)) 
# Volcano plot 

vol_met <- ggplot(data = deff.meth, aes(x = logFC, y = nlog10, fill = methylation)) +
  geom_point(size = 1, alpha = 1, pch = 21, stroke = 0.0001) + 
  scale_fill_manual(values = c("#FD841F", "#D6D5A8", "#38E54D")) + 
  geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 0.4, colour="#990000", linetype="dashed") +
  theme(legend.position="none") + 
  theme_bw() +
  theme(legend.position = "none") + 
  ggrepel::geom_text_repel(aes(label = top_myth.genes$gene), data = top_myth.genes, size = 2, force = 10, fontface = "italic", max.overlaps = Inf) 
vol_met 
ggsave(filename = "02. Volcano methylation-1.pdf", plot = vol_met, width = 6, height = 4, units = c("in")) 
