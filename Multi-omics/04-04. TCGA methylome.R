# Author : Depro Das, 3DBM Lab, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

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
