# Author : Depro Das, 3DBM Lab, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library("stringr") 
library(maftools)
library(GenVisR) 
library(data.table) 
library(ggplot2)
library(ggpubr) 
library(patchwork) 
library(viridis) 
library(maftools)

# ── Data input and classification (mutational) ──────────────────────────────── 

# Metadata (only IDH-mutant) 

meta_df <- read.csv("IDHm_metadata.csv", row.names = 1)
meta_df %>% dplyr::count(hypoxia_class) 

# Separate groups  

mh_samp <- meta_df %>% filter(hypoxia_class == "Mild_hypoxia") %>% rownames() %>% as.character()
sh_samp <- meta_df %>% filter(hypoxia_class == "Severe_hypoxia") %>% rownames() %>% as.character()

# Mutational data read and classify  

data_mut <- read.delim("data_mutations.txt", header = TRUE, stringsAsFactors = FALSE)
data_mut$Tumor_Sample_Barcode <- gsub("-[0-9A-Z]+$", "", data_mut$Tumor_Sample_Barcode)

data_mut_sh <- subset(data_mut, Tumor_Sample_Barcode %in% sh_samp)
data_mut_mh <- subset(data_mut, Tumor_Sample_Barcode %in% mh_samp)

# ── Data analysis (mutational) ────────────────────────────────────────────────  

# Convert into maf format 

maf_sh <- read.maf(maf = data_mut_sh)
maf_mh <- read.maf(maf = data_mut_mh)

# Plot summary, oncoplots, variant allele frequencies, titv, 

maf_list <- list(SH = maf_sh, MH = maf_mh)

for (name in names(maf_list)) {
  maf_obj <- maf_list[[name]] 
  
  pdf(paste0("01. ", name, " summary.pdf"), width = 10, height = 6)
  plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  pdf(paste0("02. ", name, " oncoplot.pdf"), width = 8, height = 7)
  oncoplot(maf = maf_obj, top = 20)
  dev.off() 
  
  pdf(paste0("03. ", name, " VAF boxplot.pdf"), width = 6, height = 4)
  plotVaf(maf = maf_obj, vafCol = NULL, top = 20)
  dev.off() 
  
  titv_obj <- titv(maf = maf_obj, plot = FALSE, useSyn = TRUE)
  pdf(paste0("04. ", name, " TiTv plot.pdf"), width = 8, height = 6)
  plotTiTv(res = titv_obj)
  dev.off() 
  
  pdf(paste0("05. ", name, " TCGAcompare.pdf"), width = 10, height = 8)
  tcgaCompare(maf = maf_obj, cohortName = name, logscale = TRUE, capture_size = 50)
  dev.off()
}

# Tumor mutation burden - two groups 

mh_summary <- getSampleSummary(maf_mh)
mh_summary$cohort <- "Mild_hypoxia"

sh_summary <- getSampleSummary(maf_sh)
sh_summary$cohort <- "Severe_hypoxia"

burden_dfs <- rbind(mh_summary, sh_summary)

p.tmb <- ggplot(burden_dfs, aes(x = cohort, y = total, fill = cohort)) +
  geom_violin(trim = FALSE, alpha = 1, width = 0.8) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") + 
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black", fatten = 2, size = 0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ylim(0, 120)
p.tmb
ggsave(p.tmb, file = "06. TMB violin.pdf", width = 3.5, height = 4, units = "in") 


# Differentially mutated genes 

# Need gene summary 

gene_summary_sh <- getGeneSummary(maf_sh)
gene_summary_mh <- getGeneSummary(maf_mh)

gene_compare <- merge(gene_summary_sh[, c("Hugo_Symbol", "AlteredSamples")],
                      gene_summary_mh[, c("Hugo_Symbol", "AlteredSamples")],
                      by = "Hugo_Symbol",
                      suffixes = c("_SH", "_MH"))

n_sh <- nrow(maf_sh@clinical.data)
n_mh <- nrow(maf_mh@clinical.data)

# Fisher's exact test per gene (same result with chi-square test) 

gene_compare$pvalue <- apply(gene_compare, 1, function(row) {
  tab <- matrix(c(as.numeric(row["AlteredSamples_SH"]), n_sh - as.numeric(row["AlteredSamples_SH"]),
                  as.numeric(row["AlteredSamples_MH"]), n_mh - as.numeric(row["AlteredSamples_MH"])), 
                nrow = 2, byrow = TRUE)
  fisher.test(tab)$p.value
})

gene_compare$adj_pvalue <- p.adjust(gene_compare$pvalue, method = "BH")

# Get differentially mutated genes and plot  

diff_genes <- gene_compare[gene_compare$pvalue < 0.05, ]
diff_genes <- diff_genes[order(diff_genes$pvalue), ]

diff_genes$freq_diff <- (diff_genes$AlteredSamples_SH / n_sh) - (diff_genes$AlteredSamples_MH / n_mh) 

top_diffdf <- head(diff_genes[order(abs(diff_genes$freq_diff), decreasing = TRUE), ], 20)

p.dmg <- ggplot(top_diffdf, aes(x = reorder(Hugo_Symbol, freq_diff), y = freq_diff, fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_viridis() + 
  theme_minimal(base_size = 12)
p.dmg 
ggsave(p.dmg, file = "07. DMG barplot.pdf", width = 4, height = 4, units = "in")

# Co-barplots - two groups 

genes_ord <- as.character(top_diffdf$Hugo_Symbol)
genes_ord <- factor(genes_ord, levels = rev(genes_ord))

pdf("08. DMG co-barplot.pdf", width = 6, height = 4) 
coBarplot(m1 = maf_sh, 
          m2 = maf_mh, 
          genes = genes_ord, 
          m1Name = "Severe hxpoxia", 
          m2Name = "Mild hypoxia") 
dev.off()

# Forest plot - two groups 

sh.mh_comp <- mafCompare(m1 = maf_sh, 
                         m2 = maf_mh, 
                         m1Name = "Severe hxpoxia", 
                         m2Name = "Mild hypoxia")

sh.mh_dfs <- sh.mh_comp[["results"]]
sh.mh_top <- sh.mh_dfs %>% filter(Hugo_Symbol %in% dput(as.character(top_diffdf$Hugo_Symbol)))

f.dmg <- ggplot(sh.mh_top, aes(x = Hugo_Symbol, y = or, ymin = ci.low, ymax = ci.up)) +
  geom_pointrange(aes(color = or > 1), size = 1) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  coord_flip() + 
  scale_y_log10() + 
  scale_color_manual(values = c("darkblue", "darkred"), guide = FALSE) + 
  theme_minimal(base_size = 14)
f.dmg
ggsave(f.dmg, file = "09. Forest plot.pdf", width = 4, height = 4, units = "in")

# VAFs of differentially mutated genes 

vaf_sh <- maf_sh@data %>%
  dplyr::mutate(t_vaf = t_alt_count / (t_alt_count + t_ref_count)) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, t_vaf) %>%
  mutate(group = "Severe hypoxia")

vaf_mh <- maf_mh@data %>%
  dplyr::mutate(t_vaf = t_alt_count / (t_alt_count + t_ref_count)) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, t_vaf) %>%
  mutate(group = "Mild hypoxia")

vaf_df <- bind_rows(vaf_sh, vaf_mh)
vaf_df <- vaf_df %>% filter(Hugo_Symbol %in% top_diffdf$Hugo_Symbol)

vaf.dmg <- ggplot(vaf_df, aes(x = group, y = t_vaf, fill = group)) +
  geom_violin(trim = FALSE, alpha = 1, width = 0.8) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black", size = 0.5) +
  stat_compare_means(aes(label = signif(..p.., 3)), method = "wilcox.test") +
  facet_wrap(~ Hugo_Symbol, scales = "free_y") +
  theme_classic(base_size = 12)
vaf.dmg 
ggsave(vaf.dmg, file = "10. DMG VAF.pdf", width = 8, height = 4.5, units = "in") 

# ── Copy number segmented analysis ──────────────────────────────────────────── 

seg_data <- read.delim("lgg_tcga_pan_can_atlas_2018_segments.seg", header = TRUE, stringsAsFactors = FALSE) 
seg_data$ID <- gsub("-\\d+$", "", seg_data$ID) 

seg_data <- seg_data %>%
  rename(chromosome = chrom, start = loc.start, end = loc.end, probes = num.mark, segmean = seg.mean, sample = ID) %>% 
  select(chromosome, start, end, probes, segmean, sample) 

seg_plot <- seg_data %>% mutate(category_Hpx = ifelse(sample %in% sh_samp, "Severe_hypoxia", ifelse(sample %in% mh_samp, "Mild_hypoxia", "None")))  
seg_plot <- seg_plot %>% filter(category_Hpx != "None") 
seg_plot %>% count(category_Hpx)

# Basic violin plots 

vln.hpx <- ggplot(seg_plot, aes(x = category_Hpx, y = segmean, fill = category_Hpx)) + 
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8, color = NA) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black", fatten = 2, size = 0.5) +
  facet_wrap(~ chromosome, scales = "free_y") + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", label.y = mean(seg_plot$segmean)) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "Segment mean (CNA)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")
vln.hpx
ggsave(vln.hpx, file = "11. CNA violin plot.pdf", width = 6, height = 8, units = "in")

# Frequency plots 

# seg_data$start <- seg_data$start + 1 

seg_sh <- seg_data %>% filter(sample %in% sh_samp)
seg_mh <- seg_data %>% filter(sample %in% mh_samp)

set.seed(123) 

f_sh <- cnFreq(seg_sh, genome = "hg19", CN_low_cutoff = 0, CN_high_cutoff = 0.1, plot_title = "SH whole") 
f_mh <- cnFreq(seg_mh, genome = "hg19", CN_low_cutoff = 0, CN_high_cutoff = 0.1, plot_title = "MH whole") 

f_all <- f_sh / f_mh 
f_all
ggsave(f_all, file = "12. CNA frequency.pdf", width = 10, height = 4, units = "in") 

# ── Putative copy number analysis ───────────────────────────────────────────── 

# Putative copy-number from GISTIC 2.0. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.

cna_data <- read.table(file = "data_cna.txt", sep = "\t", header = TRUE) 
colnames(cna_data) <- gsub("\\.", "-", colnames(cna_data))
colnames(cna_data) <- gsub("-\\d+$", "", colnames(cna_data)) 

common_samples <- intersect(colnames(cna_data)[-(1:2)], meta_df$patient.orig_ids)
cna_subs <- cna_data[, c("Hugo_Symbol", "Entrez_Gene_Id", common_samples)]
meta_sub <- meta_df[meta_df$patient.orig_ids %in% common_samples, ] 

results <- data.frame(Gene = cna_subs$Hugo_Symbol, p.value = NA, mean_diff = NA)

for (i in 1:nrow(cna_subs)) {
  gene_vals <- as.numeric(cna_subs[i, common_samples])
  group1 <- gene_vals[meta_sub$hypoxia_class == "Mild_hypoxia"]
  group2 <- gene_vals[meta_sub$hypoxia_class == "Severe_hypoxia"] 
  
  if (length(group1) > 2 && length(group2) > 2) {
    test <- wilcox.test(group1, group2)
    results$p.value[i] <- test$p.value
    results$mean_diff[i] <- mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)
  }
}


results$adj.p.value <- p.adjust(results$p.value, method = "fdr")
top_genes <- results[order(results$adj.p.value), ]

head(top_genes, 20)

ggplot(top_genes, aes(x = mean_diff, y = -log10(adj.p.value))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Mean CNA difference (Severe - Mild hypoxia)", y = "-log10(FDR-adjusted p-value)", title = "CNA differences between hypoxia classes")



cna_data <- cna_data %>% select(!Entrez_Gene_Id) %>% column_to_rownames(var = "Hugo_Symbol") %>% t() %>% as.data.frame() 


cna_data <- cna_data %>% rownames_to_column(var = "sample")

cna_plot <- cna_data %>% mutate(category_Hpx = ifelse(sample %in% sh_samp, "Severe_hypoxia", ifelse(sample %in% mh_samp, "Mild_hypoxia", "None")))  
cna_plot <- cna_plot %>% filter(category_Hpx != "None") 

cna_plot %>% count(category_Hpx)

# Box plots 

box.cna_hpx <- ggplot(cna_plot, aes(x = category_Hpx, y = EGFR, fill = category_Hpx)) + 
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  stat_compare_means(aes(label = "p.format"), method = "t.test") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = "EGFR copy number association") 
box.cna_hpx 

ggsave(box.cna_hpx, file = "03. Putative CNA boxplots.pdf", width = 8, height = 4, units = "in")


# Bar plots 

bar.cna_hpx <- ggplot(cna_plot, aes(x = category_Hpx, y = EGFR, fill = category_Hpx)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "EGFR copy number association")
bar.cna_hpx 

ggsave(bar.cna_hpx, file = "03. Putative CNA barplots.pdf", width = 8, height = 4, units = "in")


# Pie charts 

# Value and facet columns treated as characters 

plot_pie <- function(data, group_col = "category_Hpx", value_col = "cna_status", facet_col = NULL) {
  data[[value_col]] <- as.character(data[[value_col]]) 
  if(!is.null(facet_col)) {
    data[[facet_col]] <- as.character(data[[facet_col]])
  }
  
  df_sum <- data %>%
    group_by(.data[[group_col]], .data[[value_col]], !!!rlang::syms(facet_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[group_col]], !!!rlang::syms(facet_col)) %>%
    mutate(perc = n / sum(n) * 100, perc_label = paste0(round(perc, 1), "%"), ypos = cumsum(perc) - perc/2)
  
  p_pies <- ggplot(df_sum, aes(x = "", y = perc, fill = .data[[value_col]])) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(fill = value_col, title = paste("Distribution of", value_col, "by", group_col)) + 
    geom_text(aes(y = ypos, label = perc_label), color = "black", size = 3) 
  
  if (!is.null(facet_col)) {
    p_pies <- p_pies + facet_wrap(as.formula(paste("~", facet_col)))
  } else {
    p_pies <- p_pies + facet_wrap(as.formula(paste("~", group_col)))
  } 
  return(p_pies)
} 

pie.cna_hpx <- plot_pie(cna_plot, group_col = "category_Hpx", value_col = "EGFR")

pie.cna_com <- pie.cna_hpx
pie.cna_com 
ggsave(pie.cna_com, file = "03. Putative CNA piecharts.pdf", width = 8, height = 4, units = "in")

# ── Putative arm-level copy number analysis ─────────────────────────────────── 

# Putative arm-level copy-number from GISTIC 2.0 

cna_arms <- read.table(file = "data_armlevel_cna.txt", sep = "\t", header = TRUE) 

colnames(cna_arms) <- gsub("\\.", "-", colnames(cna_arms))
colnames(cna_arms) <- gsub("-\\d+$", "", colnames(cna_arms)) 

cna.arm_long <- cna_arms %>% pivot_longer(cols = starts_with("TCGA"), names_to = "sample", values_to = "cna_status")

cna.arm_long <- cna.arm_long %>% mutate(category_Hpx = ifelse(sample %in% sh_samp, "Severe_hypoxia", ifelse(sample %in% mh_samp, "Mild_hypoxia", "None")))  
cna.arm_long <- cna.arm_long %>% filter(category_Hpx != "None") %>% na.omit() 

colnames(cna.arm_long)
str(cna.arm_long)

pie1.arm_hpx <- plot_pie(cna.arm_long, group_col = "category_Hpx", value_col = "cna_status")

pie2.arm_hpx <- plot_pie(cna.arm_long, group_col = "category_Hpx", value_col = "ENTITY_STABLE_ID")

pie.arm_com <- pie1.arm_hpx + pie2.arm_hpx 
pie.arm_com 
ggsave(pie.arm_com, file = "Putative arm-level CNA piecharts.pdf", width = 12, height = 8, units = "in")

cna_summary <- cna.arm_long %>%
  group_by(NAME, category_Hpx, cna_status) %>%
  summarise(count = n(), .groups = "drop")

cna_summary <- cna_summary %>%
  group_by(NAME, category_Hpx) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()

ggplot(cna_summary, aes(x = "", y = percent, fill = cna_status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_grid(NAME ~ category_Hpx, scales = "free_y") +
  scale_fill_manual(values = c("Loss" = "forestgreen", "Gain" = "red", "Unchanged" = "grey80")) +
  theme_void() +
  theme(strip.text = element_text(size = 10, face = "bold"), legend.position = "bottom") 

# ── Metadata plots for supplementary ────────────────────────────────────────── 

metadata <- read.csv("LGG_HP_metaclus.csv") 
meta.res <- metadata %>% rename(PATIENT_ID = "X")
meta.res <- meta.res %>% mutate(category_mut = ifelse(PATIENT_ID %in% mut_sh_cha, "Mut.SH_C1", ifelse(PATIENT_ID %in% wld_sh_cha, "Wild.SH_C1", ifelse(PATIENT_ID %in% mh_cha, "MH_C2", "None")))) 

meta.res <- meta.res %>% select(PATIENT_ID, SUBTYPE, CANCER_TYPE_DETAILED, GRADE, AGE, Cluster_groups, category_mut)
colnames(meta.res)

meta.res <- meta.res %>% filter(category_mut != "None")

# Color palettes  

gradient_colors <- list(AGE = c("white", "lightgreen", "darkgreen"))

category_colors <- list(SUBTYPE = c("LGG_IDHmut-codel" = "tomato", "LGG_IDHmut-non-codel" = "gold", "LGG_IDHwt" = "deepskyblue", "None" = "grey70"),
                        CANCER_TYPE_DETAILED = c("Astrocytoma" = "forestgreen", "Low-Grade Glioma (NOS)" = "orange", "Oligoastrocytoma" = "red", "Oligodendroglioma" = "purple"),
                        GRADE = c("G2" = "darkblue", "G3" = "darkred", "NA" = "grey70"), 
                        Cluster_groups = c("C1" = "green", "C2" = "#00BFC4"), 
                        category_mut = c("MH_C2" = "green", "None" = "grey70", "Mut.SH_C1" = "#FF2DF1", "Wild.SH_C1" = "#FA812F")) 

# Functions for plots 

plot_numeric <- function(varname, colors) {
  ggplot(meta.res, aes(x = PATIENT_ID, y = 1, fill = .data[[varname]])) +
    geom_tile() +
    scale_fill_gradientn(colors = colors) +
    theme_void() +
    theme(legend.position = "bottom", plot.margin = margin(0, 5, 0, 5)) +
    labs(fill = varname)
}

plot_categorical <- function(varname, colors) {
  ggplot(meta.res, aes(x = PATIENT_ID, y = 1, fill = .data[[varname]])) +
    geom_tile() +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(legend.position = "bottom", plot.margin = margin(0, 5, 0, 5)) +
    labs(fill = varname)
}

meta.res <- meta.res %>%
  mutate(category_mut = factor(category_mut, levels = c("Mut.SH_C1", "Wild.SH_C1", "MH_C2", "None"))) %>%
  arrange(category_mut, SUBTYPE, GRADE) %>% 
  mutate(PATIENT_ID = factor(PATIENT_ID, levels = unique(PATIENT_ID)))

# Plot bar plot as the base 

p_meta.lgg <- ggplot(meta.res, aes(x = PATIENT_ID, y = AGE, fill = Cluster_groups)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = quantile(meta.res$AGE, 0.5, na.rm = TRUE), linetype = "dashed", color = "red", size = 0.8) +
  facet_grid(. ~ category_mut, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(5, 5, 2, 5)) 

# Combine the annotation plots

num_vars <- names(gradient_colors)
cat_vars <- names(category_colors)

num_plots <- mapply(plot_numeric, varname = num_vars, colors = gradient_colors, SIMPLIFY = FALSE) 
cat_plots <- mapply(plot_categorical, varname = cat_vars, colors = category_colors, SIMPLIFY = FALSE) 
all_plots <- list(p_meta.lgg) %>% append(cat_plots) %>% append(num_plots)

rel_height <- c(4, rep(1, length(cat_plots) + length(num_plots)))
p_combined <- wrap_plots(all_plots, ncol = 1, heights = rel_height)

f.plot <- print(p_combined)
ggsave(f.plot, file = "05. Metadata annotation plots.pdf", width = 8, height = 6, units = "in")

# Pie charts 

pie.mut <- plot_pie(meta.res, group_col = "Cluster_groups", value_col = "category_mut")
pie.grd <- plot_pie(meta.res, group_col = "Cluster_groups", value_col = "GRADE")
pie.sub <- plot_pie(meta.res, group_col = "Cluster_groups", value_col = "SUBTYPE")
pie.typ <- plot_pie(meta.res, group_col = "Cluster_groups", value_col = "CANCER_TYPE_DETAILED")

pie.meta.supp <- pie.mut + pie.grd + pie.sub + pie.typ 
pie.meta.supp 
ggsave(pie.meta.supp, file = "06. Metadata percentage pie.pdf", width = 12, height = 8, units = "in")

