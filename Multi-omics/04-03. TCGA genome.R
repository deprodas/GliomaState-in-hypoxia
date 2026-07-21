setwd("C:\\01. Research\\01. mIDH (hypoxia)\\02. Analysis\\01. Analysis (cBioportal)\\01-05. Mutation")

# Author : Depro Das, Medical Center - University of Freiburg

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
library(survival)
library(survminer)

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

# Frequency plots (https://genviz.org/module-03-genvisr/0003/05/01/cnFreq_GenVisR/) 

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

cna_subs <- cna_subs %>% 
  filter(Hugo_Symbol %in% c("EGFR")) %>% 
  select(!Entrez_Gene_Id) %>% column_to_rownames(var = "Hugo_Symbol") %>% t() %>% as.data.frame() 

cna_subs <- cna_subs %>%
  mutate(EGFR_CNA_status = ifelse(EGFR == -2, "homozygous_deletion", 
                           ifelse(EGFR == -1, "hemizygous_deletion", 
                           ifelse(EGFR == 0, "neutral_no_change", 
                           ifelse(EGFR == 1, "gain", 
                           ifelse(EGFR == 2, "high_level_amplification", NA))))))

cna_subs <- cna_subs %>% dplyr::rename(EGFR_GISTIC = "EGFR")

all(rownames(cna_subs) %in% rownames(meta_sub))
all(rownames(meta_sub) %in% rownames(cna_subs))

cna_subs <- cbind(meta_sub, cna_subs)
cna_subs %>% group_by(hypoxia_class) %>% count(EGFR_GISTIC)

# Pie charts 

# Value and facet columns treated as characters 

plot_pie <- function(data, group_col = "hypoxia_groups", value_col = "AGE", facet_col = NULL) {
  data[[value_col]] <- as.character(data[[value_col]]) 
  if(!is.null(facet_col)) {
    data[[facet_col]] <- as.character(data[[facet_col]])
  } 
  df_sum <- data %>%
    group_by(.data[[group_col]], .data[[value_col]], !!!rlang::syms(facet_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[group_col]], !!!rlang::syms(facet_col)) %>%
    mutate(perc = n / sum(n) * 100, perc_label = paste0("N=", n, "\n(", round(perc, 1), "%)"), ypos = cumsum(perc) - perc/2)
  
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

pie.cna_egfr <- plot_pie(cna_subs, group_col = "hypoxia_class", value_col = "EGFR_CNA_status")
pie.cna_egfr
ggsave(pie.cna_egfr, file = "13. EGFR CNA metapie.pdf", width = 14, height = 10, units = "in")

# Add pathway enrichment result 

res_path <- read.csv("Result ssGSEA (MPs-QAD-Progeny).csv", row.names = 1)
res_path <- res_path %>% select(EGFR_50, EGFR_all)
cna_subs <- cbind(cna_subs, res_path)

# Categories with too few observations were excluded to support reliable cox modeling 

cna_subs %>% dplyr::count(EGFR_CNA_status)
cna_filt <- cna_subs %>% dplyr::filter(EGFR_CNA_status %in% c('gain', 'neutral_no_change'))

# Cox total 

surv.obj_cna <- Surv(time = cna_filt$OS_MONTHS, event = cna_filt$OS_STATUS)

cox_cna <- coxph(surv.obj_cna ~ hypoxia_class + EGFR_CNA_status + EGFR_all, data = cna_filt)
for_cna <- ggforest(cox_cna, data = cna_filt) 
for_cna 
ggsave(file = "14. Cox regression EGFR CNA (Total).pdf", plot = for_cna, width = 6, height = 4, units = "in") 

# Cox separate 

sh_m.cna <- cna_filt %>% filter(hypoxia_class == "Severe_hypoxia") 
mh_m.cna <- cna_filt %>% filter(hypoxia_class == "Mild_hypoxia") 

sh_m.cna$EGFR_CNA_status <- relevel(factor(sh_m.cna$EGFR_CNA_status), ref = "neutral_no_change")
mh_m.cna$EGFR_CNA_status <- relevel(factor(mh_m.cna$EGFR_CNA_status), ref = "neutral_no_change")

surv.obj_sh.cna <- Surv(time = sh_m.cna$OS_MONTHS, event = sh_m.cna$OS_STATUS)
surv.obj_mh.cna <- Surv(time = mh_m.cna$OS_MONTHS, event = mh_m.cna$OS_STATUS)

cox_sh.cna <- coxph(surv.obj_sh.cna ~ EGFR_CNA_status + EGFR_50 + EGFR_all, data = sh_m.cna)
for_sh.cna <- ggforest(cox_sh.cna, data = sh_m.cna) 
for_sh.cna 
ggsave(file = "14. Cox regression EGFR CNA (SH).pdf", plot = for_sh.cna, width = 6, height = 4, units = "in") 

cox_mh.cna <- coxph(surv.obj_mh.cna ~ EGFR_CNA_status + EGFR_50 + EGFR_all, data = mh_m.cna)
for_mh.cna <- ggforest(cox_mh.cna, data = mh_m.cna) 
for_mh.cna 
ggsave(file = "14. Cox regression EGFR CNA (MH).pdf", plot = for_mh.cna, width = 6, height = 4, units = "in") 

# ── Putative arm-level copy number analysis ─────────────────────────────────── 

# Putative arm-level copy-number from GISTIC 2.0 

cna_arms <- read.table(file = "data_armlevel_cna.txt", sep = "\t", header = TRUE) 

colnames(cna_arms) <- gsub("\\.", "-", colnames(cna_arms))
colnames(cna_arms) <- gsub("-\\d+$", "", colnames(cna_arms)) 

cna.arm_long <- cna_arms %>% pivot_longer(cols = starts_with("TCGA"), names_to = "sample", values_to = "cna_status")

cna.arm_long <- cna.arm_long %>% mutate(category_Hpx = ifelse(sample %in% sh_samp, "Severe_hypoxia", ifelse(sample %in% mh_samp, "Mild_hypoxia", "None")))  
cna.arm_long <- cna.arm_long %>% filter(category_Hpx != "None") %>% na.omit() 

cna_summmary <- cna.arm_long %>%
  group_by(NAME, category_Hpx, cna_status) %>%
  summarise(count = n(), .groups = "drop")

cna_summmary <- cna_summmary %>%
  group_by(NAME, category_Hpx) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()

pie.arm_com <- ggplot(cna_summmary, aes(x = "", y = percent, fill = cna_status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_grid(NAME ~ category_Hpx, scales = "free_y") +
  scale_fill_manual(values = c("Loss" = "forestgreen", "Gain" = "red", "Unchanged" = "grey80")) +
  theme_void() +
  theme(strip.text = element_text(size = 10, face = "bold"), legend.position = "bottom") 
pie.arm_com
ggsave(pie.arm_com, file = "15. Putative arm-level CNA piecharts.pdf", width = 10, height = 60, units = "in", limitsize = FALSE)
