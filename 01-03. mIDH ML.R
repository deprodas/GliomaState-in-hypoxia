# Author : Depro Das, Department of Neurosurgery, University Hospital Freiburg, Freiburg, Germany 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(openxlsx) 
library(stringr)
library(ggplot2)
library(randomForest)
library(xgboost)
library(e1071) 
library(glmnet) 
library(pROC)
library(ggrepel)

# ── Prepare data ────────────────────────────────────────────────────────────── 

# Read limma results and metamodule genesets 

limmadata <- read.csv("Results limma (SH-MH).csv") %>%
  na.omit() %>%
  rename(genes = "X") %>%
  mutate(nlog10 = -log10(P.Value))

meta_gene <- read.xlsx("IDHm metamodule genes.xlsx") %>% dplyr::select(OD_OPC_like, OD_Astro_like, OD_Cycling, OD_RA) 
meta_labs <- meta_gene %>%
  pivot_longer(cols = everything(), names_to = "Subtype", values_to = "Gene") %>%
  filter(!is.na(Gene)) %>%
  distinct()

# Binary label: 1 = metamodule gene, 0 = non-marker

limma_bin <- limmadata %>% mutate(Marker = ifelse(genes %in% meta_labs$Gene, 1, 0))

# ── Feature matrix ──────────────────────────────────────────────────────────── 

features <- limma_bin %>% 
  transmute(logFC = logFC, AveExpr = AveExpr) %>%
  scale() %>% 
  as.data.frame()

gene_lab <- limma_bin$Marker

# ── Split and train ─────────────────────────────────────────────────────────── 

# Train-test split 

set.seed(123)
train_idx <- sample(seq_len(nrow(features)), size = 0.7 * nrow(features))

X_train <- features[train_idx, ]
y_train <- gene_lab[train_idx]
X_tests <- features[-train_idx, ]
y_tests <- gene_lab[-train_idx]

# Train models 

# Random forest

rnf_model <- randomForest(x = X_train, y = as.factor(y_train), ntree = 500, importance = TRUE)

# Xgboost 

dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
dtests <- xgb.DMatrix(data = as.matrix(X_tests), label = y_tests)

xgb_model <- xgboost(data = dtrain, nrounds = 200, objective = "binary:logistic", eval_metric = "auc", verbose = 0)

# SVM

svm_model <- svm(x = X_train, y = as.factor(y_train), probability = TRUE, kernel = "radial")

# Elastic net (glmnet)

ent_model <- cv.glmnet(as.matrix(X_train), y_train, family = "binomial", alpha = 0.5)

# ── Predict probabilities ─────────────────────────────────────────────────────  

all_features <- features

rnf_all_pred <- predict(rnf_model, all_features, type = "prob")[,2]
svm_all_pred <- attr(predict(svm_model, all_features, probability = TRUE), "probabilities")[,2]
ent_all_pred <- predict(ent_model, as.matrix(all_features), type = "response")[,1]

dens_all_mtx <- xgb.DMatrix(data = as.matrix(all_features))
xgb_all_pred <- predict(xgb_model, dens_all_mtx)

# Combine probabilities 

cons.df_ml <- limma_bin %>%
  mutate(RnF_prob = rnf_all_pred, XGB_prob = xgb_all_pred, SVM_prob = svm_all_pred, Ent_prob = ent_all_pred, consensus_prob = (RnF_prob + XGB_prob + SVM_prob + Ent_prob) / 4) %>%
  arrange(desc(consensus_prob))

# Composite scoring

cons.df_ml <- cons.df_ml %>%
  mutate(probability_SD = apply(cbind(RnF_prob, XGB_prob, SVM_prob, Ent_prob), 1, sd),
         adjusted_score = consensus_prob * (1 - probability_SD),
         composite_score = consensus_prob * 0.5 + (nlog10 / max(nlog10, na.rm = TRUE)) * 0.3 + (abs(logFC) / max(abs(logFC), na.rm = TRUE)) * 0.2)

# Annotate 

cons.df_mk <- cons.df_ml %>% 
  left_join(meta_labs, by = c("genes" = "Gene")) %>% 
  mutate(Subtype = ifelse(is.na(Subtype), "Other", Subtype))

write.csv(cons.df_mk, "Result ML-predicted markers.csv", row.names = FALSE)

# ── Plot genes ──────────────────────────────────────────────────────────────── 

cons.df_mk <- cons.df_mk %>% filter(consensus_prob > 0.01)

# Order genes 

cons.df_mk <- cons.df_mk %>% 
  arrange(desc(consensus_prob)) %>%
  mutate(genes = factor(genes, levels = unique(genes)))

p.ml_bar <- ggplot(cons.df_mk, aes(x = consensus_prob, y = genes, color = Subtype)) +
  geom_point(size = 1, alpha = 0.7) +
  coord_flip() +
  theme_bw(base_size = 10) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 
p.ml_bar 
ggsave(filename = "01. Marker probability.pdf", plot = p.ml_bar, width = 10, height = 6, units = "in")

# Zoom-in 

label_genes <- cons.df_mk %>%
  filter(consensus_prob >= 0.25 & consensus_prob <= 0.5) %>%
  arrange(desc(consensus_prob)) %>%
  slice_head(n = 10)

p.zm_bar <- ggplot(cons.df_mk, aes(x = consensus_prob, y = genes, color = Subtype)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text_repel(data = label_genes, aes(label = genes), size = 2.5, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2, segment.size = 0.2) +
  coord_flip() +
  scale_y_discrete(limits = rev) + 
  xlim(0.25, 0.5) +
  theme_bw(base_size = 10) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
p.zm_bar
ggsave(filename = "02. Marker probability zoom.pdf", plot = p.zm_bar, width = 10, height = 6, units = "in")


# ── ROC curve ───────────────────────────────────────────────────────────────── 

roc_rnf <- roc(y_tests, predict(rnf_model, X_tests, type = "prob")[,2])
roc_xgb <- roc(y_tests, predict(xgb_model, dtests))
roc_svm <- roc(y_tests, attr(predict(svm_model, X_tests, probability = TRUE), "probabilities")[,2])
roc_ent <- roc(y_tests, predict(ent_model, as.matrix(X_tests), type = "response")[,1])

roc_to_df <- function(roc_obj, model_name) {
  data.frame(TPR = rev(roc_obj$sensitivities), FPR = rev(1 - roc_obj$specificities), Model = model_name)
}

roc_df <- bind_rows(roc_to_df(roc_rnf, "Random Forest"),
                    roc_to_df(roc_xgb, "XGBoost"),
                    roc_to_df(roc_svm, "SVM"),
                    roc_to_df(roc_ent, "Elastic Net"))

auc_df <- data.frame(Model = c("Random Forest", "XGBoost", "SVM", "Elastic Net"), 
                     AUC = c(auc(roc_rnf), auc(roc_xgb), auc(roc_svm), auc(roc_ent)))

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_text(data = auc_df, aes(x = 0.6, y = seq(0.4, 0.1, length.out = nrow(auc_df)), label = paste0(Model, " AUC = ", round(AUC, 3))), inherit.aes = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom")
p_roc 
ggsave(filename = "03. ROC curve.pdf", plot = p_roc, width = 6, height = 6, units = "in")
