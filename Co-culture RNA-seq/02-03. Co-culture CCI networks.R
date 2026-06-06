# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(jetset)
library(ggplot2)
library(openxlsx) 
library(dplyr)
library(tidyverse) 
library(tidyr) 
library(viridis) 
library(icellnet)
library(gridExtra)
library(ComplexHeatmap) 
library(circlize) 
library(ggplotify) 
library(ggraph) 
library(ggsankey) 
library(RColorBrewer) 

# ── Load and select database ──────────────────────────────────────────────────

db = as.data.frame(read.csv(curl::curl(url = "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
db.name.couple = name.lr.couple(db, type = "Family")
db.name.subfam = name.lr.couple(db, type = "Subfamily")

head(db.name.couple)
head(db.name.subfam)

# ── Load bulk RNA-seq data ──────────────────────────────────────────────────── 

data = as.data.frame(read.csv("star_count_symbol.csv", sep=",", header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = "")) 
rownames(data) = data$Gene 

# ── Data manipulation and selection ───────────────────────────────────────────

target = openxlsx::read.xlsx("Sample_Identifier.xlsx")
target$Group_Condition <- paste(target$Group, target$Condition, sep = "_")

PC.target <- target %>% 
  dplyr::select(Code, Group, Group_Condition) %>% 
  rename(ID = "Code") %>% 
  rename(Cell_type = "Group") %>% 
  rename(Class = "Group_Condition")

rownames(PC.target) = PC.target$ID
data.scaled = gene.scaling(data = data, n = 1, db = db)

# Selection 

CC.data.selection = PC.target$ID[which(PC.target$Class %in% c("LGG_Hypoxia"))] # Should give a vector with name of replicates to consider for central cell
PC.data.selection = PC.target$ID[which(PC.target$Class %in% c("LGG_MICROGLIA_Hypoxia", "LGG_ASTROCYTE_Hypoxia", "LGG_Hypoxia", "LGG_Normoxia", "LGG_MICROGLIA_Normoxia", "LGG_ASTROCYTE_Normoxia"))] # Should give a vector with name of replicates to consider for partner cell 
my_Central_Cell_data = data.scaled[, CC.data.selection]
my_Partner_Cell_data = data.scaled[, PC.data.selection]

# ── Compute communication score ─────────────────────────────────────────────── 

score.computation.1 = icellnet.score(direction = "in", 
                                     PC.data = my_Partner_Cell_data, 
                                     CC.data = my_Central_Cell_data, 
                                     PC.target = PC.target,
                                     PC = c("LGG_MICROGLIA_Hypoxia", "LGG_ASTROCYTE_Hypoxia", "LGG_MICROGLIA_Normoxia", "LGG_ASTROCYTE_Normoxia", "LGG_Hypoxia", "LGG_Normoxia"), 
                                     CC.type = "RNAseq", 
                                     PC.type = "RNAseq",  
                                     db = db) 

score1 = as.data.frame(score.computation.1[[1]])
lr = score.computation.1[[2]] 
lr1 = lr %>% as.data.frame() %>% rownames_to_column(var = "Signaling")
write.xlsx(lr1, file = "All L-R interactions.xlsx")

# Select and visualize specific LR interactions 

lr1 <- lr1 %>% as.data.frame() %>% filter(if_any(everything(), ~ grepl("/ EGFR\\b", ., ignore.case = TRUE)))
lr1 <- na.omit(lr1) 
write.xlsx(lr1, file = "Specific L-R interactions.xlsx")

# ── Plot in different ways ────────────────────────────────────────────────────  

# Heatmap 

df_org <- read.xlsx("Specific L-R interactions.xlsx")
df_med <- df_org %>% remove_rownames() %>% column_to_rownames(var = "Signaling")

set.seed(123)
lr_pht <- Heatmap(as.matrix(df_med), 
                  col = colorRamp2(c(0, 30, 60), c("darkblue", "white","red")),
                  border = T) 
lr_pht
lr_pht <- as.ggplot(lr_pht) 
ggsave(filename = "1. Icellnet heatmap.pdf", plot = lr_pht, width = 4.5, height = 6, units = c("in"))

# Sankey 

df_sepr <- df_org %>% separate(Signaling, into = c("Ligand", "Receptor"), sep = " / ")
df_long <- df_sepr %>%
  pivot_longer(-c(Ligand, Receptor), names_to = "Condition", values_to = "Strength") %>%
  filter(Strength > 0) %>%
  mutate(Receptor_Condition = paste(Receptor, Condition, sep = "_")) 

df_sank <- df_long %>%
  select(Condition, Receptor, Ligand, Strength) %>%
  mutate(across(c(Condition, Ligand, Receptor), as.character), Strength = Strength * 10) %>% 
  make_long(Condition, Ligand, Receptor, value = Strength)

lr_sank <- ggplot(df_sank, aes(x = x, next_x = next_x, node = node, next_node = next_node, value = value, fill = node)) + 
  geom_sankey(flow.alpha = 0.7, node.color = "black", show.legend = FALSE, space = 1500) +
  geom_sankey_label(aes(label = node), size = 2, color = "black", fill = "white") + 
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank())
lr_sank 
ggsave(filename = "3. Icellnet sankey.pdf", plot = lr_sank, width = 3, height = 3, units = c("in"))

