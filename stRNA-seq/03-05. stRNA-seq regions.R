# Author : Depro Das, Department of Neurosurgery, University Medical Center Freiburg, Freiburg, Germany

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(Seurat)
library(dplyr)
library(igraph)
library(sf) 
library(lwgeom) 
library(ggplot2)

# ── Data input & prepare ────────────────────────────────────────────────────── 

spt_data <- readRDS("BWH23.RDS")
coord_df <- GetTissueCoordinates(spt_data)
coord_df <- data.frame(barcode = rownames(coord_df),
                       x = coord_df[, 1], 
                       y = coord_df[, 2]) 
metadata <- spt_data@meta.data
metadata$barcode <- rownames(metadata)
spt_meta <- left_join(coord_df, metadata, by = "barcode")

# Extract scaling factors

spot_pix <- spt_data@images$slice1@scale.factors$spot
spot_rad <- spot_pix / 2
um_per_pixel <- 55 / spot_pix

# Neighborhood graph (dynamically estimate the true nearest-neighbor distance in pixel space)

hpcat_df <- spt_meta %>% filter(Hypoxia_class == "True_hypoxia")
xy_value <- as.matrix(hpcat_df[, c("x", "y")])
nn_pdist <- median(apply(as.matrix(dist(xy_value)), 1, function(x) min(x[x > 0])))

# Connect spots within 1.5 spot distances to capture hexagonal adjacencies

radius <- nn_pdist * 1.5
D <- as.matrix(dist(xy_value))
adj <- D < radius & D > 0 
g <- graph_from_adjacency_matrix(adj, mode = "undirected")

# Connected components & filter

cc <- components(g)
hpcat_df$cluster <- factor(cc$membership)
hpcat_df %>% dplyr::count(cluster) 

cl_sizes <- table(hpcat_df$cluster)
large_cl <- names(cl_sizes[cl_sizes >= 50])
hp_large <- hpcat_df %>% filter(cluster %in% large_cl)
hp_large$cluster <- flavor_cluster <- factor(hp_large$cluster)

# Calculate area (for loop : first half calculates areas per cluster, 2nd half for circular polygon visualization) 

unq_large_cl <- unique(hp_large$cluster)
cl_area_list <- list()
circles_list <- list() 

for (cl in unq_large_cl) {
  cluster_df <- hp_large %>% filter(cluster == cl) 
  cluster_sf <- sf::st_as_sf(cluster_df, coords = c("x", "y"))
  circles_sp <- sf::st_buffer(cluster_sf, dist = spot_rad)
  cluster_un <- sf::st_union(circles_sp)
  py.area_px <- as.numeric(sf::st_area(cluster_un))
  py.area_mm <- py.area_px * (um_per_pixel ^ 2) / 1e6
  cl_area_list[[as.character(cl)]] <- data.frame(Clust_ids = paste("Cluster", cl),
                                                 Tot_spots = nrow(cluster_df),
                                                 Areas_mm2 = round(py.area_mm, 4)) 
  
  cl_polygon <- lwgeom::st_minimum_bounding_circle(st_union(cluster_sf))
  polygon_sf <- sf::st_sf(cluster = factor(cl), geometry = cl_polygon)
  circles_list[[as.character(cl)]] <- polygon_sf
}

# Combine results

cl_summary <- do.call(rbind, cl_area_list)
rownames(cl_summary) <- NULL
print(cl_summary) 
all_circle <- do.call(rbind, circles_list)

# Plot

p.spt_meta <- sf::st_as_sf(spt_meta, coords = c("x", "y"), remove = FALSE)
p.hp_large <- sf::st_as_sf(hp_large, coords = c("x", "y"), remove = FALSE)

pt_spatial <- ggplot(data = p.spt_meta) +
  geom_sf(color = "#E5E7E9", size = 0.5, alpha = 0.4, inherit.aes = FALSE) +
  geom_sf(data = all_circle, aes(fill = factor(cluster), col = factor(cluster)), linewidth = 0.8, alpha = 0.15, inherit.aes = FALSE) + 
  geom_sf(data = p.hp_large, aes(col = factor(cluster)), size = 1, alpha = 0.9, inherit.aes = FALSE) + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") + 
  coord_sf(ylim = c(max(spt_meta$y), min(spt_meta$y))) + 
  theme_void() + 
  ggtitle(paste0("Sample: ", unique(spt_meta$sample_ID), "_niches-", length(unq_large_cl)))
pt_spatial 
ggsave(filename = paste0("1. ", unique(spt_meta$sample_ID), " niches_mm2.pdf"), plot = pt_spatial, width = 4, height = 5)
