############################################################
# 1. Load libraries
############################################################

library(celldex)
library(dplyr)
library(ggplot2)
library(anndata)
library(grid)
library(devtools)
library(SeuratData)
library(Seurat)
library(homologene)
library(Seurat)
library(reticulate)
library(anndata)
library(dplyr)
library(stringr)
library(tidyverse)
library(fastplyr)
library(Matrix)
library(stringr)
library(biomaRt)
library(anndata)
library(tibble)
library(dplyr)
library(edgeR)


############################################################
# 2. Load custom CytoSig helper functions
############################################################

cyto <- new.env()
source("lib/R-functions.R", local = cyto)

attach(cyto)


############################################################
# 3. Load pseudobulk AnnData object
############################################################

# pseudobulk_data <- read_h5ad(
#   "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/pdata.h5ad"
# )

pseudobulk_data <- read_h5ad(
  "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/single_cell_normal/pdata_cytosig_normal.h5ad"
)


############################################################
# 4. Convert counts and normalize with edgeR
############################################################

# Convert AnnData matrix to genes x samples
data_ex <- t(as.matrix(pseudobulk_data$X))

# Create DGE object
dge <- DGEList(counts = data_ex)

# Normalize for library size using TMM
dge <- calcNormFactors(dge)

# Compute log2 CPM values
logCPM <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)


############################################################
# 5. Convert mouse genes to human orthologs
############################################################

# homologene() queries homologous genes between species.

# The resulting table contains multiple columns.
# Keep only:
#   - mouse gene symbol
#   - human ortholog symbol
# ------------------------------------------------------------------
# Step 3: Keep only genes present in the expression matrix
# ------------------------------------------------------------------

# Remove mappings where the mouse gene is not present
# in the expression matrix
# Remove mappings with missing human orthologs


results <- list()

# logCPM is genes x samples
# Make sure rownames are mouse genes
rownames(logCPM) <- rownames(data_ex)
colnames(logCPM) <- colnames(data_ex)

orth <- homologene(
  rownames(logCPM),
  inTax = 10090,   # mouse
  outTax = 9606    # human
)

orth <- orth[, c("10090", "9606")]
colnames(orth) <- c("mouse_gene", "human_gene")

orth <- orth[orth$mouse_gene %in% rownames(logCPM), ]
orth <- orth[orth$human_gene != "", ]

# Keep only clean 1:1 mouse-human mappings
orth_unique <- orth[
  !duplicated(orth$human_gene) &
    !duplicated(orth$human_gene, fromLast = TRUE) &
    !duplicated(orth$mouse_gene) &
    !duplicated(orth$mouse_gene, fromLast = TRUE),
]

data_ex_human <- logCPM[orth_unique$mouse_gene, , drop = FALSE]
rownames(data_ex_human) <- orth_unique$human_gene

dim(data_ex_human)
head(rownames(data_ex_human))




############################################################
# 6. Run CytoSig on pseudobulk logCPM data
############################################################

results$CS_pseudobulk <- runCytoSig(
  data_ex_human,
  expand_signature = 0,
  resultsdir = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/cytosig_pseudobulk_logCPM",
  condaenv = "cytosig.v0.1",
  overwrite = TRUE
)


############################################################
# 7. Extract CytoSig Z-scores and p-values
############################################################

cs_mat <- as.matrix(results[["CS_pseudobulk"]][["Zscore"]])
p_mat  <- as.matrix(results[["CS_pseudobulk"]][["Pvalue"]])


############################################################
# 8. Create sample condition annotation
############################################################

ann_col <- data.frame(
  condition = case_when(
    grepl("^ctrl", colnames(cs_mat)) ~ "ctrl",
    grepl("^GF_noICI", colnames(cs_mat)) ~ "GF_noici",
    grepl("^GF", colnames(cs_mat)) ~ "GF",
    grepl("^effector|^11mix", colnames(cs_mat)) ~ "effector",
    TRUE ~ "unknown"
  )
)

rownames(ann_col) <- colnames(cs_mat)

table(ann_col$condition, useNA = "ifany")
table(ann_col$condition, useNA = "ifany")

all(rownames(ann_col) == colnames(cs_mat))

rownames(ann_col) <- colnames(cs_mat)


############################################################
# 9. Create significance matrix from p-values
############################################################

# Make sure dimensions/names match
p_mat <- p_mat[rownames(cs_mat), colnames(cs_mat)]

sig_mat <- ifelse(
  p_mat < 0.001, "***",
  ifelse(
    p_mat < 0.01, "**",
    ifelse(p_mat < 0.05, "*", "")
  )
)


############################################################
# 10. Reorder samples manually
############################################################

# col_order <- c("GF_noici1_all","GF_noici2_all","GF1_all","GF2_all","ctrl1_all","ctrl2_all", "effector1_all", "effector2_all")

col_order <- c(

  "GF1_all",
  "GF2_all",
  "ctrl1_all",
  "ctrl2_all",
  "effector1_all",
  "effector2_all"
)

# Reorder matrix and significance annotation
cs_mat  <- cs_mat[, col_order]
sig_mat <- sig_mat[, col_order]

# Annotation rows must match column names
ann_col <- ann_col[col_order, , drop = FALSE]


############################################################
# 11. Plot CytoSig heatmap for individual samples
############################################################

q <- pheatmap::pheatmap(
  cs_mat,
  scale = "row",
  annotation_col = ann_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize = 14,
  border_color = NA,
  
  display_numbers = sig_mat,
  number_color = "black",
  fontsize_number = 14
)

# ggsave(
#   "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/figures/cytosig/cytosig_heatmap_cytokine_activity_group_sig_zscore.svg",
#   plot = q,
#   width = 8,
#   height = 14,
#   dpi = 300,
#   bg = "white"
# )


############################################################
# 12. Define biological replicate groups
############################################################

group_map <- c(
  "GF1_all"       = "GF",
  "GF2_all"       = "GF",
  "ctrl1_all"     = "ctrl",
  "ctrl2_all"     = "ctrl",
  "effector1_all" = "effector",
  "effector2_all" = "effector"
)


############################################################
# 13. Average only direction-consistent significant signals
############################################################

# Convert significance matrix to logical TRUE/FALSE
sig_logical <- sig_mat != ""

# Unique biological groups
groups <- unique(group_map)

# Empty result matrix for averaged CytoSig scores
cs_avg <- matrix(
  NA,
  nrow = nrow(cs_mat),
  ncol = length(groups),
  dimnames = list(rownames(cs_mat), groups)
)

# Empty result matrix for significance labels
sig_avg <- matrix(
  "",
  nrow = nrow(cs_mat),
  ncol = length(groups),
  dimnames = list(rownames(cs_mat), groups)
)

for(g in groups){
  
  cols <- names(group_map[group_map == g])
  
  # Extract the two replicates for this group
  x1 <- cs_mat[, cols[1]]
  x2 <- cs_mat[, cols[2]]
  
  # Extract significance for the two replicates
  s1 <- sig_logical[, cols[1]]
  s2 <- sig_logical[, cols[2]]
  
  # Keep cytokines where both replicates have the same direction
  same_dir <- sign(x1) == sign(x2)
  
  # Keep cytokines significant in at least one replicate
  sig_keep <- s1 | s2
  
  # Alternative stricter option:
  # sig_keep <- s1 & s2
  
  # Final filter: same direction and significant
  keep <- same_dir & sig_keep
  
  # Average the two replicates
  avg <- rowMeans(cbind(x1, x2), na.rm = TRUE)
  
  # Set inconsistent/non-significant cytokines to NA
  avg[!keep] <- NA
  
  # Store averaged scores
  cs_avg[, g] <- avg
  
  # Store significance annotation
  sig_avg[keep, g] <- "*"
}

cs_avg <- as.data.frame(cs_avg)
sig_avg <- as.data.frame(sig_avg)


############################################################
# 14. Prepare grouped heatmap matrix
############################################################

library(pheatmap)

cs_plot <- as.matrix(cs_avg)

# Remove cytokines with no retained signal in any group
cs_plot <- cs_plot[rowSums(!is.na(cs_plot)) > 0, ]

# Replace remaining NA with 0 for clustering/plotting
cs_plot[is.na(cs_plot)] <- 0


############################################################
# 15. Create grouped heatmap annotation
############################################################

annotation_col <- data.frame(
  Group = factor(
    colnames(cs_avg),
    levels = c("GF", "ctrl", "effector")
  )
)

rownames(annotation_col) <- colnames(cs_plot)


############################################################
# 16. Plot grouped CytoSig heatmap
############################################################

q <- pheatmap(
  mat = cs_plot,
  
  scale = "none",
  
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  
  annotation_col = annotation_col,
  
  color = colorRampPalette(c("blue", "white", "red"))(100),
  
  display_numbers = sig_avg[
    rownames(cs_plot),
    colnames(cs_plot)
  ],
  
  number_color = "black",
  fontsize_number = 10,
  
  fontsize = 14,
  fontsize_row = 12,
  fontsize_col = 14,
  
  border_color = NA,
  
  main = "CytoSig grouped pseudobulk"
)