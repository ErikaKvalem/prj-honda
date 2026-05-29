############################################################
# Load required libraries and resolve function conflicts
############################################################

library("conflicted")
library("readr")
library("dplyr")

# Prefer dplyr functions when there are package conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("first", "dplyr")


library("tibble")
library("stringr")


# Prefer base/dplyr functions for additional conflicts
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
conflict_prefer("as.factor", "base")


############################################################
# Load differential expression results
############################################################

# Load all genes from effector vs ctrl DE result table
resIHW <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWallGenes.tsv")

# Remove antisense genes
resIHW <- resIHW[resIHW$gene_type != "antisense", ]



############################################################
# Load count matrix and significant gene table
############################################################

# Load raw/count matrix
counts <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/counts_matrix.tsv")

# Load significant genes and remove antisense genes
sig_genes_df <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWsigGenes.tsv") %>%
  dplyr::filter(gene_type != "antisense")

############################################################
# Define significance thresholds
############################################################

fc_cutoff = 0.58
fdr_cutoff= 0.1

############################################################
# Load libraries for ComplexHeatmap
############################################################

library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(colorRamp2)

# Prefer base setdiff and colorRamp2 function
conflicts_prefer(base::setdiff)
conflicts_prefer(colorRamp2::colorRamp2)

############################################################
# Extract significant genes
############################################################

# Get unique significant gene names
sig_genes <- sig_genes_df$gene_name %>%
  na.omit() %>%
  unique()

############################################################
# Define gene modules for row annotation
############################################################

module_list <- list(
  "TCR repertoire / clonal signature" = c("Trav7-1", "Trav8-2", "Trav10n", "Tcrg-C3","Tcrg-V3","Trav9-4","Trbj1-5"),
  "Effector CD8 T cell module" = c("Cd7", "Cxcr3"),
  "Interferon module" = c("Ifit1", "Ifit3", "Ifit1bl1", "Ifit3b", "Ifi47", "Rsad2"),
  "Stress response / proliferation" = c("Hist1h1c", "Hist1h2ap", "Hspa1b", "Hspa1a")
)

############################################################
# Select samples to include in the heatmap
############################################################

# Define sample columns from count matrix
expr_cols <- c("effector1", "effector2", "ctrl1", "ctrl2")

############################################################
# Subset count matrix to significant genes
############################################################

counts_sub <- counts %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::filter(gene_name %in% sig_genes) %>%
  dplyr::select(gene_name, all_of(expr_cols))

# Stop if none of the significant genes are present in the count matrix
if (nrow(counts_sub) == 0) {
  stop("No genes from effector_ctrl_IHWsigGenes.tsv were found in counts$gene_name.")
}

# Remove duplicated gene names, keeping the first occurrence
counts_sub <- counts_sub %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

############################################################
# Build expression matrix
############################################################

# Move gene names to row names
mat_plot <- counts_sub %>%
  tibble::column_to_rownames("gene_name")

# Convert to numeric matrix
mat_plot <- as.matrix(mat_plot)
mode(mat_plot) <- "numeric"

# Keep genes in the same order as the significant gene file
sig_genes_present <- sig_genes[sig_genes %in% rownames(mat_plot)]
mat_plot <- mat_plot[sig_genes_present, , drop = FALSE]

# Log2-transform counts
mat_plot <- log2(mat_plot + 1)

# Remove genes with no variation across samples
keep <- apply(mat_plot, 1, function(x) sd(x, na.rm = TRUE) > 0)
mat_plot <- mat_plot[keep, , drop = FALSE]

# Stop if all genes have zero variance
if (nrow(mat_plot) == 0) {
  stop("All selected genes have zero variance across samples.")
}

# Scale expression per gene using row z-score
mat_scaled <- t(scale(t(mat_plot)))

############################################################
# Build row annotation dataframe
############################################################

# Store all genes present in scaled matrix
all_genes <- rownames(mat_scaled)

# Create dataframe assigning selected genes to predefined modules
row_annot_df <- bind_rows(lapply(names(module_list), function(mod) {
  data.frame(
    gene = module_list[[mod]],
    module = mod,
    stringsAsFactors = FALSE
  )
})) %>%
  dplyr::filter(gene %in% all_genes)

# Identify significant genes not assigned to predefined modules
other_genes <- base::setdiff(all_genes, unlist(module_list))

# Assign remaining genes to "Other"
other_df <- data.frame(
  gene = other_genes,
  module = "Other",
  stringsAsFactors = FALSE
)

# Combine module genes and other genes
row_annot_df <- bind_rows(row_annot_df, other_df)

############################################################
# Define gene order in the heatmap
############################################################

# Put module genes first, then other significant genes
module_gene_order <- unlist(module_list)
other_gene_order <- sig_genes_present[sig_genes_present %in% other_genes]

desired_gene_order <- c(
  module_gene_order[module_gene_order %in% row_annot_df$gene],
  other_gene_order
)

desired_gene_order <- desired_gene_order[desired_gene_order %in% row_annot_df$gene]

# Apply gene order to annotation dataframe
row_annot_df$gene <- factor(row_annot_df$gene, levels = desired_gene_order)
row_annot_df <- row_annot_df %>%
  dplyr::arrange(gene)

# Define order of module categories
row_annot_df$module <- factor(
  row_annot_df$module,
  levels = c(
    "TCR repertoire / clonal signature",
    "Effector CD8 T cell module",
    "Interferon module",
    "Stress response / proliferation",
    "Other"
  )
)

# Reorder matrix rows to match row annotation
mat_scaled <- mat_scaled[as.character(row_annot_df$gene), , drop = FALSE]

# Reverse sample column order
mat_scaled <- mat_scaled[, rev(colnames(mat_scaled))]

############################################################
# Create column annotation
############################################################

# Define sample condition annotation
annot_col <- data.frame(
  condition = c( "ctrl_ici", "ctrl_ici","effector_ici", "effector_ici"),
  row.names = c("ctrl1", "ctrl2","effector1", "effector2")
)

# Reorder annotation to match matrix columns
annot_col <- annot_col[colnames(mat_scaled), , drop = FALSE]

############################################################
# Define heatmap colors
############################################################

# Colors for gene modules
module_colors <- c(
  "TCR repertoire / clonal signature" = "#4DAF4A",
  "Effector CD8 T cell module" = "#377EB8",
  "Interferon module" = "#E41A1C",
  "Stress response / proliferation" = "#984EA3",
  "Other" = "gray70"
)

# Colors for sample conditions
condition_colors <- c(
  "ctrl_ici" = "#00D9DE",
  "effector_ici" = "#FF9289"
)

# Color scale for z-score expression
col_fun <- colorRamp2(c(-2, 0, 2), c("#6A99C5", "#FAFDC7", "#DC3F2E"))

############################################################
# Create heatmap annotations
############################################################

# Left annotation shows gene module per row
left_annot <- rowAnnotation(
  Module = row_annot_df$module,
  col = list(Module = module_colors),
  show_annotation_name = FALSE,
  width = unit(6, "mm")
)

# Recreate column annotation
annot_col <- data.frame(
  condition = c( "ctrl_ici", "ctrl_ici","effector_ici", "effector_ici"),
  row.names = c("ctrl1", "ctrl2","effector1", "effector2")
)

# Top annotation shows sample condition
top_annot <- HeatmapAnnotation(
  condition = annot_col$condition,
  col = list(condition = condition_colors),
  show_annotation_name = TRUE
)

############################################################
# Build and draw heatmap
############################################################

ht <- Heatmap(
  mat_scaled,
  name = "z-score",
  col = col_fun,
  left_annotation = left_annot,
  top_annotation = top_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = row_annot_df$module,
  cluster_row_slices = FALSE,
  row_title_rot = 90,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 14, fontface = "italic"),
  column_names_gp = gpar(fontsize = 14),
  border = TRUE,
  column_title = "effector vs ctrl"
)

# Draw heatmap with legends on the right
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

############################################################
# Optional: save heatmap as SVG
############################################################

library(svglite)

#svglite(
#  "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/heatmap_dge_vertical.svg",
#  width = 8,
#  height = 14,
#  bg = "white"
#)

# Draw heatmap again, intended for SVG export
draw(ht)

# Close graphics device
dev.off()