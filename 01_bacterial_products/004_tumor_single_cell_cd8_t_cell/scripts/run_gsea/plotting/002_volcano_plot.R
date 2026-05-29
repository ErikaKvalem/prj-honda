############################################################
# Load required libraries and resolve function conflicts
############################################################

library("conflicted")
library("BiocParallel")
library("DESeq2")
library("IHW")
library("ggplot2")
library("pcaExplorer")
library("topGO")
library("clusterProfiler")
library("ReactomePA")
library("writexl")
library("readr")
library("dplyr")

# Prefer specific functions from dplyr when conflicts exist
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("first", "dplyr")

library("EnhancedVolcano")
library("ggpubr")
library("tibble")
library("stringr")
library("ggrepel")
library("grateful")

# Prefer base/dplyr versions of additional conflicting functions
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
conflict_prefer("as.factor", "base")

# Function to remove version suffixes from Ensembl gene IDs
remove_ensg_version <- function(x) gsub("\\.[0-9]*$", "", x)

############################################################
# Load differential expression results
############################################################

# Load all genes from effector vs ctrl DESeq2/IHW result table
resIHW <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWallGenes.tsv")

# Remove antisense genes
resIHW <- resIHW[resIHW$gene_type != "antisense", ]

############################################################
# Filter significant DE genes
############################################################

# Keep genes with abs(log2FC) > 0.58 and adjusted p-value < 0.1
resIHW_filtered <- resIHW[
  abs(resIHW$log2FoldChange) > 0.58 &
    resIHW$padj < 0.1,
]

############################################################
# Load count matrix and significant genes table
############################################################

# Load count matrix
counts <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/counts_matrix.tsv")

# Load significant genes table and remove antisense genes
sig_genes_df <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWsigGenes.tsv") %>%
  dplyr::filter(gene_type != "antisense")

############################################################
# Define volcano plot cutoffs and genes of interest
############################################################

# Define log2 fold-change and FDR cutoffs
fc_cutoff = 0.58
fdr_cutoff= 0.1

############################################################
# Volcano plot: second version labeling significant genes
############################################################

# Extract genes passing fold-change and adjusted p-value thresholds
genes_to_label <- resIHW_filtered$gene_name[
  abs(resIHW$log2FoldChange) > 0.58 &
    resIHW$padj < 0.1
]



# Create volcano plot labeling significant genes
p <- EnhancedVolcano(
  resIHW,
  lab = resIHW$gene_name,
  #  selectLab = goi$gene_name,
  selectLab =  genes_to_label,
  labSize = 8,
  drawConnectors = TRUE,
  colConnectors = "black",
  x = "log2FoldChange",
  y = "padj",
  pCutoff = fdr_cutoff,
  FCcutoff = fc_cutoff,
  subtitle = "",
  legendPosition = "right",
  caption = paste0(
    "fold change cutoff: ",
    round(2**fc_cutoff, 1),
    ", adj.p-value cutoff: ",
    fdr_cutoff
  ),
  maxoverlapsConnectors = Inf,
  title = "effector vs ctrl",
  labFace = "italic",
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
)

# Limit volcano plot y-axis
p <- p + coord_cartesian(ylim = c(0, 10))

# Display final volcano plot
p

############################################################
# Optional: save volcano plot as SVG
############################################################

#ggsave(
#  "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_ICI_vs_ctrl_ICI/volcano_log2foldchange058.svg",
#  plot = p,
#  width = 12,
#  height = 12,
#  dpi = 300,
#  bg = "white"
#)

