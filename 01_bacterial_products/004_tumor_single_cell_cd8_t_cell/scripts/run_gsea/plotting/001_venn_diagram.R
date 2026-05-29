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

# Prefer dplyr versions of commonly conflicting functions
conflict_prefer("lag", "dplyr")
conflict_prefer("union", "base")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("first", "dplyr")

library("EnhancedVolcano")
library("ggpubr")
library("tibble")
library("stringr")
library("ggrepel")
library("grateful")

# Additional conflict preferences
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
conflict_prefer("as.factor", "base")

# Helper function to remove Ensembl version numbers
remove_ensg_version <- function(x) gsub("\\.[0-9]*$", "", x)

library(conflicted)

# Explicitly prefer base set operations
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::union)

############################################################
# Load differential expression result tables
############################################################

# Effector ICI vs Ctrl ICI
resIHW_e_vs_c <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_ICI_vs_ctrl_ICI/effector_ICI_ctrl_ICI_IHWsigGenes.tsv")

# Effector ICI vs GF ICI
resIHW_e_vs_gf <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_ICI_vs_GF_ICI/effector_ICI_GF_ICI_IHWsigGenes.tsv")

# GF ICI vs GF no ICI
resIHW_gf_vs_gfnoici <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/GF_ICI_vs_GF_noICI/GF_ICI_GF_noICI_IHWsigGenes.tsv")

# Ctrl ICI vs GF ICI
resIHW_c_vs_gf <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/ctrl_ICI_vs_GF_ICI/ctrl_ICI_GF_ICI_IHWsigGenes.tsv")

############################################################
# Remove unwanted transcript biotypes
############################################################

# Remove processed transcripts and antisense genes
resIHW_e_vs_c <- resIHW_e_vs_c %>%
  filter(!gene_type %in% c("processed_transcript", "antisense"))

resIHW_e_vs_gf <- resIHW_e_vs_gf %>%
  filter(!gene_type %in% c("processed_transcript", "antisense"))

resIHW_gf_vs_gfnoici <- resIHW_gf_vs_gfnoici %>%
  filter(!gene_type %in% c("processed_transcript", "antisense"))

resIHW_c_vs_gf <- resIHW_c_vs_gf %>%
  filter(!gene_type %in% c("processed_transcript", "antisense"))

############################################################
# Load count matrix and sample annotation
############################################################

counts <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/run_gsea/tables/cd8_all_effector_ctrl/counts_matrix.tsv")

sample_anno <- read_csv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/run_gsea/tables/cd8_all_effector_ctrl/samplesheet.csv")

############################################################
# Create gene labels indicating direction of regulation
############################################################

# Function to append "_up" or "_down" to gene names
add_gene_label <- function(df) {
  df %>%
    mutate(
      gene_name_label = case_when(
        log2FoldChange > 0 ~ paste0(gene_name, "_up"),
        log2FoldChange < 0 ~ paste0(gene_name, "_down"),
        TRUE ~ paste0(gene_name, "_zero")
      )
    )
}

# Apply labeling function to all DEG tables
resIHW_e_vs_c <- add_gene_label(resIHW_e_vs_c)
resIHW_e_vs_gf <- add_gene_label(resIHW_e_vs_gf)
resIHW_gf_vs_gfnoici <- add_gene_label(resIHW_gf_vs_gfnoici)
resIHW_c_vs_gf <- add_gene_label(resIHW_c_vs_gf)

############################################################
# Create gene sets for Venn diagram analysis
############################################################

venn_list <- list(
  effector_vs_ctrl = unique(resIHW_e_vs_c$gene_name_label),
  effector_vs_GF = unique(resIHW_e_vs_gf$gene_name_label),
  GF_ICI_vs_GF_noICI = unique(resIHW_gf_vs_gfnoici$gene_name_label),
  ctrl_vs_GF = unique(resIHW_c_vs_gf$gene_name_label)
)

############################################################
# Plot Venn diagram
############################################################

library(ggVennDiagram)

ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue")

############################################################
# Function to extract specific Venn regions
############################################################

library(purrr)

get_venn_region <- function(venn_list, include, exclude = NULL) {
  
  # Find intersection of included sets
  res <- Reduce(base::intersect, venn_list[include])
  
  # Remove genes present in excluded sets
  if (!is.null(exclude)) {
    res <- base::setdiff(res, Reduce(base::union, venn_list[exclude]))
  }
  
  res
}

############################################################
# Extract specific overlaps from Venn diagram
############################################################

# Genes unique to effector_vs_ctrl
only_effector_vs_ctrl <- get_venn_region(
  venn_list,
  include = "effector_vs_ctrl",
  exclude = c("effector_vs_GF", "ctrl_vs_GF")
)

# Shared between effector_vs_ctrl and effector_vs_GF
# but absent from ctrl_vs_GF
effector_ctrl__and__effector_gf <- get_venn_region(
  venn_list,
  include = c("effector_vs_ctrl", "effector_vs_GF"),
  exclude = "ctrl_vs_GF"
)

# Shared across all three comparisons
all_three <- get_venn_region(
  venn_list,
  include = c("effector_vs_ctrl", "effector_vs_GF", "ctrl_vs_GF")
)

############################################################
# Store all Venn regions in a list
############################################################

venn_genes <- list(
  
  only_effector_vs_ctrl = get_venn_region(
    venn_list,
    include = "effector_vs_ctrl",
    exclude = c("effector_vs_GF", "ctrl_vs_GF")
  ),
  
  only_effector_vs_GF = get_venn_region(
    venn_list,
    include = "effector_vs_GF",
    exclude = c("effector_vs_ctrl", "ctrl_vs_GF")
  ),
  
  only_ctrl_vs_GF = get_venn_region(
    venn_list,
    include = "ctrl_vs_GF",
    exclude = c("effector_vs_ctrl", "effector_vs_GF")
  ),
  
  effector_ctrl_and_effector_gf = get_venn_region(
    venn_list,
    include = c("effector_vs_ctrl", "effector_vs_GF"),
    exclude = "ctrl_vs_GF"
  ),
  
  effector_ctrl_and_ctrl_gf = get_venn_region(
    venn_list,
    include = c("effector_vs_ctrl", "ctrl_vs_GF"),
    exclude = "effector_vs_GF"
  ),
  
  effector_gf_and_ctrl_gf = get_venn_region(
    venn_list,
    include = c("effector_vs_GF", "ctrl_vs_GF"),
    exclude = "effector_vs_ctrl"
  ),
  
  gf_ici_vs_gf_noici = get_venn_region(
    venn_list,
    include = c("GF_ICI_vs_GF_noICI"),
    exclude = 
  ),
  
  all_three = get_venn_region(
    venn_list,
    include = c("effector_vs_ctrl", "effector_vs_GF", "ctrl_vs_GF")
  )
)

############################################################
# Inspect selected gene sets
############################################################

venn_genes$only_effector_vs_ctrl
venn_genes$only_effector_vs_GF
venn_genes$only_ctrl_vs_GF
venn_genes$gf_ici_vs_gf_noici
