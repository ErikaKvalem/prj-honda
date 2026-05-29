#!/usr/bin/env Rscript


library("conflicted")
library("docopt")
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
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
conflict_prefer("as.factor", "base")
remove_ensg_version <- function(x) gsub("\\.[0-9]*$", "", x)
anno_db <- "org.Mm.eg.db"
library(anno_db, character.only = TRUE)


org_kegg <- "mmu"
org_reactome <- "mouse"
org_wp <- "Mus musculus"
fc_cutoff = 1
fdr_cutoff= 0.1
prefix = "effector vs ctrl"


resIHW <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWallGenes.tsv")
counts <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/counts_matrix.tsv")
sample_sheet <- read_csv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/samplesheet.csv")
resIHWsig <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWsigGenes.tsv") %>%
  dplyr::filter(gene_type != "antisense")
resIHWsig_fc <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWsigFCgenes_2_fold.tsv") 
results_dir <- "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway"


###### Run TOPGO analysis
de_symbols <- resIHWsig_fc$gene_id


##### Pathway enrichment analysis
hgnc_to_entrez <- AnnotationDbi::select(get(anno_db), resIHW %>% pull("gene_name") %>% unique(), keytype = "SYMBOL", columns = c("ENTREZID"))

# full list with ENTREZIDs added
resIHW_entrez <- resIHW %>% inner_join(hgnc_to_entrez, by = c("gene_name" = "SYMBOL"))
universe <- resIHW_entrez %>%
  pull("ENTREZID") %>%
  unique()

# list of significant genes with ENTREZIDs added
resIHWsig_fc_entrez <- resIHWsig_fc %>% inner_join(hgnc_to_entrez, by = c("gene_name" = "SYMBOL"))
de_foldchanges <- resIHWsig_fc_entrez$log2FoldChange
names(de_foldchanges) <- resIHWsig_fc_entrez$ENTREZID


## GSEA
# for GSEA use genes ranked by test statistic
res_ihw_ranked <- resIHW_entrez %>%
  arrange(-stat) %>%
  select(ENTREZID, stat) %>%
  na.omit() %>%
  distinct(ENTREZID, .keep_all = TRUE)

ranked_gene_list <- res_ihw_ranked$stat

names(ranked_gene_list) <- res_ihw_ranked$ENTREZID

gsea_tests <- list(
  "KEGG" = function(ranked_gene_list) {
    gseKEGG(geneList = ranked_gene_list, organism = org_kegg, pvalueCutoff = 1)
  },
  "Reactome" = function(ranked_gene_list) {
    gsePathway(geneList = ranked_gene_list, organism = org_reactome, pvalueCutoff = 1)
  },
  "WikiPathway" = function(ranked_gene_list) {
    gseWP(geneList = ranked_gene_list, organism = org_wp, pvalueCutoff = 1)
  },
  "GO_BP" = function(ranked_gene_list) {
    gseGO(
      geneList = ranked_gene_list,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      minGSSize = 10
    )
  },
  "GO_MF" = function(ranked_gene_list) {
    gseGO(
      geneList = ranked_gene_list,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      minGSSize = 10
    )
  },
  "GO_CC" = function(ranked_gene_list) {
    gseGO(
      geneList = ranked_gene_list,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      minGSSize = 10
    )
  }
)

names(gsea_tests)

################################################################### Reactome
gsea_name <- "Reactome"

test_fun <- gsea_tests[[gsea_name]]
gsea_res <- test_fun(ranked_gene_list)

plot_df <- gsea_res@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    Direction = ifelse(NES > 0, "Activated", "Suppressed")
  ) %>%
  group_by(Direction) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Description_wrapped = str_wrap(Description, width = 50),
    Description_wrapped = forcats::fct_reorder(Description_wrapped, NES)
  )

nrow(plot_df)

p <- ggplot(
  plot_df,
  aes(
    x = p.adjust,
    y = Description_wrapped,
    size = setSize,
    color = NES
  )
) +
  geom_point(alpha = 0.9) +
  facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Enrichment Score"
  ) +
  scale_x_reverse() +
  labs(
    title = "Reactome",
    x = "Adjusted p-value",
    y = NULL,
    size = "Set size"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

p



################################################################### Wikipathway 


names(gsea_tests)
gsea_name="WikiPathway"

test_fun <- gsea_tests[[gsea_name]]
gsea_res <- test_fun(ranked_gene_list)

y_labels <- gsea_res$Description # Replace with the actual way to extract your labels

# Find the maximum length of the labels
max_length <- max(nchar(y_labels))

# Set a wrapping width based on the maximum label length
wrap_width <- ifelse(max_length > 50, 50, max_length) # Adjust as needed

# Filter significant results (p.adjust < 0.05)
gsea_res_sig <- gsea_res
gsea_res_sig@result <- gsea_res@result[gsea_res@result$p.adjust < 0.05, ]


gsea_res_activated <- gsea_res_sig
gsea_res_activated@result <- gsea_res_sig@result %>%
  dplyr::filter(NES > 0)


p <- dotplot(gsea_res_activated, showCategory = 15) +
  ggtitle("WikiPathway - Activated") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text.y = element_text(size = 13, hjust = 1),
    axis.title.y = element_blank()
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width))

############################################ GO_BP ACTIVATED


names(gsea_tests)
gsea_name="GO_BP"

test_fun <- gsea_tests[[gsea_name]]
gsea_res <- test_fun(ranked_gene_list)

y_labels <- gsea_res$Description # Replace with the actual way to extract your labels

# Find the maximum length of the labels
max_length <- max(nchar(y_labels))

# Set a wrapping width based on the maximum label length
wrap_width <- ifelse(max_length > 50, 50, max_length) # Adjust as needed

# Filter significant results (p.adjust < 0.05)
gsea_res_sig <- gsea_res
gsea_res_sig@result <- gsea_res@result[gsea_res@result$p.adjust < 0.05, ]


gsea_res_activated <- gsea_res_sig
gsea_res_activated@result <- gsea_res_sig@result %>%
  dplyr::filter(NES > 0)


p <- dotplot(gsea_res_activated, showCategory = 15) +
  ggtitle("GO BP - Activated") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text.y = element_text(size = 13, hjust = 1),
    axis.title.y = element_blank()
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width))

############################################ GO_MF ACTIVATED


names(gsea_tests)
gsea_name="GO_MF"

test_fun <- gsea_tests[[gsea_name]]
gsea_res <- test_fun(ranked_gene_list)

y_labels <- gsea_res$Description # Replace with the actual way to extract your labels

# Find the maximum length of the labels
max_length <- max(nchar(y_labels))

# Set a wrapping width based on the maximum label length
wrap_width <- ifelse(max_length > 50, 50, max_length) # Adjust as needed

# Filter significant results (p.adjust < 0.05)
gsea_res_sig <- gsea_res
gsea_res_sig@result <- gsea_res@result[gsea_res@result$p.adjust < 0.05, ]


gsea_res_activated <- gsea_res_sig
gsea_res_activated@result <- gsea_res_sig@result %>%
  dplyr::filter(NES > 0)


p <- dotplot(gsea_res_sig, showCategory = 8, split = ".sign") + facet_grid(. ~ .sign) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Maintain the left margin
    strip.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 13, hjust = 1)
  ) + # Keep the text size
  scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) # Dynamic wrapping

p


