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


resIHW <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWallGenes.tsv")
counts <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/counts_matrix.tsv")
sample_sheet <- read_csv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/samplesheet.csv")
resIHWsig <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWsigGenes.tsv") %>%
  dplyr::filter(gene_type != "antisense")
resIHWsig_fc <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_vs_ctrl_v0/effector_ctrl_IHWsigFCgenes_2_fold.tsv") 
results_dir <- "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway"
prefix = "example1"

save_plot <- function(filename, p, width = NULL, height = NULL) {
  if (!is.null(width) && !is.null(height)) {
    ggsave(file.path(paste0(filename, ".png")), plot = p, width = width, height = height)
    ggsave(file.path(paste0(filename, ".svg")), plot = p, width = width, height = height)
  } else {
    ggsave(file.path(paste0(filename, ".png")), plot = p)
    ggsave(file.path(paste0(filename, ".svg")), plot = p)
  }
}

#### Get parameters from docopt

# Input and output
#sampleAnnotationCSV <- arguments$sample_sheet
#readCountFile <- arguments$count_table
#results_dir <- arguments$result_dir
#dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
#paired_grp <- arguments$paired_grp


###### Run TOPGO analysis
de_symbols <- resIHWsig_fc$gene_id
#bg_symbols <- rownames(dds)[rowSums(counts(dds)) > 0]

#lapply(c("BP", "MF", "CC"), function(ontology) {
#  topgoDE <- topGOtable(de_symbols, bg_symbols,
#    ontology = ontology,
#    mapping = anno_db,
#    geneID = gene_id_type
#  )
#  write_tsv(topgoDE, file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".tsv")))
#  write_xlsx(topgoDE %>% select(-genes), file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".xlsx")))
#})


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

## ORA
ora_tests <- list(
  "KEGG" = function(genes, universe) {
    enrichKEGG(
      gene         = genes,
      universe     = universe,
      organism     = org_kegg,
      pvalueCutoff = 0.05
    )
  },
  "Reactome" = function(genes, universe) {
    enrichPathway(
      gene = genes,
      organism = org_reactome,
      universe = universe,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  },
  "WikiPathway" = function(genes, universe) {
    enrichWP(
      gene = genes,
      universe = universe,
      organism = org_wp,
      pvalueCutoff = 0.05
    )
  },
  "GO_BP" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  },
  "GO_MF" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "MF",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  },
  "GO_CC" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "CC",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  }
)

# Warmup GO database - work around https://github.com/YuLab-SMU/clusterProfiler/issues/207
._ <- enrichGO(universe[1], OrgDb = get(anno_db), keyType = "ENTREZID", ont = "BP", universe = universe)

get_heatplot_dims <- function(p) {
  nr_gene <- length(unique(p$data$Gene))
  nr_cat <- length(unique(p$data$categoryID))

  hp_width <- min(nr_gene * 0.25, 40)
  hp_height <- min(nr_cat * 0.25, 40)

  return(c(hp_width, hp_height))
}

bplapply(names(ora_tests), function(ora_name) {
  message(paste0("Performing ", ora_name, " ORA-test..."))

  test_fun <- ora_tests[[ora_name]]
  ora_res <- test_fun(resIHWsig_fc_entrez$ENTREZID, universe)

  if (!is.null(ora_res)) {
    ora_res <- setReadable(ora_res, OrgDb = anno_db, keyType = "ENTREZID")
    res_tab <- as_tibble(ora_res@result)
    write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_ORA_", ora_name, ".tsv")))

    if (min(res_tab$p.adjust) < 0.05 & length(unique(res_tab$geneID)) > 1) {
      # Get the y-axis labels from the ora_res object
      y_labels <- ora_res$Description # Replace with the actual way to extract your labels

      # Find the maximum length of the labels
      max_length <- max(nchar(y_labels))

      # Set a wrapping width based on the maximum label length
      wrap_width <- ifelse(max_length > 50, 50, max_length) # Adjust as needed

      # Create the dotplot with automatic label adjustment
      p <- dotplot(ora_res, showCategory = 30) +
        theme(
          plot.margin = unit(c(1, 1, 1, 1), "cm"), # Maintain the left margin
          axis.text.y = element_text(size = 13, hjust = 1)
        ) + # Keep the text size
        scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) # Dynamic wrapping

      save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_dotplot")), p, width = 15, height = 10)

      p <- cnetplot(ora_res,
        #categorySize = "pvalue",
        showCategory = min(5, length(ora_res@result$ID)),
        foldChange = de_foldchanges,
        #vertex.label.font = 6
      ) + scale_color_gradient2(name = "log2FoldChange", low = "blue3", high = "firebrick")

      save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_cnetplot")), p, width = 15, height = 12)

      if (!skip_heatplots) {
        p <- heatplot(ora_res, foldChange = de_foldchanges, showCategory = 40) +
          scale_fill_gradient2(midpoint = 0, low = "blue4", mid = "white", high = "red4")
        hp_dims <- get_heatplot_dims(p)

        save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_heatplot")), p, width = hp_dims[1], height = hp_dims[2])
      }
    } else {
      message(paste0("Warning: No significant enrichment in ", ora_name, " ORA analysis. "))
    }
  } else {
    message(paste0("Warning: No gene can be mapped in ", ora_name, " ORA analysis. "))
  }
})


## GSEA
skip_gsea = FALSE
if (!skip_gsea) {
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

  bplapply(names(gsea_tests), function(gsea_name) {
    message(paste0("Performing ", gsea_name, " GSEA-test..."))

    test_fun <- gsea_tests[[gsea_name]]
    gsea_res <- test_fun(ranked_gene_list)

    if (!is.null(gsea_res)) {
      gsea_res <- setReadable(gsea_res, OrgDb = get(anno_db), keyType = "ENTREZID")
      res_tab <- gsea_res@result %>% as_tibble()

      write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, ".tsv")))
      if (min(res_tab$p.adjust) < 0.05) {
        # Get the y-axis labels from the gsea_res object
        y_labels <- gsea_res$Description # Replace with the actual way to extract your labels

        # Find the maximum length of the labels
        max_length <- max(nchar(y_labels))

        # Set a wrapping width based on the maximum label length
        wrap_width <- ifelse(max_length > 50, 50, max_length) # Adjust as needed

        # Filter significant results (p.adjust < 0.05)
        gsea_res_sig <- gsea_res
        gsea_res_sig@result <- gsea_res@result[gsea_res@result$p.adjust < 0.05, ]

        # Check if any significant results remain
        if (nrow(gsea_res_sig@result) == 0) {
          print("No significantly enriched gene sets found with p.adjust < 0.1")
        } else {
          # Create the dotplot with automatic label adjustment
          p <- dotplot(gsea_res_sig, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) +
            theme(
              plot.margin = unit(c(1, 1, 1, 1), "cm"), # Maintain the left margin
              strip.text.x = element_text(size = 16, face = "bold"),
              axis.text.y = element_text(size = 13, hjust = 1)
            ) + # Keep the text size
            scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) # Dynamic wrapping

          save_plot(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_dotplot")), p, width = 20, height = 15)
        }

        p <- cnetplot(gsea_res,
          #categorySize = "pvalue",
          showCategory = 5,
          foldChange = de_foldchanges,
          #vertex.label.font = 6
        ) + scale_color_gradient2(name = "log2FoldChange", low = "blue3", high = "firebrick")

        save_plot(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_cnetplot")), p, width = 15, height = 12)

        # GSEA generates to long gene lists so that the heatplot gets to overloaded
        # p <- heatplot(gsea_res, foldChange=de_foldchanges, showCategory=40) +
        #   scale_fill_gradient2(midpoint=0, low="blue4", mid="white", high="red4" )
        #
        # hp_dims <- get_heatplot_dims(p)
        #
        # ggsave(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_heatplot.png")), plot = p, width = hp_dims[1], height = hp_dims[2])
      } else {
        message(paste0("Warning: No significant enrichment in ", gsea_name, " GSEA analysis. "))
      }
    } else {
      message(paste0("Warning: No gene can be mapped in ", gsea_name, " GSEA analysis. "))
    }
  })
}


############################################ REACTOME ACTIVATED
names(gsea_tests)
gsea_name="Reactome"

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
  ggtitle("Reactome - Activated") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text.y = element_text(size = 13, hjust = 1),
    axis.title.y = element_blank()
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width))

p

#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/Reactome_Activated_effector_vs_ctrl.png",
#  plot = p,
#  width = 8,
#  height = 6,
#  dpi = 300
#)
#
## SVG (vector format, dpi not needed but can be set)
#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/Reactome_Activated_effector_vs_ctrl.svg",
#  plot = p,
#  width = 8,
#  height = 6
#)

############################################ WIKIPATHWAY ACTIVATED

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

#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/Wikipathway_Activated_effector_vs_ctrl.png",
#  plot = p,
#  width = 8,
#  height = 6,
#  dpi = 300
#)
#
## SVG (vector format, dpi not needed but can be set)
#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/Wikipathway_Activated_effector_vs_ctrl.svg",
#  plot = p,
#  width = 8,
#  height = 6
#)



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



#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/GO_BP_Activated_effector_vs_ctrl.png",
#  plot = p,
#  width = 8,
#  height = 6,
#  dpi = 300
#)
#
## SVG (vector format, dpi not needed but can be set)
#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/GO_BP_Activated_effector_vs_ctrl.svg",
#  plot = p,
#  width = 8,
#  height = 6
#)
#


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


#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/GO_MF_effector_vs_ctrl.png",
#  plot = p,
#  width = 10,
#  height = 6,
#  dpi = 300
#)
#
## SVG (vector format, dpi not needed but can be set)
#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/GO_MF_effector_vs_ctrl.svg",
#  plot = p,
#  width = 10,
#  height = 6
#)


##########################


library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(stringr)

gsea_name = "Reactome"

test_fun <- gsea_tests[[gsea_name]]
gsea_res <- test_fun(ranked_gene_list)

# Filter significant results
gsea_res_sig <- gsea_res
gsea_res_sig@result <- gsea_res@result %>%
  filter(p.adjust < 0.05)

# Keep activated pathways
gsea_res_activated <- gsea_res_sig
gsea_res_activated@result <- gsea_res_sig@result %>%
  filter(NES > 0)

# Wrap labels
y_labels <- gsea_res_activated@result$Description
max_length <- max(nchar(y_labels))
wrap_width <- ifelse(max_length > 50, 50, max_length)

# Plot
p <- dotplot(
  gsea_res_activated,
  showCategory = 30,
  x = "p.adjust",
  color = "NES"
) +
  ggtitle("Reactome") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text.y = element_text(size = 13, hjust = 1),
    axis.title.y = element_blank()
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Enrichment Score"
  )

p

ggsave(
  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/pathway_NES/go_bp_effector_vs_ctrl.png",
  plot = p,
  width = 8.5,
  height = 6,
  dpi = 300
)

# SVG (vector format, dpi not needed but can be set)
ggsave(
  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/pathway_NES/go_bp_effector_vs_ctrl.svg",
  plot = p,
  width = 8.5,
  height = 6
)

library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)

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

#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/pathway_NES/wikipathway_effector_vs_ctrl.png",
#  plot = p,
#  width = 9,
#  height = 6,
#  dpi = 300
#)

#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/pathway_NES/wikipathway_effector_vs_ctrl.svg",
#  plot = p,
#  width = 9,
#  height = 6,
#  dpi = 300
#)


###

gomf <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_ICI_vs_ctrl_ICI/effector_ICI_ctrl_ICI_GSEA_GO_MF.tsv")
gobp <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_ICI_vs_ctrl_ICI/effector_ICI_ctrl_ICI_GSEA_GO_BP.tsv")
reactome <- read_tsv("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/effector_ICI_vs_ctrl_ICI/effector_ICI_ctrl_ICI_GSEA_Reactome.tsv")
top_terms <- gsea_res@result %>%
  filter(p.adjust < 0.05) %>%
  slice_min(order_by = p.adjust, n = 10)


subset <- reactome %>%
  filter(Description %in% top_terms$Description)

gene_freq <- subset %>%
  select(Description, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/") %>%
  count(core_enrichment, sort = TRUE)

head(gene_freq, 30)

library(dplyr)

keep_terms <- c(
  "Cross-presentation of soluble exogenous antigens (endosomes)",
  "Antigen processing-Cross presentation",
  "Interleukin-1 family signaling"
)


gsea_subset@result <- gsea_res_sig@result %>%
  filter(Description %in% keep_terms)

gsea_subset_readable <- setReadable(
  gsea_subset,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID"
)

cnetplot(
  gsea_subset_readable,
  showCategory = 10,
  foldChange = ranked_gene_list
)


###############

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(pheatmap)

# -----------------------------
# 1. Extract unique core enrichment genes
# -----------------------------
core_genes <- gsea_subset_readable@result %>%
  dplyr::select(core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::pull(core_enrichment) %>%
  unique()

core_genes <- core_genes[core_genes %in% resIHWsig$gene_name]

length(core_genes)
head(core_genes)

# Use your real condition/sample columns here
expr_cols <- c(
  "GF_noICI1", "GF_noICI2",
  "GF1", "GF2",
  "ctrl1", "ctrl2",
  "effector1", "effector2"
)

counts_sub <- counts %>%
  filter(!is.na(gene_name)) %>%
  filter(gene_name %in% core_genes) %>%
  select(gene_name, all_of(expr_cols)) %>%
  distinct(gene_name, .keep_all = TRUE)

nrow(counts_sub)

mat_plot <- counts_sub %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

mat_plot_log <- log2(mat_plot + 1)

mat_plot_scaled <- t(scale(t(mat_plot_log)))

mat_plot_scaled[is.na(mat_plot_scaled)] <- 0

mode(mat_plot) <- "numeric"
annotation_col <- data.frame(
  Condition = c(
    "GF_noici", "GF_noici",
    "GF", "GF",
    "ctrl", "ctrl",
    "effector", "effector"
  )
)

rownames(annotation_col) <- colnames(mat_plot_scaled)
library(pheatmap)

p <- pheatmap(
  mat_plot_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  main = "Core enrichment genes from top GSEA pathways"
)


#ggsave(
#  filename = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/results/analysis/pathway/pathway_NES/heatmap_core_gomf.svg",
#  plot = p,
#  width = 5,
#  height = 6,
#  dpi = 300
#)


##########

library(dplyr)
library(tibble)
library(pheatmap)

genes_to_plot <- c("Hspa1a", "Hist1h2ap", "Hspa1b", "Hist1h1c")

counts_sub_genes <- counts %>%
  filter(gene_name %in% genes_to_plot)

mat_plot <- counts_sub_genes %>%
  select(-gene_id) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

mode(mat_plot) <- "numeric"

# enforce sample order
sample_order <- c(
  "GF_noICI1", "GF_noICI2",
  "GF1", "GF2",
  "ctrl1", "ctrl2",
  "effector1", "effector2"
)

mat_plot <- mat_plot[, sample_order]

# log transform
mat_plot_log <- log2(mat_plot + 1)

# row scaling
mat_plot_scaled <- t(scale(t(mat_plot_log)))
mat_plot_scaled[is.na(mat_plot_scaled)] <- 0

# annotation
annotation_col <- data.frame(
  Condition = c(
    "GF_noICI", "GF_noICI",
    "GF", "GF",
    "ctrl", "ctrl",
    "effector", "effector"
  )
)

rownames(annotation_col) <- sample_order

# heatmap
pheatmap(
  mat_plot_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  main = "Heat Shock related genes p.adj < 0.1"
)
