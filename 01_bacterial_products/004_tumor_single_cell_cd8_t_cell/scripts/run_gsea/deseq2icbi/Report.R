#!/usr/bin/env Rscript

' Report.R

# Generate a LaTeX report of DESeq2 analysis results.

Usage:
  Report.R <deseq2_res_dir> <comparison> <bib_file> <template_file> [options]

Arguments:
  <deseq2_res_dir>  Path to the directory containing DESeq2 results.
  <comparison>      Name of the comparison, e.g., "groupA_vs_groupB".
  <bib_file>        Path to the BibTeX file for references.
  <template_file>   Path to the LaTeX template file.

Options:
  --paired_grp=<paired_grp>     Column name in the sample sheet indicating paired samples (optional).
  --remove_batch_effect         Whether batch effect correction was applied (TRUE/FALSE).
  --batch_col=<batch_col>       Optional: column in sample annotation that contains the batch
  --covariate=<covariate>       Covariate formula used in the DESeq2 analysis (optional).
' -> doc

library("conflicted") # Manage function name conflicts
library("docopt") # Command-line argument parser
library("dplyr") # Data manipulation
library("tools") # File path manipulation
library("whisker") # Templating engine

# Parse command-line arguments
arguments <- docopt(doc, version = "0.1")

# Extract arguments
deseq2_res_dir <- arguments$deseq2_res_dir
comparison <- arguments$comparison
paired_grp <- ifelse(is.null(arguments$paired_grp), "No", arguments$paired_grp)
remove_batch_effect <- ifelse(!arguments$remove_batch_effect, "No", "Yes")
batch_col <- ifelse(is.null(arguments$batch_col), "None", arguments$batch_col)
covariate <- ifelse(is.null(arguments$covariate), "None", arguments$covariate)
bib_file <- tools::file_path_sans_ext(arguments$bib_file) # Extract base name without extension
template_file <- arguments$template_file

# Function to escape LaTeX special characters
escape_underscores <- function(instring) {
  gsub("_", "\\\\_", instring) # Escape underscores
}


# Function to generate LaTeX image inclusion code, only if the image file exists
generate_image_code_if_exists <- function(deseq2_res_dir, filename, caption) {
  filepath <- file.path(deseq2_res_dir, filename)
  if (file.exists(filepath)) {
    list(path = escape_underscores(filepath), caption = caption) # Return a list with path and caption
  } else {
    NULL # Return NULL if file doesn't exist
  }
}

#' Generate the LaTeX report by populating the template with data
#'
#' @param resultsDir Output directory for the report.
#' @param deseq2_res_dir Input directory containing the DESeq2 results.
#' @param comparison The raw comparison string.
#' @param comparison_clean The comparison string with underscores escaped for LaTeX.
#' @param paired  "Yes" or "No" indicating if a paired analysis was done.
#' @param covariates  String describing the covariates used, if any.
#' @param batchEffects  String describing the batch effects handling.
#' @param bib_file Base name of the BibTeX file.
#' @param template_file  Path to the Mustache template.
#'
#' @return The path to the generated LaTeX report file.

generate_latex_report <- function(resultsDir, deseq2_res_dir, comparison, comparison_clean, paired, covariates, batchEffects, batch_col, bib_file, template_file) {
  reportFile <- file.path(resultsDir, paste0("DESeq2_Analysis_Report_", comparison, ".tex"))

  comparison <- gsub("_vs_", "_", comparison) # Clean up comparison string
  has_deg = if (file.exists(file.path(deseq2_res_dir, paste0(comparison, "_IHWsigGenes.tsv")))) { TRUE } else { FALSE }
  no_deg = !has_deg

  # Create a list of data to populate the template.  `figures` is now a list of lists,
  # to accommodate an arbitrary number of figures. Each sublist represents a row in the report,
  # containing one or two figures. If has_pair is true, then a second figure will be generated.
  data <- list(
    comparison_clean = escape_underscores(comparison_clean),
    paired = escape_underscores(paired),
    covariates = escape_underscores(covariates),
    batchEffects = batchEffects,
    batch_col = escape_underscores(batch_col),
    has_deg = has_deg,
    no_deg = no_deg,
    pca_plots = list(
      list( # First row of plots
        image_path = escape_underscores(file.path(deseq2_res_dir, paste0(comparison, "_PCA.png"))),
        caption = "PCA plot of samples.",
        has_pair = ifelse(batchEffects == "Yes", TRUE, FALSE),
        pair_path = ifelse(batchEffects == "Yes", escape_underscores(file.path(deseq2_res_dir, paste0(comparison, "_PCA_after_batch_effect_correction.png"))), ""),
        pair_caption = ifelse(batchEffects == "Yes", "PCA after batch effect correction.", "")
      )
    ),
    expr_plots = if (file.exists(file.path(deseq2_res_dir, paste0(comparison, "_biotype_counts.png")))) {
      list(
          list( # Second row of plots
            image_path = escape_underscores(file.path(deseq2_res_dir, paste0(comparison, "_volcano.png"))),
            caption = "Volcano plot of DE genes: p-value",
            has_pair = TRUE,
            pair_path = escape_underscores(file.path(deseq2_res_dir, paste0(comparison, "_volcano_padj.png"))),
            pair_caption = "Volcano plot of DE genes: p-adjusted"
          ),
          list( # Third row of plots
            image_path = escape_underscores(file.path(deseq2_res_dir, paste0(comparison, "_biotype_counts.png"))),
            caption = "Biotype counts plot.",
            has_pair = TRUE,
            pair_path = escape_underscores(file.path(deseq2_res_dir, paste0(comparison, "_number_of_detected_genes.png"))),
            pair_caption = "Number of detected genes."
          )
        )
    } else {
      NULL
    },
    # Remaining figures (ORA and GSEA) in two separate lists
    figures_ORA = lapply(c("KEGG", "Reactome", "WikiPathway", "GO_BP", "GO_MF", "GO_CC"), function(x) {
      generate_image_code_if_exists(deseq2_res_dir, paste0(comparison, "_ORA_", x, "_dotplot.png"), paste("Over Representation Analysis of", escape_underscores(x)))
    }),
    figures_GSEA = lapply(c("Reactome", "WikiPathway", "GO_BP", "GO_MF", "GO_CC"), function(x) {
      generate_image_code_if_exists(deseq2_res_dir, paste0(comparison, "_GSEA_", x, "_dotplot.png"), paste("Gene Set Enrichment Analysis of", escape_underscores(x)))
    }),
    # BibTeX file
    bib_file = bib_file
  )

  # Render the template using whisker
  tex_content <- whisker::whisker.render(template = readLines(template_file), data = data)

  # Write the LaTeX report to file
  writeLines(tex_content, con = reportFile)

  return(reportFile)
}


resultsDir <- "./"
comparison_clean <- gsub("_vs_", " vs ", comparison)


reportFile <- generate_latex_report(resultsDir, deseq2_res_dir, comparison, comparison_clean, paired_grp, covariate, remove_batch_effect, batch_col, bib_file, template_file)
