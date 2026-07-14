library(devtools)
#install_github("azkajavaid/MousseR-package")

library(mousseR)
#vignette("mousse-vignette", package = "mousseR")

#install.packages("remotes")
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)

library(Seurat)


library(Seurat)
library(reticulate)
library(anndata)
library(dplyr)
library(stringr)
library(tidyverse)
library(fastplyr)
library(Matrix)
library(stringr)
##############

data <- read_h5ad("/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/mdata_gex_12022026_paga_anno.h5ad")


#####
data_seurat <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 500, min.cells = 30)

DefaultAssay(data_seurat) <- "RNA"

data_seurat <- NormalizeData(data_seurat)
data_seurat <- FindVariableFeatures(data_seurat)
data_seurat <- ScaleData(data_seurat)

Idents(data_seurat) <- "sample_id"



######## Apply the mousse function

counts.matrix <- as.data.frame(t(as.matrix(data_seurat@assays$RNA@layers$data)))
rownames(counts.matrix) <- colnames(data_seurat)
colnames(counts.matrix) <- str_to_title(rownames(data_seurat))
mousse.scores <- mousse(counts.matrix = counts.matrix, numGenes = 60)
#steps below only necessary for the heatmap display formatting in pdf using Latex
colnames(mousse.scores) <- gsub("α", "a", colnames(mousse.scores))
colnames(mousse.scores) <- gsub("β", "b", colnames(mousse.scores))
colnames(mousse.scores) <- gsub("IFNγ", "IFNg", colnames(mousse.scores))

data_seurat <- AddMetaData(
  data_seurat,
  metadata = mousse.scores
)

data_seurat[["Mousse"]] <- CreateAssayObject(data = (t(mousse.scores)))
DefaultAssay(data_seurat) <- "Mousse"
data_seurat <- ScaleData(data_seurat)
Idents(data_seurat) <- "sample_id"
data_seurat[["Mousse"]] <- CreateAssayObject(
  counts = as.matrix(t(mousse.scores))
)


markers <- FindAllMarkers(
  data_seurat,
  assay = "Mousse",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0,
  only.pos = FALSE
)

topmarkers <- markers %>% 
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::filter(abs(avg_log2FC) > 0.3) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::ungroup()

pbm <- AggregateExpression(
  data_seurat,
  assays = "Mousse",
  group.by = c("sample_id"),
  return.seurat = TRUE
)




Idents(pbm) <- pbm@active.ident

rename_samples <- c(
  "GF-ICI1" = "GF_noICI1",
  "GF-ICI2" = "GF_noICI2",
  "GF-ICI1-plus" = "GF1",
  "GF-ICI2-plus" = "GF2",
  "g10mix-ICI1" = "ctrl1",
  "g10mix-ICI2" = "ctrl2",
  "g11mix-ICI1" = "effector1",
  "g11mix-ICI2" = "effector2"
)

pbm <- RenameCells(
  pbm,
  new.names = rename_samples[colnames(pbm)]
)

sample_order <- c(
  "GF_noICI1",
  "GF_noICI2",
  "GF1",
  "GF2",
  "ctrl1",
  "ctrl2",
  "effector1",
  "effector2"
)

# Reorder columns
pbm <- pbm[, sample_order]

# Replace metadata sample_id with renamed names
pbm$sample_id <- factor(
  colnames(pbm),
  levels = sample_order
)

# Set identities
Idents(pbm) <- "sample_id"

# Check
pbm@active.ident
levels(Idents(pbm))

p <- DoHeatmap(
  pbm,
  features = unique(topmarkers$gene),
  draw.lines = FALSE
) + labs(
  subtitle = "p.adj < 0.01 & avg_log2FC > 0.25"
) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 20), angle=0,  plot.subtitle = element_text(size = 12),
    axis.text.x.top = element_text(size = 20, angle = 0, hjust = 0.5),
    
    legend.position = "right"
  )

print(p)

#ggsave(
#  "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/figures/heatmap_mousse_cytokines.svg",
#  plot = p,
#  width = 15,
#  height = 20,
#  dpi = 300
#)
