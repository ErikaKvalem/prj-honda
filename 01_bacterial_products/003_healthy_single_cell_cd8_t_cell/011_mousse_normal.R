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

data <- read_h5ad("/data/projects/2021/MicrobialMetabolites/single-cell-sorted-cd8/results/40_gex_surface_prot/003_trajectory_annotated_mudata_normal.h5ad")


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
  dplyr::filter(abs(avg_log2FC) > 0.5) %>%
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




sample_order <- c(

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
  subtitle = "p.adj < 0.01 & avg_log2FC > 0.5"
) + 
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme(
    axis.text.y = element_text(size = 40),
    axis.text.x = element_text(size = 40), angle=0,  plot.subtitle = element_text(size = 20),
    axis.text.x.top = element_text(size = 40, angle = 0, hjust = 1),
    
    legend.position = "right"
  )

p <- DoHeatmap(
  pbm,
  features = unique(topmarkers$gene),
  draw.lines = FALSE
) + labs(
  subtitle = "p.adj < 0.01 & avg_log2FC > 0.5"
) + 
  scale_fill_gradient2(
    name = "Z-score",   # <-- This changes the legend title cleanly
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme(
    axis.text.y = element_text(size = 40),
    axis.text.x = element_text(size = 40, angle = 0), # Fixed angle placement here
    plot.subtitle = element_text(size = 20),
    axis.text.x.top = element_text(size = 40, angle = 0, hjust = 1),
    legend.position = "right"
  )

print(p)

ggsave(
  "/data/projects/2021/MicrobialMetabolites/single-cell-sorted-cd8/results/prj_honda_15062026/publication_figures/healthy/heatmap_mousse_cytokines_normal.svg",
  plot = p,
  width = 15,
  height = 20,
  dpi = 300
)
