library(tidyverse)

# -------------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------------
tsv_path <- "/data/projects/2021/MicrobialMetabolites/bacterial-supernatant/20_deseq2icbi/paired_grp/deseq2_11mix_vs_10mix/redo_deseq2_04032025/11mix_10mix_GSEA_GO_BP.tsv"

df <- read.delim(tsv_path, sep = "\t", stringsAsFactors = FALSE)
colnames(df) <- trimws(colnames(df))

# Standardize score column name if named 'enrichmentScore' or 'NES'
score_col <- if ("NES" %in% colnames(df)) "NES" else "enrichmentScore"

# -------------------------------------------------------------
# 2. Select Activated Terms (Neutrophil + Targets)
# -------------------------------------------------------------
target_terms <- c(
  "positive regulation of cell migration",
  "G protein-coupled receptor signaling pathway",
  "cilium organization",
  "chemotaxis",
  "leukocyte migration",
  "regulation of chemotaxis",
  "leukocyte chemotaxis"
)

df_activated <- df %>%
  filter(
    p.adjust < 0.1,
    .data[[score_col]] > 0,
    grepl("neutrophil", Description, ignore.case = TRUE) | Description %in% target_terms
  ) %>%
  distinct(Description, .keep_all = TRUE) %>%
  mutate(Direction = "Activated")

# -------------------------------------------------------------
# 3. Select Top 11 Suppressed Terms (Negative score, Highest setSize)
# -------------------------------------------------------------
df_suppressed <- df %>%
  filter(
    p.adjust < 0.1,
    .data[[score_col]] < 0
  ) %>%
  arrange(desc(setSize)) %>%
  slice_head(n = 10) %>%
  distinct(Description, .keep_all = TRUE) %>%
  mutate(Direction = "Suppressed")

# -------------------------------------------------------------
# 4. Combine and Format for Side-by-Side Plotting
# -------------------------------------------------------------
df_plot <- bind_rows(df_activated, df_suppressed) %>%
  mutate(
    Direction = factor(Direction, levels = c("Activated", "Suppressed")),
    score = .data[[score_col]]
  ) %>%
  group_by(Direction) %>%
  arrange(desc(p.adjust)) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# -------------------------------------------------------------
# 5. Create Plot (Horizontal Facets)
# -------------------------------------------------------------
p <- ggplot(df_plot, aes(x = p.adjust, y = Description)) +
  geom_point(
    aes(size = setSize, fill = score),
    shape = 21,
    color = "black",
    stroke = 0.8
  ) +
  
  # Horizontal side-by-side layout: columns by Direction
  facet_grid(. ~ Direction, scales = "free", space = "free_y") +
  
  # Reverse X-axis so smaller p-values are placed to the right
  scale_x_reverse(
    labels = scales::number_format(accuracy = 0.001),
    expand = expansion(mult = c(0.04, 0.1))
  ) +
  
  # Divergent fill gradient for points
  scale_fill_gradient2(
    low = "#5533ff",
    mid = "white",
    high = "#ff1a1a",
    midpoint = 0,
    name = "Enrichment Score",
    guide = guide_colorbar(
      barheight = unit(8, "cm"),
      barwidth  = unit(1.8, "cm"),
      title.vjust = 1.2
    )
  ) +
  
  # Dot sizing
  scale_size_continuous(
    range = c(5, 14),
    name = "Set size",
    guide = guide_legend(
      override.aes = list(shape = 21, fill = "gray70", color = "black")
    )
  ) +
  
  labs(
    title = "GO BP",
    x = "Adjusted p-value",
    y = NULL
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,  size = 32),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.grid.major = element_line(color = "#ebebeb"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#d9d9d9", color = "black"),
    strip.text.x = element_text( size = 32, color = "black"),
    axis.text.y = element_text(color = "black", size = 32),
    axis.text.x = element_text(color = "black", size = 30),
    axis.title.x = element_text(color = "black", size = 30, margin = margin(t = 15)),
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20, color = "black"),
    legend.key   = element_blank(),
    legend.spacing.y = unit(1.2, "cm")
  )


print(p)

out_dir <- "/data/scratch/kvalem/projects/2021/prj-honda_og/01_bacterial_products/001_organoid_bacterial_supernatant/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = file.path(out_dir, "neutrophil_GO_BP.svg"),
  plot = p,
  width = 30,
  height = 12
)

ggsave(
  filename = file.path(out_dir, "neutrophil_GO_BP.png"),
  plot = p,
  width = 30,
  height = 12,
  dpi = 300
)
