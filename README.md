# prj-honda

![Version Badge](https://img.shields.io/badge/Version-1.0.2-brightgreen?style=for-the-badge)

## Introduction

Bulk RNA-seq analysis of mourine organoids and in vivo CD8+ sorted single-cell RNA seq analysis



01_bacterial_products/
 -  001_organoid_bacterial_supernatant:  Nextflow / nf-core configuration profiles and DESeq2 downstream analysis scripts
 -  003_healthy_single_cell_cd8_t_cell: Cellranger multi configuration files and preprocessing, MuData setup, annotation, TCR analysis
 -  004_tumor_single_cell_cd8_t_cell: Cellranger multi config files (GEX + V(D)J + ADT/CITE-seq), Multimodal integration, trajectory, scCODA, metabolism

02_molecular_mimicry: External submodule / pointer to mimicry prediction pipeline

## Authors

Erika Kvalem, Gregor Sturm
