#!/bin/bash
nextflow run main.nf -profile icbi --genes_of_interest="/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/run_gsea/deseq2icbi/genes_of_interest.txt" -resume "$@"
