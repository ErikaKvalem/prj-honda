# R script to run for "Bioinformatic Service Unit" projects

## Installation

### Get the code

To install the DESeq2ICBI pipeline clone the script as submodule into your data
analysis dir.

```
git submodule add git@gitlab.i-med.ac.at:icbi-lab/tools/deseq2icbi.git deseq2icbi
```

### Deploy

Run the `deploy.sh` script from the `deseq2icbi` directory

```
cd deseq2icbi
./deploy.sh
```

### Configure

Edit the `nextflow.config` file in the data analysis dir to configure your comparisons and
other settings.

## Usage

After editing the `nextflow.config` file, run the nextflow pipeline from within the data analysis dir as follows:

```
bash run_deseq_nf.sh 
```

### R Script options
```
runDESeq2_ICBI.R

Usage:
  runDESeq2_ICBI.R <sample_sheet> <count_table> --result_dir=<res_dir> --c1=<c1> --c2=<c2> [options]
  runDESeq2_ICBI.R --help

Arguments:
  <sample_sheet>                CSV file with the sample annotations.
  <count_table>                 TSV file with the read counts

Mandatory options:
  --result_dir=<res_dir>        Output directory
  --c1=<c1>                     Contrast level 1 (perturbation). Needs to be contained in condition_col.
  --c2=<c2>                     Contrast level 2 (baseline). Needs to be contained in condition_col.

Optional options:
  --nfcore                      Indicate that the input samplesheet is from the nf-core RNA-seq ppipeline.
                                Will merge entries from the same sample and infer the sample_id from `group` and `replicate`.
                                If set, this option overrides `sample_col`.
  --condition_col=<cond_col>    Column in sample annotation that contains the condition [default: group]
  --sample_col=<sample_col>     Column in sample annotation that contains the sample names
                                (needs to match the colnames of the count table). [default: sample]
  --paired_grp=<paired_grp>     Column that conatins the name of the paired samples, when dealing with
                                paired data.
  --remove_batch_effect         Indicate that batch effect correction should be applied [default: FALSE]
                                If batch effect correction should be performed, a batch column is needed in the
                                samplesheet (see also --batch_col below)
  --batch_col=<batch_col>       Optional: column in sample annotation that contains the batch
  --covariate_formula=<formula> Formula to model additional covariates (need to be columns in the samplesheet)
                                that will be appended to the formula built from `condition_col`.
                                E.g. `+ age + sex`. Per default, no covariates are modelled.
  --plot_title=<title>          Title shown above plots. Is built from contrast per default.
  --prefix=<prefix>             Results file prefix. Is built from contrasts per default.
  --fdr_cutoff=<fdr>            False discovery rate for GO analysis and volcano plots [default: 0.1]
  --fc_cutoff=<log2 fc cutoff>  Fold change (log2) cutoff for volcano plots [default: 1]
  --gtf_file=<gtf>              Path to the GTF file used for featurecounts. If specified, a Biotype QC
                                will be performed.
  --gene_id_type=<id_type>      Type of the identifier in the `gene_id` column compatible with AnnotationDbi [default: ENSEMBL]
  --n_cpus=<n_cpus>             Number of cores to use for DESeq2 [default: 1]
  --skip_gsea                   Skip Gene-Set-Enrichment-Analysis step
  --genes_of_interest=<genes>   File (tsv) containing a list of genes to highlight in the volcano plot (column must be named
                                "gene_name").
                                If an optional column named "group" is present, each gene will be associated with the corresponding
                                gene group and a separate volcanon plot for each gene group will be generated (e.g. cytokines).
  --organism=<human|mouse>      Ensebml annotation db [default: human]
  --save_workspace              Save R workspace for this analysis [default: FALSE]
  --save_init_workspace         Save R workspace before analysis for manual step by step debugging [default: FALSE]
  --save_sessioninfo            Save R sessionInfo() to keep info about library version [default: TRUE]
  --skip_heatplots              Skip the generation of Heatplots for the Dataset [default: TRUE]
  --merge_duplicate_genes        Merge duplicate genes by summing up the counts. Additional ensembl IDs will be added as new column
                                "secondary_gene_ids" in the result tables [default: FALSE]
```

The R script takes the merged count table from the nf-core rnaseq pipeline as input (featurecounts.merged.counts.tsv, see example).

As second input file you need the sample table in csv format. See "sampleTableN.tsv" as example. We suggested to use the sample table that was used for running the nf-core rnaseq pipeline.

| group | replicate | pair_group | fastq_1 | fastq_2 | strandedness |
| ------| --------- | ---------- | ------- | ------- | ------------ |
| WT | 1 | | wt1_1.fq.gz | wt1_2.fq.gz | forward |
| WT | 2 | | wt2_1.fq.gz | wt2_2.fq.gz | forward |
| WT | 3 | | wt3_1.fq.gz | wt3_2.fq.gz | forward |
| KO | 1 | | ko1_1.fq.gz | ko1_2.fq.gz | forward |
| KO | 2 | | ko2_1.fq.gz | ko2_2.fq.gz | forward |
| KO | 3 | | ko3_1.fq.gz | ko3_2.fq.gz | forward |


## Note:
You need to run the 3.x version of the RNA-seq pipeline. The current release of
3.0 contains a bug, so use the dev version currently! 

## Methods summary

A summary of the methods and a result file description can be found in the 
`DESeq2_Analysis_Report_[TEST]_vs_[CTRL].pdf` in each `[TEST]_vs_[CTRL]`
comparison result folder. In this folder you will also find detailed
information about the R package version in use (**_sessionInfo.json`).
