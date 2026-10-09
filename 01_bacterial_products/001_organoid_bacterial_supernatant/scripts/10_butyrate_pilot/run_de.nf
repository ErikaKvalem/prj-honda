#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process DESEQ2 {
    tag { meta.id }
    publishDir "../data/20_deseq2icbi", mode: "copy"

    cpus 4
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-de2"
    errorStrategy 'finish'

    input:
        tuple val(meta), path(expr), path(samplesheet)

    output:
        path("deseq2_res*")

    script:
    def covariates = meta['covariates'] ? "--covariate_formula=${meta['covariate_formula']}" : ""
    def sc_rm = meta['singlecell'] ? "rm deseq2_res_sc_**/*detectedGenes*_min_10_reads*.tsv" : ""
    def paired = meta['paired_grp'] ? "--paired_grp=${meta['paired_grp']}" : ""
    """
    mkdir deseq2_res_${meta.id}
    runDESeq2_ICBI.R $samplesheet $expr \\
        --result_dir=deseq2_res_${meta.id} \\
        --c1=${meta.c1} \\
        --c2=${meta.c2} \\
        --sample_col=${meta.sample_col} \\
        --condition_col=${meta.condition_col} \\
        --n_cpus=${task.cpus} \\
        ${paired} \\
        ${covariates}

    ${sc_rm}
    """
}

workflow {
    def input_path = "../../data/70_de_analysis/71_make_pseudobulk/"
    DESEQ2(
        Channel.from(
            [
                [id: "basal_vs_apical", c1: "apical", c2: "basal", sample_col: "sample", condition_col: "orientation", covariate_formula:"+organoid+treatment"],
                [id: "treatment_vs_control", c1: "10mM_NaB", c2: "ctrl", sample_col: "sample", condition_col: "treatment", covariate_formula:"+organoid+orientation"],
            ]
        ).map {
            it -> [
                it,
                file("../data/10_rnaseq_pipeline/star_salmon/salmon.merged.gene_counts.tsv", checkIfExists: true),
                file("../tables/samplesheet_rnaseq.csv", checkIfExists: true)
            ]
        }
    )
}
