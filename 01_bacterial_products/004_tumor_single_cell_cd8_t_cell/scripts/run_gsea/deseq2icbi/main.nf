#!/usr/bin/env nextflow

/**
 * Workflow for performing differential expression analysis using DESeq2 and generating a report.
 */

// Process to copy raw counts to the results directory.
process add_raw_counts {
    publishDir "${params.resultDir}", mode: "link" // Publish results by linking, avoids copying large files

    input:
      path raw_counts // Input raw counts file

    output:
        file "raw_counts.tsv" // Output raw counts file (renamed for consistency)

    """
    cp -L ${raw_counts} raw_counts.tsv
    """
}

// Process to perform differential expression analysis using DESeq2.
process de_analysis {

  errorStrategy "ignore" // Ignore errors in this process (for individual comparisons)
  publishDir "${params.resultDir}", mode: "link" // Publish results by linking

  cpus = params.cpus // Use the specified number of CPUs

  input:
    path sample_sheet // Input sample sheet file
    path raw_counts   // Input raw counts file
    path goi_f        // Input file with genes of interest (optional)
    path gtf          // Input GTF file
    each contrast    // Each contrast definition (from params.contrasts)

  output:
    path out_dir, emit: deseq2_res_dir // Output directory containing DESeq2 results
    val contrast, emit: contrast       // Emit the contrast value for downstream use
    val out_dir, emit: deseq2_subdir   // Emit the subdirectory name

  script:
    def test = contrast.TEST       // Extract test condition from contrast
    def ref = contrast.REF        // Extract reference condition from contrast
    def condcol = contrast.CONDCOL != null ? contrast.CONDCOL : "group" //
    def goi = goi_f.name != 'NO_FILE' ? "--genes_of_interest=$goi_f" : '' // Include genes of interest if provided
    out_dir = (contrast.CONDCOL) ? contrast.CONDCOL + "/" + test + "_vs_" + ref : test + "_vs_" + ref   // Construct output directory name

    def remove_batch_effect = params.remove_batch_effect ? "--remove_batch_effect" : '' // Include batch effect removal if specified
    def batch_col = params.batch_col != "" ? "--batch_col=${params.batch_col}" : '' // Include batch column if specified
     def sample_col = (params.sample_col && params.sample_col != "") ? "--sample_col=${params.sample_col}" : '' // Set sample column 
    def paired_grp = params.paired_grp != "" ? "--paired_grp=${params.paired_grp}" : '' // Include paired group info if specified
    def covariate_formula = params.covariate_formula != "" ? '--covariate_formula="' + params.covariate_formula + '"' : '' // Include covariate formula if specified 
    def ext_opt = params.ext_opt != "" ? params.ext_opt : "" // Include any extra options

    """
    mkdir -p ${out_dir}
    /data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/run_gsea/deseq2icbi/runDESeq2_ICBI.R ${sample_sheet} ${raw_counts} \\
      --result_dir=${out_dir} \\
      --c1=${test} \\
      --c2=${ref} \\
      --condition_col=${condcol} \\
      ${sample_col} \\
      ${remove_batch_effect} \\
      ${batch_col} \\
      ${paired_grp} \\
      ${covariate_formula} \\
      ${ext_opt} \\
      --plot_title="${test} vs ${ref}" \\
      --gtf_file=${gtf} \\
      --organism=${params.organism} \\
      ${goi} \\
      --n_cpus=${task.cpus} \\
      --save_workspace \\
      --save_sessioninfo
    """
}


// Process to generate a PDF report from the DESeq2 results.
process generate_report {

    publishDir "${params.resultDir}/${subdir}", mode: "copy" // Publish the report
    cache false // Disable caching for this process

    input:
      path analysis        // Path to the DESeq2 results directory
      val contrast         // The contrast used for this analysis
      val subdir
      path bib_file       // Path to the BibTeX file
      path template_file   // Path to the report template

    output:
       path "DESeq2_Analysis_Report_${contrast_name}.pdf", emit: pdf_file // Output PDF report

    script:
      contrast_name = contrast.TEST + "_vs_" + contrast.REF
      paired_grp = params.paired_grp != "" ? "--paired_grp=${params.paired_grp}" : ""
      remove_batch_effect = params.remove_batch_effect ? "--remove_batch_effect" : ""
      batch_col = params.batch_col != "" ? "--batch_col=${params.batch_col}" : ''
      covariate_formula = params.covariate_formula != "" ? '--covariate="' + params.covariate_formula + '"' : ''

    """
    module purge
    module load texlive
    /data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/28012026/run_gsea/deseq2icbi/Report.R ${analysis} ${contrast_name} ${bib_file}  ${template_file} \\
      ${paired_grp} \\
      ${remove_batch_effect} \\
      ${batch_col} \\
      ${covariate_formula}

    pdflatex -interaction=nonstopmode DESeq2_Analysis_Report_${contrast_name}.tex
    bibtex DESeq2_Analysis_Report_${contrast_name}
    pdflatex -interaction=nonstopmode DESeq2_Analysis_Report_${contrast_name}.tex
    pdflatex -interaction=nonstopmode DESeq2_Analysis_Report_${contrast_name}.tex
    """
}

// Process to generate a BibTeX file citing used R packages.
process add_cite {
    publishDir "${params.resultDir}", mode: "link"

    input:
      val packages // List of R packages to cite

    output:
        path "grateful-refs.bib", emit: bib_file  // Output BibTeX file
        path "grateful-report.html"              // Output HTML report (optional)


    script:

    """
    #!/usr/bin/env Rscript
    library("grateful")
    cite_packages(cite.tidyverse = FALSE, pkgs = ${packages}, out.dir = "./")
    """
}


// Main workflow definition.
workflow {
    sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists:true) // Input sample sheet
    raw_counts = Channel.fromPath(params.raw_counts, checkIfExists:true)   // Input raw counts
    goi = file(params.goi, checkIfExists:false)                            // Input genes of interest (optional)
    gtf = Channel.fromPath(params.gtf, checkIfExists:true)                  // Input GTF file
    contrasts = Channel.fromList(params.contrasts)                        // Input contrast definitions
    packages_to_cite = "c(${params.cite_packages.collect { "\"$it\"" }.join(', ')})" // Format packages for R
    report_template_file = file(workflow.projectDir + "/report_template.tex", checkIfExists:true) // Path to report template

    add_raw_counts(raw_counts)                  // Run the `add_raw_counts` process
    deseq2_out = de_analysis(sample_sheet, raw_counts, goi, gtf, contrasts)  // Run the `de_analysis` process
    citation_out = add_cite(packages_to_cite)            // Run the `add_cite` process
    report_out = generate_report(deseq2_out.deseq2_res_dir, deseq2_out.contrast, deseq2_out.deseq2_subdir, citation_out.bib_file, report_template_file) // Run the `generate_report` process
}
