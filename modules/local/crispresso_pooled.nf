// Runs CRISPRessoPooled once per fastq-manifest sample, against the
// deduplicated amplicons of that sample's project_name (from CHECK_SAMPLES).

process CRISPRESSO_POOLED {
    tag "${sample}"
    label 'medium'
    publishDir "${params.outdir}/crispresso", mode: 'copy'

    input:
        tuple val(sample), path(fastq_1), path(fastq_2), path(amplicons)
        val n_processes
        val min_reads_to_use_region

    output:
        path "CRISPRessoPooled_on_${sample}", emit: results
        path "*.html", emit: html_reports

    script:
    """
    CRISPRessoPooled \\
        --fastq_r1 ${fastq_1} \\
        --fastq_r2 ${fastq_2} \\
        --amplicons_file ${amplicons} \\
        --n_processes ${n_processes} \\
        --min_reads_to_use_region ${min_reads_to_use_region} \\
        --output_folder . \\
        --name ${sample}
    """
}
