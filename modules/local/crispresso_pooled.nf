// Runs CRISPRessoPooled once per fastq-manifest sample, against the
// deduplicated amplicons of that sample's project_name (from CHECK_SAMPLES).

include { asBool } from '../../lib/utils.nf'

process CRISPRESSO_POOLED {
    tag "${sample}"
    label 'medium'
    publishDir "${params.outdir}/crispresso", mode: 'copy'

    input:
        tuple val(sample), path(fastq_1), path(fastq_2), path(amplicons)
        val n_processes
        val min_reads_to_use_region
        val skip_failed

    output:
        path "CRISPRessoPooled_on_${sample}", emit: results
        path "*.html", emit: html_reports

    script:
    // CRISPResso's own --skip_failed: one amplicon that errors is dropped and the
    // rest of the sample's panel still runs, instead of the whole CRISPRessoPooled
    // task failing. Amplicons missing from the output are reported as not_tested
    // (missing_from_output) by RENDER_REPORT, so they are never silently lost.
    def skip_failed_arg = asBool(skip_failed) ? '--skip_failed' : ''
    """
    CRISPRessoPooled \\
        --fastq_r1 ${fastq_1} \\
        --fastq_r2 ${fastq_2} \\
        --amplicons_file ${amplicons} \\
        --n_processes ${n_processes} \\
        --min_reads_to_use_region ${min_reads_to_use_region} \\
        ${skip_failed_arg} \\
        --output_folder . \\
        --name ${sample}
    """
}
