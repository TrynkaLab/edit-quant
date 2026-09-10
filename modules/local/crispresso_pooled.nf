// Runs CRISPRessoPooled once per fastq-manifest sample, against the single
// shared deduplicated amplicons file produced by DEDUP_AMPLICONS. There is no
// per-sample linking to a specific amplicon manifest row -- every sample is
// screened against the full designed amplicon panel, matching real usage.

process CRISPRESSO_POOLED {
    tag "${sample}"
    label 'medium'
    publishDir "${params.outdir}/crispresso", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)
    path amplicons_unique
    val n_processes

    output:
    path "CRISPRessoPooled_on_${sample}", emit: results

    script:
    """
    CRISPRessoPooled \\
        --fastq_r1 ${fastq_1} \\
        --fastq_r2 ${fastq_2} \\
        --amplicons_file ${amplicons_unique} \\
        --n_processes ${n_processes} \\
        --output_folder . \\
        --name ${sample}
    """
}
