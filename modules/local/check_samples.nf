// Matches each fastq-manifest sample to the amplicons designed for its
// project_name (bin/check_samples.py). Emits one amplicons file per project,
// with technical replicates (wells sharing an amplicon) merged, plus
// reports/sample_check.tsv and reports/amplicon_replicates.tsv; main.nf warns about and skips any
// sample whose project has no amplicons, rather than screening it against
// another project's panel.

process CHECK_SAMPLES {
    tag "${fastq_manifest.name}"
    label 'tiny'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'reports/**'
    publishDir "${params.outdir}/crispresso", mode: 'copy', pattern: 'amplicons_by_project/**'

    input:
    path fastq_manifest
    path amplicon_manifest
    path crispresso_amplicons

    output:
    path 'reports/sample_check.tsv', emit: sample_check
    path 'reports/amplicon_replicates.tsv', emit: amplicon_replicates
    path 'amplicons_by_project/*.amplicons.txt', optional: true, emit: amplicons

    script:
    """
    check_samples.py \\
        --amplicon-manifest ${amplicon_manifest} \\
        --crispresso-amplicons ${crispresso_amplicons} \\
        --fastq-manifest ${fastq_manifest} \\
        --outdir .
    """
}
