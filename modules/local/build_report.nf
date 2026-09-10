// Wraps bin/build_report.py (unchanged): turns the two BWA alignments into
// validation/amplicon/off-target reports, extracted amplicon FASTAs, and the
// CRISPRessoPooled amplicons description file. See mapping_amplicons/README.md
// for the full column/status-code documentation.

process BUILD_REPORT {
    tag "${amplicon_manifest.name}"
    label 'normal_plus'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'reports/**'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'amplicons/**'

    input:
    path amplicon_manifest
    tuple val(exact_label), path(exact_sam)
    tuple val(mm_label), path(mm_sam)
    path ref_fasta
    path ref_fai
    path gene_bed
    path exon_bed
    val gene_window
    val max_mismatch
    val high_copy_threshold
    val reassign_primary

    output:
        path 'reports/validation_report.tsv', emit: validation_report
        path 'reports/amplicon_report.tsv', emit: amplicon_report
        path 'reports/off_targets.tsv', emit: off_targets
        path 'reports/crispresso_amplicons.txt', emit: crispresso_amplicons
        path 'reports/reassigned_primaries.tsv', emit: reassigned_log
        path 'amplicons/*.fa', emit: amplicon_fastas
        path 'amplicons/all_amplicons.fasta', emit: all_amplicons_fasta
        path '*.reassigned.sam', optional: true, emit: reassigned_sam

    script:
    def exon_arg     = (exon_bed.name != 'NO_FILE.bed') ? "--exon-bed ${exon_bed}" : ''
    def reassign_arg = reassign_primary ? '--reassign-primary' : ''
    """
    build_report.py \\
        --tsv ${amplicon_manifest} \\
        --exact-sam ${exact_sam} \\
        --sam ${mm_sam} \\
        --ref ${ref_fasta} \\
        --fai ${ref_fai} \\
        --gene-bed ${gene_bed} \\
        --gene-window ${gene_window} \\
        --outdir . \\
        --max-mismatch ${max_mismatch} \\
        --high-copy-threshold ${high_copy_threshold} \\
        ${exon_arg} \\
        ${reassign_arg}
    """
}
