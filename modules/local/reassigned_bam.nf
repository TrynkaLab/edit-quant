// Sorts + indexes the corrected SAM that BUILD_REPORT writes when
// reassign_primary is enabled and at least one position was corrected --
// mirrors run_pipeline.sh's final conditional block. Runs zero or one times:
// BUILD_REPORT's reassigned_sam output is optional, so this process simply
// never fires when there's nothing to sort.

process REASSIGNED_BAM {
    tag 'reassigned'
    label 'small'
    publishDir "${params.outdir}/align", mode: 'copy'

    input:
    path reassigned_sam

    output:
    path 'queries.exact.reassigned.sorted.bam', emit: bam
    path 'queries.exact.reassigned.sorted.bam.bai', emit: bai

    script:
    """
    samtools sort -@ ${task.cpus} -o queries.exact.reassigned.sorted.bam ${reassigned_sam}
    samtools index queries.exact.reassigned.sorted.bam
    """
}
