// Single BWA aln + samse process, invoked twice from main.nf:
//   prefix='queries.exact' n_mismatches=0            -> source of truth for position/gene-overlap
//   prefix='queries'       n_mismatches=params.bwa_n  -> off-target search (mismatch-tolerant)
// Filenames match the original run_pipeline.sh outputs exactly
// (align/queries.exact.sorted.bam vs align/queries.sorted.bam).

process BWA_ALN {
    tag "${prefix}"
    label 'normal'
    publishDir "${params.outdir}/align", mode: 'copy'

    input:
    tuple val(prefix), path(fasta)
    path index_files
    val index_prefix_name
    val n_mismatches
    val max_xa

    output:
    tuple val(prefix), path("${prefix}.sam"), emit: sam
    tuple val(prefix), path("${prefix}.sorted.bam"), path("${prefix}.sorted.bam.bai"), emit: bam

    script:
    """
    bwa aln -t ${task.cpus} -n ${n_mismatches} -N ${index_prefix_name} ${fasta} \\
        > ${prefix}.sai 2> bwa_aln.${prefix}.log
    bwa samse -n ${max_xa} ${index_prefix_name} ${prefix}.sai ${fasta} \\
        > ${prefix}.sam 2> bwa_samse.${prefix}.log
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${prefix}.sam
    samtools index ${prefix}.sorted.bam
    """
}
