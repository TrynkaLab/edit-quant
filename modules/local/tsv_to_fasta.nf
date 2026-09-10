// Converts the amplicon design manifest into a FASTA of guide/primer query
// sequences to align, via bin/tsv_to_fasta.py (unchanged from the original
// mapping_amplicons pipeline).

process TSV_TO_FASTA {
    tag "${amplicon_manifest.name}"
    label 'tiny'
    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:
    path amplicon_manifest

    output:
    path 'queries.fasta', emit: fasta

    script:
    """
    tsv_to_fasta.py ${amplicon_manifest} > queries.fasta
    """
}
