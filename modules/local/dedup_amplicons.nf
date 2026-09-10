// CRISPRessoPooled requires unique amplicon sequences. Deduplicates
// crispresso_amplicons.txt by amplicon_seq (column 2), keeping one row per
// unique sequence -- replicates the user's real analysis script exactly:
//   sort -t$'\t' -k2,2 -u crispresso_amplicons.txt > amplicons__unique.txt

process DEDUP_AMPLICONS {
    tag "${crispresso_amplicons.name}"
    label 'tiny'
    publishDir "${params.outdir}/crispresso", mode: 'copy'

    input:
    path crispresso_amplicons

    output:
    path 'amplicons_unique.txt', emit: amplicons_unique

    script:
    """
    sort -t\$'\\t' -k2,2 -u ${crispresso_amplicons} > amplicons_unique.txt
    """
}
