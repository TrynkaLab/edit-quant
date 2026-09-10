// One-time reference prep: samtools faidx, bwa index, and (optionally) a
// gene/exon BED from a GTF -- the Nextflow equivalent of setup_reference.sh.
// main.nf only invokes these when the corresponding output isn't already
// present on disk, matching setup_reference.sh's idempotent skip-if-exists
// behavior. Auto-built artifacts are always written under params.outdir
// rather than next to ref_fasta, so a read-only/shared reference directory
// is never written to.

process FAIDX {
    tag "${ref_fasta.name}"
    label 'tiny'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "${ref_fasta.name}.fai", emit: fai

    script:
    """
    samtools faidx --fai-idx ${ref_fasta.name}.fai ${ref_fasta}
    """
}

process BWA_INDEX {
    tag "${ref_fasta.name}"
    label 'medium'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "${ref_fasta.name}.{amb,ann,bwt,pac,sa}", emit: index_files

    script:
    """
    bwa index -p ${ref_fasta.name} ${ref_fasta}
    """
}

process GTF_TO_GENEBED {
    tag "${genes_gtf.name}"
    label 'small'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path genes_gtf

    output:
    path 'genome.genes.bed', emit: gene_bed
    path 'genome.exons.bed', emit: exon_bed

    script:
    """
    gtf_to_genebed.py ${genes_gtf} genome.genes.bed genome.exons.bed
    """
}
