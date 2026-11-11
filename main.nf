#!/usr/bin/env nextflow
// Amplicon mapping + CRISPRessoPooled editing-efficiency pipeline.
//
// Translated from OTAR3086/analysis/mapping_amplicons/scripts: aligns each
// amplicon-manifest row's guides/primers to a reference genome with BWA,
// extracts the reference amplicon sequence, builds validation/QC reports and
// a CRISPRessoPooled amplicons file, then runs CRISPRessoPooled per fastq
// sample against the (deduplicated) amplicons of that sample's project_name.

include { validateParameters }              from 'plugin/nf-schema'
include { checkParamsMain }                  from './lib/checks.nf'
include { slugifySample; asBool }            from './lib/utils.nf'
include { FAIDX; BWA_INDEX; GTF_TO_GENEBED } from './modules/local/prepare_reference'
include { TSV_TO_FASTA }                     from './modules/local/tsv_to_fasta'
include { BWA_ALN as BWA_ALN_EXACT }         from './modules/local/bwa_align'
include { BWA_ALN as BWA_ALN_MISMATCH }      from './modules/local/bwa_align'
include { BUILD_REPORT }                     from './modules/local/build_report'
include { REASSIGNED_BAM }                   from './modules/local/reassigned_bam'
include { CHECK_SAMPLES }                    from './modules/local/check_samples'
include { CRISPRESSO_POOLED }                from './modules/local/crispresso_pooled'
include { RENDER_REPORT }                    from './modules/local/render_report'

workflow {

    // ------------------------------------------------------------
    // Check parameters
    if (asBool(params.skip_checks)) {
        log.warn("Skipping parameter checks")
    } else {

        // Validate parameters plugin
        validateParameters()

        // Custom parameter checks
        checkParamsMain(params)
    }

    ch_amplicon_manifest = channel.fromPath(params.amplicon_manifest, checkIfExists: true)
    ch_ref_fasta         = channel.fromPath(params.ref_fasta, checkIfExists: true)

    // --- .fai: build only if missing, matching setup_reference.sh's idempotency ---
    def ref_fai_path = params.ref_fai ?: "${params.ref_fasta}.fai"
    if (file(ref_fai_path).exists()) {
        ch_ref_fai = channel.fromPath(ref_fai_path, checkIfExists: true)
    } else {
        FAIDX(ch_ref_fasta)
        ch_ref_fai = FAIDX.out.fai
    }

    // --- BWA index: build only if missing ---
    def bwa_prefix_path    = params.bwa_index_prefix ?: params.ref_fasta
    def bwa_index_complete = ['amb', 'ann', 'bwt', 'pac', 'sa'].every { ext ->
        file("${bwa_prefix_path}.${ext}").exists()
    }
    def bwa_index_prefix_name
    if (bwa_index_complete) {
        ch_bwa_index_files = channel.fromPath(
            ['amb', 'ann', 'bwt', 'pac', 'sa'].collect { ext -> "${bwa_prefix_path}.${ext}" },
            checkIfExists: true
        ).toList()
        bwa_index_prefix_name = file(bwa_prefix_path).name
    } else {
        BWA_INDEX(ch_ref_fasta)
        ch_bwa_index_files    = BWA_INDEX.out.index_files
        bwa_index_prefix_name = file(params.ref_fasta).name
    }

    // --- gene / exon BED: use given paths, else build from a GTF ---
    if (params.gene_bed) {
        ch_gene_bed = channel.fromPath(params.gene_bed, checkIfExists: true)
        ch_exon_bed = params.exon_bed
            ? channel.fromPath(params.exon_bed, checkIfExists: true)
            : channel.fromPath("${projectDir}/assets/NO_FILE.bed")
    } else {
        GTF_TO_GENEBED(channel.fromPath(params.genes_gtf, checkIfExists: true))
        ch_gene_bed = GTF_TO_GENEBED.out.gene_bed
        ch_exon_bed = GTF_TO_GENEBED.out.exon_bed
    }

    // --- amplicon manifest -> query FASTA (guides + primers) ---
    TSV_TO_FASTA(ch_amplicon_manifest)

    // --- two BWA passes: exact-match (source of truth) + mismatch-tolerant (off-targets) ---
    ch_exact_input = TSV_TO_FASTA.out.fasta.map { fasta -> tuple('queries.exact', fasta) }
    ch_mm_input    = TSV_TO_FASTA.out.fasta.map { fasta -> tuple('queries', fasta) }

    exact_bwa = BWA_ALN_EXACT(ch_exact_input, ch_bwa_index_files, bwa_index_prefix_name, 0, 1000)
    mm_bwa    = BWA_ALN_MISMATCH(ch_mm_input, ch_bwa_index_files, bwa_index_prefix_name, params.bwa_n, params.max_xa)

    // --- build all reports + extract amplicons + write crispresso_amplicons.txt ---
    build_report_result = BUILD_REPORT(
        ch_amplicon_manifest,
        exact_bwa.sam,
        mm_bwa.sam,
        ch_ref_fasta,
        ch_ref_fai,
        ch_gene_bed,
        ch_exon_bed,
        params.gene_window,
        params.max_mismatch,
        params.high_copy_threshold,
        params.reassign_primary,
        params.guide_pam,
        params.guide_pam_side,
        params.quantification_window_center,
    )

    // runs 0 or 1 times: reassigned_sam is only emitted when reassign_primary
    // corrected at least one position
    REASSIGNED_BAM(build_report_result.reassigned_sam)

    // --- split amplicons by project_name (deduplicated by sequence, as
    //     CRISPRessoPooled requires) and match each fastq sample to its project ---
    CHECK_SAMPLES(
        channel.fromPath(params.fastq_manifest, checkIfExists: true),
        ch_amplicon_manifest,
        build_report_result.crispresso_amplicons,
    )

    ch_project_amplicons = CHECK_SAMPLES.out.amplicons
        .flatten()
        .map { amplicons -> tuple(amplicons.name, amplicons) }

    // Samples without amplicons for their project are skipped with a warning,
    // rather than screened against another project's panel
    def skip_reasons = [
        NO_PROJECT_IN_AMPLICON_MANIFEST: "its project_name has no rows in the amplicon manifest",
        NO_AMPLICONS_EXTRACTED         : "none of its project's amplicons could be extracted from the reference or supplied in the manifest",
    ]

    ch_fastq = CHECK_SAMPLES.out.sample_check
        .splitCsv(header: true, sep: '\t')
        .filter { row ->
            if (row.status != 'OK') {
                log.warn("Skipping sample '${row.sample}' (project_name=${row.project_name}): ${skip_reasons[row.status] ?: row.status}")
                return false
            }
            return true
        }
        .map { row ->
            tuple(row.amplicons_file, slugifySample(row.sample), file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }
        .combine(ch_project_amplicons, by: 0)
        .map { _amplicons_name, sample, fastq_1, fastq_2, amplicons -> tuple(sample, fastq_1, fastq_2, amplicons) }

    // --- CRISPRessoPooled per fastq sample, against its project's amplicons ---
    // cores/memory come from the CRISPRESSO_POOLED block in conf/processes.config
    // (driven by crispresso_n_processes / crispresso_memory); the process passes
    // task.cpus to --n_processes so the request and the flag can't drift
    CRISPRESSO_POOLED(
        ch_fastq,
        params.crispresso_min_reads_to_use_region,
        params.crispresso_skip_failed,
        params.crispresso_suppress_plots,
    )

    // --- summary report over all samples ---
    RENDER_REPORT(
        CRISPRESSO_POOLED.out.results.collect().ifEmpty([]),
        CHECK_SAMPLES.out.sample_check,
        CHECK_SAMPLES.out.amplicon_replicates,
        build_report_result.amplicon_report,
        build_report_result.guide_qc,
        ch_amplicon_manifest,
        file("${projectDir}/assets/report.html.j2"),
        workflow.manifest.version ?: 'unknown',
    )
}
