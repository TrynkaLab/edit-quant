#!/usr/bin/env nextflow

include { validateAmpliconManifestContent; validateFastqManifestContent; formatValidationErrors } from './validateInputs.nf'

// Sanity check the run parameters, erroring out early on invalid combinations
def checkParamsMain(params) {

    // Input manifests
    checkParamsManifests(params)

    // Reference / annotation
    checkParamsReference(params)

    // Summary report
    checkParamsReport(params)
}


def checkParamsReport(params) {
        try {
            java.util.regex.Pattern.compile(params.report_wt_pattern.toString())
        } catch (java.util.regex.PatternSyntaxException e) {
            error("report_wt_pattern '${params.report_wt_pattern}' is not a valid regex - ${e.description}")
        }
}


def checkParamsManifests(params) {
        //------------------------------------------------------------
        // Defaults & sanity checks
        //------------------------------------------------------------
        if (params.amplicon_manifest == null) {
            error "amplicon_manifest file parameter is required: --amplicon_manifest"
        }

        if (!file(params.amplicon_manifest).exists()) {
            error("Specified amplicon manifest does not exist or is not readable: ${params.amplicon_manifest}")
        }

        def amplicon = validateAmpliconManifestContent(params.amplicon_manifest.toString())
        if (amplicon.errors) {
            error(formatValidationErrors("amplicon_manifest", amplicon.errors))
        }
        // rows that can't produce an amplicon are reported, not fatal -- see
        // validateAmpliconManifestContent
        amplicon.warnings.each { warning -> log.warn("amplicon_manifest: ${warning}") }

        if (params.fastq_manifest == null) {
            error "fastq_manifest file parameter is required: --fastq_manifest"
        }

        if (!file(params.fastq_manifest).exists()) {
            error("Specified fastq manifest does not exist or is not readable: ${params.fastq_manifest}")
        }

        def fastqErrors = validateFastqManifestContent(params.fastq_manifest.toString())
        if (fastqErrors) {
            error(formatValidationErrors("fastq_manifest", fastqErrors))
        }
}


def checkParamsReference(params) {
        if (params.ref_fasta == null) {
            error "ref_fasta file parameter is required: --ref_fasta"
        }

        if (params.gene_bed == null && params.genes_gtf == null) {
            error "One of gene_bed or genes_gtf is required: --gene_bed / --genes_gtf"
        }

        // gene_bed takes precedence in main.nf - warn so it's not mistaken for both being used
        if (params.gene_bed != null && params.genes_gtf != null) {
            log.warn "Both gene_bed and genes_gtf are set - gene_bed is used and genes_gtf is ignored"
        }

        if (params.gene_bed == null && params.exon_bed != null) {
            log.warn "exon_bed is set without gene_bed - exon_bed is ignored and both BEDs are built from genes_gtf"
        }
}
