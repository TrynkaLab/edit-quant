#!/usr/bin/env nextflow

// Content validation for the user-supplied TSV manifests (amplicon manifest,
// fastq manifest). Each validator returns a list of problems, so a malformed
// file fails fast with one aggregated error instead of falling over deep in a
// downstream process.

include { slugifySample; sanitizeName; normaliseHeader } from './utils.nf'

// The only accepted "not set" marker for optional fields.
def SENTINEL() { return "none" }

// Values people commonly mistype in place of the "none" sentinel - these are
// deliberately rejected everywhere rather than silently tolerated
def BLANKISH() { return ["", "na", "n/a", SENTINEL(), "null", "nan", "false"] }

def isBlankish(String value) {
    return value == null || BLANKISH().contains(value.trim().toLowerCase())
}

// value must be exactly the sentinel, or satisfy pattern - used for optional columns
def checkOptional(String value, def pattern, String colName, String hint) {
    if (value == SENTINEL()) return null
    if (isBlankish(value)) {
        return "${colName} '${value}' is invalid - use ${hint} or '${SENTINEL()}'"
    }
    if (pattern != null && !(value ==~ pattern)) {
        return "${colName} '${value}' is invalid - use ${hint} or '${SENTINEL()}'"
    }
    return null
}

// value must satisfy pattern (or just be non-blank if pattern is null) - used for required columns
def checkRequired(String value, def pattern, String colName, String hint) {
    if (isBlankish(value)) {
        return "required field '${colName}' is missing/blank"
    }
    if (pattern != null && !(value ==~ pattern)) {
        return "${colName} '${value}' is invalid - use ${hint}"
    }
    return null
}

// required file path: non-blank and must exist
def checkFile(String value, String colName) {
    def err = checkRequired(value, null, colName, "a file path")
    if (err) return err
    if (!file(value).exists()) {
        return "${colName} '${value}' does not exist or is not readable"
    }
    return null
}

// Reads a TSV into [header, rows], stripping a UTF-8 BOM and skipping blank lines.
// Each row is [rowNum, fields].
def readTsv(String path) {
    def lines = new File(path).readLines()
    if (lines.isEmpty()) return null
    def header = lines[0].replace('﻿', '').split('\t', -1)*.trim()
    def rows = []
    lines.drop(1).eachWithIndex { line, i ->
        if (line.trim().isEmpty()) return
        rows << [i + 1, line.split('\t', -1)*.trim()]
    }
    return [header, rows]
}

// Amplicon design manifest. Headers are matched loosely (see
// bin/pipeline_common.py), extra columns (GC_percent, final_* etc.) are allowed.
def validateAmpliconManifestContent(String path) {
    // logical name -> normalised header key, as in bin/pipeline_common.py FIELD_KEYS
    def AMPLICON_COLUMNS = [
        project_name: "projectname", plate_name: "platename", well_position: "wellposition",
        gene_name: "genename", gene_id: "geneid", species: "species",
        guide1: "guide1", guide2: "guide2", guide3: "guide3",
        amplicon_size: "ampliconsize", fwd_primer: "fwdprimer", rev_primer: "revprimer"
    ]
    // Optional: absent in manifests written before the column existed, so a missing
    // header is not an error. Mirrors bin/pipeline_common.py OPTIONAL_FIELD_KEYS.
    def OPTIONAL_AMPLICON_COLUMNS = [amplicon_seq: ["ampliconseq", "ampliconsequence"]]
    // guides/primers may legitimately be NA (e.g. non-targeting controls have no amplicon)
    def REQUIRED_VALUES = ["project_name", "plate_name", "well_position", "gene_name"]

    def tsv = readTsv(path)
    if (tsv == null) return ["file is empty"]
    def (header, rows) = tsv

    def normIndex = header.withIndex().collectEntries { name, i -> [(normaliseHeader(name)): i] }
    def missing = AMPLICON_COLUMNS.findAll { logical, key -> !normIndex.containsKey(key) }.keySet()
    if (missing) {
        return ["header is missing required column(s) ${missing} (matched ignoring case/punctuation). Found: ${header}".toString()]
    }
    def colIndex = AMPLICON_COLUMNS.collectEntries { logical, key -> [(logical): normIndex[key]] }
    OPTIONAL_AMPLICON_COLUMNS.each { logical, keys ->
        def key = keys.find { k -> normIndex.containsKey(k) }
        if (key != null) colIndex[logical] = normIndex[key]
    }

    def errors = []
    def warnings = []
    def sanitizedProjects = [:].withDefault { [] as Set }

    rows.each { rowNum, fields ->
        if (fields.size() < header.size()) {
            errors << "row ${rowNum}: expected ${header.size()} tab-separated columns, got ${fields.size()}".toString()
            return
        }
        def project = fields[colIndex.project_name]
        def context = isBlankish(project) ? "row ${rowNum}" : "row ${rowNum} (project_name=${project})"

        def rowErrors = REQUIRED_VALUES.collect { col -> checkRequired(fields[colIndex[col]], null, col, "") }
        rowErrors.findAll { it != null }.each { errors << "${context}: ${it}".toString() }

        // A supplied amplicon_seq is what every edit % for the row gets measured
        // against, so a mis-pasted cell is worth failing on rather than aligning.
        def suppliedSeq = colIndex.amplicon_seq == null ? "" : fields[colIndex.amplicon_seq]
        def hasSuppliedSeq = !isBlankish(suppliedSeq)
        if (hasSuppliedSeq && !(suppliedSeq.toUpperCase() ==~ /^[ACGTN]+$/)) {
            def shown = suppliedSeq.size() > 40 ? "${suppliedSeq.take(40)}... (${suppliedSeq.size()} chars)" : suppliedSeq
            errors << "${context}: amplicon_seq '${shown}' contains characters other than ACGTN".toString()
        }

        // Warn rather than error: a row with no primers and no sequence yields no
        // amplicon, which is how non-targeting controls are legitimately expressed.
        // It just can't be quantified, so say so up front instead of failing the run.
        def hasPrimers = !isBlankish(fields[colIndex.fwd_primer]) && !isBlankish(fields[colIndex.rev_primer])
        if (!hasPrimers && !hasSuppliedSeq) {
            warnings << "${context}: no usable primer pair and no amplicon_seq - no amplicon can be built for this row".toString()
        }

        if (!isBlankish(project)) sanitizedProjects[sanitizeName(project)] << project
    }

    if (rows.isEmpty()) {
        errors << "has a header but no data rows"
    }

    // project_name is used to name the per-project amplicon files, so two
    // projects must not collapse to the same file name
    sanitizedProjects.findAll { name, projects -> projects.size() > 1 }.each { name, projects ->
        errors << "project_name values ${projects.sort()} all map to the same name '${name}' once special characters are replaced - rename one of them".toString()
    }

    return [errors: errors, warnings: warnings]
}

// Fastq manifest: one row per sample, linked to the amplicon manifest by project_name
def validateFastqManifestContent(String path) {
    def REQUIRED_COLUMNS = ["sample", "project_name", "fastq_1", "fastq_2"] as Set
    def OPTIONAL_COLUMNS = ["donor"] as Set

    def tsv = readTsv(path)
    if (tsv == null) return ["file is empty"]
    def (header, rows) = tsv

    def missing = REQUIRED_COLUMNS - (header as Set)
    def extra = (header as Set) - REQUIRED_COLUMNS - OPTIONAL_COLUMNS
    if (missing || extra) {
        def msg = "header does not match expected columns ${REQUIRED_COLUMNS} (optional: ${OPTIONAL_COLUMNS})."
        if (missing) msg += " Missing: ${missing}."
        if (extra) msg += " Unexpected: ${extra}."
        return [msg.toString()]
    }
    if (header.size() != (header as Set).size()) {
        return ["header has duplicated column names: ${header}".toString()]
    }

    def colIndex = header.withIndex().collectEntries { name, i -> [(name): i] }
    def hasDonor = colIndex.containsKey("donor")

    def errors = []
    def slugs = [:].withDefault { [] }

    rows.each { rowNum, fields ->
        if (fields.size() != header.size()) {
            errors << "row ${rowNum}: expected ${header.size()} tab-separated columns, got ${fields.size()}".toString()
            return
        }
        def sample = fields[colIndex.sample]
        def context = isBlankish(sample) ? "row ${rowNum}" : "row ${rowNum} (sample=${sample})"

        def rowErrors = []
        rowErrors << checkRequired(sample, null, "sample", "a sample name")
        rowErrors << checkRequired(fields[colIndex.project_name], null, "project_name", "a project_name from the amplicon manifest")
        rowErrors << checkFile(fields[colIndex.fastq_1], "fastq_1")
        rowErrors << checkFile(fields[colIndex.fastq_2], "fastq_2")
        if (hasDonor) {
            rowErrors << checkOptional(fields[colIndex.donor], null, "donor", "a donor name")
        }
        rowErrors.findAll { it != null }.each { errors << "${context}: ${it}".toString() }

        if (!isBlankish(sample)) slugs[slugifySample(sample)] << sample
    }

    if (rows.isEmpty()) {
        errors << "has a header but no data rows"
    }

    // Sample names become CRISPResso output folder names, so they must stay
    // unique after CRISPResso's own slugification
    slugs.findAll { slug, samples -> samples.size() > 1 }.each { slug, samples ->
        def msg = samples.toSet().size() == 1 ?
            "sample '${samples[0]}' appears on ${samples.size()} rows - sample names must be unique" :
            "samples ${samples} all become '${slug}' in CRISPResso output names - sample names must be unique"
        errors << msg.toString()
    }

    return errors
}

def formatValidationErrors(String label, List errors) {
    if (!errors) return null
    def lines = errors.collect { "  ${it}" }.join('\n')
    return "${label} has ${errors.size()} problem(s):\n${lines}"
}
