#!/usr/bin/env nextflow

// CRISPResso slugifies --name internally (collapses runs of '_' etc. down to a
// single '_'), so a sample name containing e.g. a double underscore would make
// it write "CRISPRessoPooled_on_<slugified>" while Nextflow still expects the
// raw "CRISPRessoPooled_on_<sample>" -- pre-slugify so both sides agree.
def slugifySample(String sample) {
    return sample.replaceAll(/[\s'*"\/\\\[\]:;|,<>?]/, '_').replaceAll(/_{2,}/, '_')
}

// Mirrors bin/pipeline_common.py's sanitize(), which names the per-project
// amplicon files
def sanitizeName(String name) {
    return name.trim().replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

// Column names in the amplicon manifest are matched after stripping
// punctuation/case, mirroring bin/pipeline_common.py
def normaliseHeader(String header) {
    return header.toLowerCase().replaceAll(/[^a-z0-9]/, '')
}
