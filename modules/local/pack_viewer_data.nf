// Packs one sample's CRISPRessoPooled per-amplicon tables into a single small
// gzipped JSON for the report's amplicon viewer (bin/pack_viewer_data.py).
//
// Runs once per sample rather than inside RENDER_REPORT so the thousands of
// per-amplicon files are parsed in parallel and never staged into the report task.
// The output is named after the sample: Nextflow stages inputs by basename and does
// not de-duplicate, so a fixed name would collide when these are collected.

process PACK_VIEWER_DATA {
    tag "${sample}"
    label 'small'

    input:
        tuple val(sample), path(pooled_dir)
        val min_reads
        val quantification_window_center

    output:
        path "${sample}.viewer.json.gz", emit: data

    script:
    """
    pack_viewer_data.py \\
        --pooled ${pooled_dir} \\
        --sample ${sample} \\
        --min-reads ${min_reads} \\
        --quantification-window-center ${quantification_window_center} \\
        --out ${sample}.viewer.json.gz
    """
}
