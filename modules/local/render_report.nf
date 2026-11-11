// Summary report aggregated over all samples (bin/render_report.py +
// assets/report.html.j2), laid out like the tglow-pipeline QC report. Runs once,
// after every CRISPRessoPooled task has finished; samples without CRISPResso
// output (failed tasks under an 'ignore' error strategy) are reported as
// crispresso_failed rather than failing the report.

process RENDER_REPORT {
    tag 'report'
    label 'small'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
        path crispresso_results, stageAs: "crispresso/*"
        path viewer_data, stageAs: "viewer/*"
        path sample_check
        path amplicon_replicates
        path amplicon_report
        path guide_qc
        path amplicon_manifest
        path template
        path viewer_assets
        val pipeline_version

    output:
        path 'edit_quant_report.html', emit: html
        path 'pair_summary.tsv', emit: pairs

    script:
    // viewer/ is empty when the viewer is disabled or no sample produced a payload;
    // render_report.py then simply omits the tab.
    // viewer_assets stages amplicon_viewer.css/.body.html/.js into the task directory:
    // the template is staged on its own, so its directory is not the assets directory.
    """
    export MPLCONFIGDIR=\$PWD/.mplconfig
    mkdir -p viewer
    render_report.py \\
        --crispresso_dir crispresso \\
        --viewer_dir viewer \\
        --viewer_assets . \\
        --sample_check ${sample_check} \\
        --amplicon_replicates ${amplicon_replicates} \\
        --amplicon_report ${amplicon_report} \\
        --guide_qc ${guide_qc} \\
        --amplicon_manifest ${amplicon_manifest} \\
        --template ${template} \\
        --edit_threshold ${params.report_edit_threshold} \\
        --wt_max_edit ${params.report_wt_max_edit} \\
        --min_reads ${params.report_min_reads} \\
        --wt_pattern '${params.report_wt_pattern}' \\
        --plate_format ${params.report_plate_format} \\
        --pam '${params.guide_pam}' \\
        --pam_side ${params.guide_pam_side} \\
        --quantification_window_center ${params.quantification_window_center} \\
        --pipeline_version '${pipeline_version}' \\
        --output edit_quant_report.html \\
        --pairs_output pair_summary.tsv
    """
}
