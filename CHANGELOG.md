# v0.0.2
- fastq_manifest now requires a `project_name` column and takes an optional `donor` column
- Each sample is screened only against amplicons with a matching amplicon-manifest `project_name`
- New CHECK_SAMPLES step: samples without amplicons for their project are warned about and skipped (`reports/sample_check.tsv`)
- Input manifests are validated up front (`lib/checks.nf`, `lib/validateInputs.nf`); `--skip_checks` to bypass
- Replaces DEDUP_AMPLICONS / `crispresso/amplicons_unique.txt` with `crispresso/amplicons_by_project/`
- New RENDER_REPORT step: `report/edit_quant_report.html` (Summary, QC, Plates, Amplicon x donor tabs, tglow QC report layout) and `report/pair_summary.tsv`; `report_*` params
- WT samples (`report_wt_pattern`) get their own QC status: tested WT pairs with edit % <= `report_wt_max_edit` (default 10) are `success_wt` instead of `fail`; WT pairs above it fail
- New `crispresso_min_reads_to_use_region` param (default 100), passed to CRISPRessoPooled; its built-in default of 1000 skipped amplicons with 1-1000 reads (37% of own-project sample x amplicon pairs not analysed in the v0.0.1 run)
- Technical replicates (wells in a project sharing an amplicon) are merged into one named amplicon with the union of their guides, instead of keeping only the first well; mapping in `reports/amplicon_replicates.tsv`

# v0.0.1
Initial implementation