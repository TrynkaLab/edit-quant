# unreleased
- CRISPResso guide sequences now come from the design manifest (U→T) instead of being re-extracted from the reference at each guide's aligned coordinates, which returned the + strand and so silently reverse-complemented every minus-strand guide — CRISPResso matches either orientation without warning but reads the 3' end as PAM-proximal, putting the quantification window ~14 bp out for a 20mer
- New guide orientation/PAM QC: `reports/guide_qc.tsv`, one row per designed guide, with its match position and strand in the extracted amplicon, PAM-proximal flank, the cut index CRISPResso will derive, and a status (`OK`, `PAM_ON_WRONG_END`, `NO_PAM`, `PAM_OUTSIDE_AMPLICON`, `MULTIPLE_MATCHES`, `NOT_IN_AMPLICON`, `NO_AMPLICON`)
- Guides CRISPResso can't find in their amplicon are no longer passed to it (it would warn, ignore them and widen the window to the whole amplicon); they're reported as `NOT_IN_AMPLICON` instead
- New report tab "Design QC" between Summary and QC, with an amplicon-extraction section and a guide section, each with status counts and a filterable table
- New `guide_pam` (default `NGG`), `guide_pam_side` (default `3prime`) and `quantification_window_center` (default `-3`) params driving the orientation check and reported cut index; set `TTTV`/`5prime`/`1` for Cas12a/Cpf1

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