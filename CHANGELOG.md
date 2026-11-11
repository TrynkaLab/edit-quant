# unreleased
- Fixed: boolean params given on the command line (`--reassign_primary false`, `--skip_checks false`) arrived as the String `"false"`, which is truthy in Groovy, so they behaved exactly like `true`. Coerced via a new `asBool` helper in `lib/utils.nf`
- New `crispresso_skip_failed` param (default `true`), passing CRISPResso's own `--skip_failed` to CRISPRessoPooled: an amplicon that errors is dropped and the rest of that sample's panel is still quantified, instead of the whole task failing. Skipped amplicons are reported as `not_tested` (`missing_from_output`) rather than lost
- Amplicon sequences can now be supplied directly via an optional `amplicon_seq` manifest column (also accepted as `amplicon_sequence`), for amplicons with no genomic locus (transgenes, plasmids) or whose true sequence differs from the reference. Manifests without the column are unaffected
- A supplied sequence always wins, but the reference is still extracted where the primers map and the two are compared: new `amplicon_source`, `reference_check` (`match`/`match_revcomp`/`contains_reference`/`mismatch`/`no_reference`), `primers_in_supplied` (`bounds`/`interior`/`no`/`NA`) and `reference_length` columns in `reports/amplicon_report.tsv`, surfaced in the report's Design QC tab; new `SUPPLIED_DIFFERS_FROM_REFERENCE` and `SUPPLIED_PRIMERS_NOT_FOUND` statuses
- Primers may be `NA` on rows supplying a sequence; a row with neither a usable primer pair nor an `amplicon_seq` is warned about rather than failing the run
- Fixed: reference sequences are now uppercased on extraction. `samtools faidx` reproduces the reference's own case and the standard human references are soft-masked, so any guide landing in a repeat-masked region was reported `NOT_IN_AMPLICON` and dropped from the CRISPResso amplicons file, widening its quantification window to the whole amplicon
- Fixed: `NA` in a guide/primer cell no longer becomes a 2 bp `NA` query sequence in `queries.fasta` (both N and A are valid bases, so it passed the emptiness check and aligned); `none`/`null`/`n/a` are now honoured as absent everywhere, matching the manifest validator
- New `EMPTY_AMPLICON` status guards against a zero-length extraction reaching CRISPResso
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