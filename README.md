# edit_quant

Nextflow (DSL2) pipeline: aligns a designed amplicon panel's guides + primers
to a reference genome with BWA, extracts each amplicon's reference sequence,
QC's primer orientation/gene overlap, then runs
[CRISPRessoPooled](https://github.com/pinellolab/CRISPResso2) on each FASTQ
sample against the (deduplicated) amplicons of that sample's project.

Translated from the bash/Python pipeline in
`OTAR3086/analysis/mapping_amplicons/scripts` — see that directory's
`README.md` for the full rationale behind the two-BWA-pass design, primer
reassignment, gene-overlap window, high-copy-number flagging, etc. All of
that logic is reused unchanged (`bin/*.py`), just re-wired as Nextflow
processes.

## Inputs

### `--amplicon_manifest` (required)

Tab-separated, header row required. Column names are matched after
stripping punctuation/case (`guide_.1`/`guide_1`/`Guide1` all work):
`project_name, plate_name, well_position, gene_name, gene_id, species,
guide1, guide2, guide3, amplicon_size, fwd_primer, rev_primer`. Guides are
given as RNA (`U`) and converted to DNA; primers are the plain
`fwd_primer`/`rev_primer` columns (not adaptor-tagged `final_*` columns).

One optional column is recognised: **`amplicon_seq`** (also accepted as
`amplicon_sequence`). Supply it to use a known amplicon sequence instead of
extracting one from the reference — for amplicons with no genomic locus
(transgenes, plasmids, integrated constructs), or whose true sequence differs
from the reference. A supplied sequence always wins, but where the primers
also map the reference sequence is still extracted and compared, and any
disagreement is flagged in `reports/amplicon_report.tsv` and the report's
Design QC tab. Primers may be `NA` on rows that supply a sequence; a row with
neither a usable primer pair nor an `amplicon_seq` is warned about and
produces no amplicon.

### `--fastq_manifest` (required)

Tab-separated, header row required. Columns `sample`, `project_name`,
`fastq_1`, `fastq_2` are required; `donor` is optional (use `none` where
not applicable). Column names must match exactly.

```
sample	project_name	donor	fastq_1	fastq_2
sample_a	proj-1	A	/path/to/sample_a_R1.fastq.gz	/path/to/sample_a_R2.fastq.gz
```

`project_name` links each sample to the `amplicon_manifest` rows with the
same `project_name`: each sample is screened with CRISPRessoPooled against
only its own project's amplicons. Samples whose project has no amplicons
(not in the amplicon manifest, or none could be extracted) are reported
with a warning and skipped; see `reports/sample_check.tsv`.

### Input checks

Both manifests are checked before anything runs (`lib/checks.nf`,
`lib/validateInputs.nf`): headers, required values, that fastq files exist,
and that sample names are unique (also after CRISPResso's own renaming).
All problems in a file are reported together. Blank-like values (`NA`,
`null`, empty, ...) are rejected; use `none` for an optional value that
doesn't apply. `--skip_checks` turns these checks off.

### Reference / config params

| param | required? | notes |
|---|---|---|
| `ref_fasta` | yes | reference genome FASTA |
| `ref_fai` | no | default `${ref_fasta}.fai`; built with `samtools faidx` if missing |
| `bwa_index_prefix` | no | default `ref_fasta`; built with `bwa index` if missing |
| `gene_bed` | one of `gene_bed`/`genes_gtf` | prebuilt gene BED (chrom,start0,end,gene_id,gene_name,strand) |
| `exon_bed` | no | prebuilt exon BED; adds `overlapping_exons` to `validation_report.tsv` |
| `genes_gtf` | one of `gene_bed`/`genes_gtf` | builds `gene_bed`/`exon_bed` from an Ensembl GTF if `gene_bed` isn't given |
| `bwa_n` | no | default `2` — off-target pass max edit distance |
| `max_xa` | no | default `1000` — max BWA alt hits reported in the off-target pass |
| `max_mismatch` | no | default `2` — edit-distance threshold for `off_targets.tsv` |
| `gene_window` | no | default `500` — bp padding around genes for overlap checks |
| `high_copy_threshold` | no | default `50` — flags `HIGH_COPY_NUMBER` above this many tied-best hits |
| `reassign_primary` | no | default `true` — promote gene-matching tied alternates to primary |
| `guide_pam` | no | default `NGG` — PAM adjacent to each guide (IUPAC) for the `guide_qc.tsv` orientation check; `TTTV` for Cas12a/Cpf1 |
| `guide_pam_side` | no | default `3prime` — which end of the protospacer that PAM sits at; `5prime` for Cas12a/Cpf1 |
| `quantification_window_center` | no | default `-3` — CRISPResso's `-wc`, mirrored to report each guide's cut index; `1` for Cas12a/Cpf1. Reporting only, not passed to CRISPResso |
| `crispresso_n_processes` | no | default `4` — cores per CRISPRessoPooled task. Also sets the process's `cpus`, and is passed to CRISPResso as `--n_processes ${task.cpus}` so the scheduler request and the flag can't drift |
| `crispresso_memory` | no | default `16.GB` — memory per CRISPRessoPooled task, multiplied by the attempt number on retry. CRISPResso's memory use rises with `crispresso_n_processes`, so raise the two together |
| `crispresso_suppress_plots` | no | default `true` — pass `--suppress_plots` to CRISPRessoPooled. A large speed-up on big panels; this pipeline's report is built from `SAMPLES_QUANTIFICATION_SUMMARY.txt`, not CRISPResso's plots, so nothing in it is lost. CRISPResso's own HTML report is gated by a separate flag and is still written |
| `crispresso_skip_failed` | no | default `true` — pass `--skip_failed` to CRISPRessoPooled, so an amplicon that errors is dropped and the rest of the sample's panel still runs. Skipped amplicons show up as `not_tested` (`missing_from_output`) in the report |
| `crispresso_min_reads_to_use_region` | no | default `100` — CRISPRessoPooled skips amplicons with fewer assigned reads (CRISPResso's own default, 1000, left ~1/3 of amplicons unanalysed); keep in line with `report_min_reads` |
| `amplicon_viewer` | no | default `true` — embed the Amplicons tab in the HTML report: per-base read composition against the reference, with guides, cut sites and the quantification window drawn underneath. Data is packed once per sample by `PACK_VIEWER_DATA`, adding roughly 1.8 KB per amplicon to the report |
| `amplicon_viewer_min_reads` | no | default `0` — skip amplicons below this aligned-read count in the viewer tab only; the QC tab and `pair_summary.tsv` always report every amplicon |
| `report_edit_threshold` | no | default `50` — edit % at or above which a tested pair is a successful edit |
| `report_wt_max_edit` | no | default `10` — edit % at or below which a tested WT pair is `success_wt` (no editing); above it the pair fails |
| `report_min_reads` | no | default `100` — min aligned reads for a pair to count as tested |
| `report_wt_pattern` | no | default `(^\|_)WT$` — regex on `donor` (or `sample` if donor is `none`); matching samples are reported separately |
| `report_plate_format` | no | default `auto` — plate heatmap layout (`auto`, `24`, `96`, `384`, `1536`) |
| `skip_checks` | no | default `false` — skip the input checks above |
| `outdir` | no | default `./results` |

If `ref_fai`/BWA index/`gene_bed`+`exon_bed` already exist on disk, the
pipeline uses them as-is (no rebuild). Anything auto-built is written under
`${outdir}/reference/`, never into the original reference directory (which
may be shared/read-only).

## Environment

All processes run inside one conda environment auto-built by Nextflow from
`environment.yml` (`conda.enabled = true`, `process.conda` points at the
file) — `bwa`, `samtools`, `crispresso2`, `python`. No manual env setup
required.

## Running

```bash
nextflow run main.nf \
  --amplicon_manifest /path/to/amplicon_manifest.tsv \
  --fastq_manifest    /path/to/fastq_manifest.tsv \
  --ref_fasta         /path/to/genome.fasta \
  --genes_gtf         /path/to/genes.gtf.gz \
  --outdir            ./results
```

On the Sanger LSF farm, add `-profile farm` (see `conf/farm.config` — adjust
the placeholder queue name to your actual queues before relying on it).

### `edit_quant` runner script

A `tglow-pipeline`-style wrapper (`./edit_quant`) is included for driving runs
via LSF (or locally) with sensible logging/`-resume` defaults baked in:

```bash
# local (no LSF submission)
./edit_quant -l -- \
  --amplicon_manifest /path/to/amplicon_manifest.tsv \
  --fastq_manifest    /path/to/fastq_manifest.tsv \
  --ref_fasta         /path/to/genome.fasta \
  --genes_gtf         /path/to/genes.gtf.gz \
  --outdir            ./results

# submit the nextflow driver to LSF (per-process resources still come from
# conf/farm.config's executor settings)
./edit_quant -p farm -g <your_lsf_group> -- \
  --amplicon_manifest /path/to/amplicon_manifest.tsv \
  --fastq_manifest    /path/to/fastq_manifest.tsv \
  --ref_fasta         /path/to/genome.fasta \
  --gene_bed          /path/to/genome.genes.bed \
  --outdir            ./results
```

Everything after `--` is passed straight through to `nextflow run` (so any
pipeline `--param` works there), an optional `-c run.config` can supply
further overrides, and `-w` sets the Nextflow work directory (default
`../workdir`, one level above this pipeline folder). Each run gets its own
`logs/<timestamp>/` with the Nextflow log, HTML report, and trace. See
`./edit_quant -h` for all flags. The proxy/module-load lines near the top of
the script are Sanger-farm-specific — adjust for your own installation.

## Outputs (under `--outdir`)

- `fasta/queries.fasta` — guide/primer query sequences
- `align/queries.exact.sorted.bam`(+`.bai`) — exact-match alignment (source
  of truth for position/gene-overlap/amplicon extraction)
- `align/queries.sorted.bam`(+`.bai`) — mismatch-tolerant alignment (off-target search)
- `align/queries.exact.reassigned.sorted.bam`(+`.bai`) — corrected BAM, only
  when `reassign_primary` corrected at least one position
- `reports/validation_report.tsv` — per guide/primer: position, gene/exon
  overlap, status
- `reports/amplicon_report.tsv` — per sample: amplicon coordinates, length,
  orientation QC. When a primer has several equally-good hits that all pass
  the gene/orientation checks (e.g. it sits in a tandem repeat), the hit
  giving the length closest to `amplicon_size` is used; any remaining
  mismatch with `amplicon_size` is flagged `LENGTH_MISMATCH` instead of `OK`.
  `amplicon_source` says whether the sequence came from the `reference` or the
  manifest's `amplicon_seq` column. For a supplied sequence, `reference_check`
  compares it with the reference at the primers' coordinates — `match`,
  `match_revcomp` (a minus-strand design, which reads as the reverse complement
  of the plus-strand extraction), `contains_reference` (a construct carrying the
  amplicon), `mismatch`, or `no_reference` — alongside `reference_length` so a
  mismatch can be sized; `primers_in_supplied` says whether the sequence *is*
  the PCR product (`bounds`) or merely contains it (`interior`), and the
  declared-size check only applies to the former. Supplied-sequence statuses:
  `SUPPLIED_DIFFERS_FROM_REFERENCE`, `SUPPLIED_PRIMERS_NOT_FOUND`
- `reports/guide_qc.tsv` — per designed guide: where (and on which strand) it
  matches its extracted amplicon, the PAM-proximal flank, the cut index
  CRISPResso will derive, whether it was passed to CRISPResso, and a status.
  CRISPResso accepts a guide on either strand of the amplicon, but reads it as
  the protospacer 5'->3' with the PAM immediately 3' of it, so the sequence
  must be the designed spacer and not its reverse complement — it matches
  either way round, without warning, and silently places the quantification
  window at the wrong end of the guide (~14 bp out for a 20mer at `-wc -3`).
  Statuses: `OK`, `PAM_ON_WRONG_END` (the manifest holds the reverse complement
  of the designed protospacer), `NO_PAM`, `PAM_OUTSIDE_AMPLICON` (guide runs up
  against the amplicon edge), `MULTIPLE_MATCHES`, `NOT_IN_AMPLICON` (dropped
  from the CRISPResso amplicons file, since CRISPResso would warn and ignore
  it anyway, widening the window to the whole amplicon) and `NO_AMPLICON`
- `reports/off_targets.tsv` — additional low-edit-distance hits
- `reports/reassigned_primaries.tsv` — log of every tie-break correction
- `reports/crispresso_amplicons.txt` — CRISPRessoPooled amplicons file (one
  row per sample whose amplicon was extracted). `guide_seq` holds the designed
  guides from the manifest (U→T), passed through as written so the orientation
  the design carries reaches CRISPResso intact; guides flagged
  `NOT_IN_AMPLICON` in `guide_qc.tsv` are left out
- `amplicons/<sample_id>.fa`, `amplicons/all_amplicons.fasta` — the amplicon
  sequence used for each sample, extracted from the reference or taken from the
  manifest's `amplicon_seq`. The header's coordinate field reads `no_locus` for
  a supplied sequence with no reference position
- `reports/sample_check.tsv` — per fastq-manifest sample: project, donor,
  number of amplicons and status (`OK`, `NO_PROJECT_IN_AMPLICON_MANIFEST`,
  `NO_AMPLICONS_EXTRACTED`); only `OK` samples are run
- `crispresso/amplicons_by_project/<project_name>.amplicons.txt` — the
  above split by `project_name`. CRISPRessoPooled requires unique sequences,
  and technical replicates (wells in a project sharing an amplicon) can't be
  separated in a sample's pooled reads, so each group is merged into one
  amplicon named after all its wells (e.g. `IL2RA_G4-D8-B9_<plate>_<project>`)
  with the union of their guides
- `reports/amplicon_replicates.tsv` — maps every well's amplicon name to the
  merged CRISPResso amplicon name, with the number of replicates and whether
  their guides differed
- `crispresso/CRISPRessoPooled_on_<sample>/` — one per `OK` sample in
  `sample_check.tsv`
- `report/edit_quant_report.html` — self-contained summary report over all
  samples, laid out like the tglow-pipeline QC report:
  - *Summary*: sample/amplicon/pair counts, successful edits and mean
    edit %, with WT samples (`report_wt_pattern`) in their own column so they
    don't contribute; histogram of aligned reads per pair
  - *QC*: success/success_wt/fail/not-tested pairs per sample, and a filterable table
    of every pair
  - *Plates*: edit % heatmap per design plate and sample
  - *Amplicon x donor*: edit % heatmap per project
- `report/pair_summary.tsv` — the report's pair table: one row per sample x
  amplicon with reads, edit % and status. For edited samples, `success`/`fail`
  are tested pairs (&ge; `report_min_reads` aligned reads) above/below
  `report_edit_threshold`. WT samples are expected to be unedited, so their
  tested pairs are `success_wt` at or below `report_wt_max_edit` and `fail`
  above it;
  `not_tested` pairs give a reason (`no_reads`, `not_aligned`, `low_depth`,
  `crispresso_failed`)
