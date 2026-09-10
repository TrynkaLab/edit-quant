# edit_quant

Nextflow (DSL2) pipeline: aligns a designed amplicon panel's guides + primers
to a reference genome with BWA, extracts each amplicon's reference sequence,
QC's primer orientation/gene overlap, then runs
[CRISPRessoPooled](https://github.com/pinellolab/CRISPResso2) on each FASTQ
sample against the resulting (deduplicated) amplicons panel.

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

### `--fastq_manifest` (required)

Tab-separated, header row required, 3 columns:

```
sample	fastq_1	fastq_2
sample_a	/path/to/sample_a_R1.fastq.gz	/path/to/sample_a_R2.fastq.gz
```

Every sample is screened with CRISPRessoPooled against the *same* shared,
deduplicated amplicons panel — there's no per-sample link to a specific
`amplicon_manifest` row required or used.

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
| `crispresso_n_processes` | no | default `1` |
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
  orientation QC
- `reports/off_targets.tsv` — additional low-edit-distance hits
- `reports/reassigned_primaries.tsv` — log of every tie-break correction
- `reports/crispresso_amplicons.txt` — CRISPRessoPooled amplicons file (one
  row per sample whose amplicon was extracted)
- `amplicons/<sample_id>.fa`, `amplicons/all_amplicons.fasta` — extracted
  reference amplicon sequences
- `crispresso/amplicons_unique.txt` — the above, deduplicated by
  `amplicon_seq` (CRISPRessoPooled requires unique sequences)
- `crispresso/CRISPRessoPooled_on_<sample>/` — one per fastq-manifest sample
