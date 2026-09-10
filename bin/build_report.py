#!/usr/bin/env python3
"""Turn a BWA alignment of guides/primers into:

  - reports/validation_report.tsv    one row per guide/primer: position, gene overlap
  - reports/amplicon_report.tsv      one row per sample: amplicon coordinates & QC
  - reports/off_targets.tsv          additional high-quality (low edit-distance) hits
  - reports/reassigned_primaries.tsv every position corrected to a gene-matching tie
  - reports/crispresso_amplicons.txt CRISPRessoPooled amplicon description file
  - amplicons/<sample_id>.fa         extracted reference sequence of each amplicon
  - amplicons/all_amplicons.fasta    all of the above concatenated
  - <exact-sam>.reassigned.sam       (only with --reassign-primary) corrected SAM

Two alignments are expected:
  --exact-sam   bwa aln -n 0 (exact matches only) -- the source of truth for
                position/gene-overlap/amplicon-extraction. Restricting to 0
                mismatches keeps tie counts naturally small (a handful, not
                hundreds), so BWA's XA tag never gets silently suppressed the
                way it can under a mismatch-tolerant search, and no
                per-read re-alignment is needed to resolve ties.
  --sam         the mismatch-tolerant search (BWA_N/MAX_XA) -- used only for
                off_targets.tsv. If a read has no exact hit at all (e.g. a
                real SNP at the site), its position falls back to this run.

Usage:
  build_report.py --tsv input.tsv --exact-sam queries.exact.sam --sam queries.sam \
      --ref genome.fa --gene-bed genes.bed [--exon-bed exons.bed] --outdir OUTDIR \
      [--max-mismatch 2] [--reassign-primary]
"""
import argparse
import bisect
import os
import subprocess
import sys
from dataclasses import dataclass, field
from typing import List, Optional

sys.path.insert(0, os.path.dirname(__file__))
from pipeline_common import read_rows, SEQ_TYPES  # noqa: E402

REF_CONSUMING_OPS = set("MDN=X")
CIGAR_TOKEN_RE = __import__("re").compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class Hit:
    chrom: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str
    mapq: Optional[int]
    nm: Optional[int]
    source: str  # "primary" or "XA"
    cigar: str = ""


@dataclass
class QueryResult:
    sample_id: str
    seq_type: str
    gene_name_hdr: str
    gene_id_hdr: str
    primary: Optional[Hit] = None
    alt: List[Hit] = field(default_factory=list)
    n_best_hits: Optional[int] = None  # from X0 tag: count of hits tied at the best edit distance
    n_suboptimal_hits: Optional[int] = None  # from X1 tag: count of additional worse-but-within-threshold hits


def cigar_ref_len(cigar):
    if cigar == "*":
        return 0
    total = 0
    for length, op in CIGAR_TOKEN_RE.findall(cigar):
        if op in REF_CONSUMING_OPS:
            total += int(length)
    return total


def parse_qname(qname):
    """sample_id and seq_type are the only load-bearing fields (used as the
    lookup key everywhere below); gene_name/gene_id in the header are purely
    informational, so a header with too few/many fields (e.g. from a
    differently-versioned tsv_to_fasta.py) still parses instead of silently
    dropping the alignment. Returns None if even sample_id|seq_type can't be
    recovered."""
    parts = qname.split("|")
    if len(parts) < 2:
        return None
    sample_id, seq_type = parts[0], parts[1]
    gene_name_hdr = parts[2] if len(parts) > 2 else ""
    gene_id_hdr = parts[3] if len(parts) > 3 else ""
    return sample_id, seq_type, gene_name_hdr, gene_id_hdr


def parse_sam(sam_path):
    """Return dict[(sample_id, seq_type)] -> QueryResult."""
    results = {}
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            qname, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
            tags = {}
            for tag in f[11:]:
                tparts = tag.split(":", 2)
                if len(tparts) == 3:
                    tags[tparts[0]] = tparts[2]

            parsed = parse_qname(qname)
            if parsed is None:
                sys.stderr.write(
                    "WARNING: cannot parse query name %r (expected at least "
                    "sample_id|seq_type), skipping\n" % qname
                )
                continue
            sample_id, seq_type, gene_name_hdr, gene_id_hdr = parsed

            key = (sample_id, seq_type)
            qr = results.setdefault(
                key, QueryResult(sample_id, seq_type, gene_name_hdr, gene_id_hdr)
            )

            unmapped = bool(flag & 4)
            if not unmapped:
                strand = "-" if (flag & 16) else "+"
                nm = int(tags["NM"]) if "NM" in tags else None
                qr.primary = Hit(rname, pos, pos + cigar_ref_len(cigar) - 1, strand, mapq, nm, "primary", cigar)
            if "X0" in tags:
                try:
                    qr.n_best_hits = int(tags["X0"])
                except ValueError:
                    pass
            if "X1" in tags:
                try:
                    qr.n_suboptimal_hits = int(tags["X1"])
                except ValueError:
                    pass

            xa = tags.get("XA", "")
            for entry in xa.split(";"):
                if not entry:
                    continue
                try:
                    xchrom, xpos_s, xcigar, xnm_s = entry.split(",")
                except ValueError:
                    continue
                xstrand = "-" if xpos_s.startswith("-") else "+"
                xpos = int(xpos_s.lstrip("+-"))
                xnm = int(xnm_s) if xnm_s != "" else None
                qr.alt.append(
                    Hit(xchrom, xpos, xpos + cigar_ref_len(xcigar) - 1, xstrand, None, xnm, "XA", xcigar)
                )
    return results


def format_xa(alt_hits):
    """Rebuild a BWA-style XA tag value ("chrom,+pos,cigar,nm;...") from a list of Hits."""
    if not alt_hits:
        return None
    return ";".join(
        "%s,%s%d,%s,%d" % (h.chrom, h.strand, h.start, h.cigar, h.nm) for h in alt_hits
    ) + ";"


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def write_reassigned_sam(orig_sam_path, out_sam_path, reassignments):
    """Copy orig_sam_path (the exact-match SAM) to out_sam_path, rewriting
    FLAG/RNAME/POS/CIGAR/NM/XA/SEQ/QUAL for any (sample_id, seq_type) present
    in `reassignments` so a BAM built from the output reflects the corrected
    primary alignment. MAPQ is left untouched -- it already reflects the
    read's overall mapping ambiguity, which doesn't change just because a
    different equally-good location was chosen."""
    with open(orig_sam_path) as fh, open(out_sam_path, "w") as out:
        for line in fh:
            if line.startswith("@"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            parsed = parse_qname(f[0])
            r = reassignments.get((parsed[0], parsed[1])) if parsed else None
            if r is None:
                out.write(line)
                continue

            flag = int(f[1])
            old_strand = "-" if (flag & 0x10) else "+"
            f[1] = str((flag & ~0x10) | (0x10 if r["strand"] == "-" else 0))
            f[2] = r["chrom"]
            f[3] = str(r["pos"])
            f[5] = r["cigar"]
            if r["strand"] != old_strand:
                # SEQ/QUAL are stored relative to FLAG's strand bit (SAM convention) --
                # a strand-flipping reassignment must revcomp them to match, or the
                # record claims one strand while its bases represent the other.
                f[9] = revcomp(f[9])
                if f[10] != "*":
                    f[10] = f[10][::-1]

            rest, found_nm, found_xa = [], False, False
            for tag in f[11:]:
                name = tag.split(":", 1)[0]
                if name == "NM":
                    rest.append("NM:i:%d" % r["nm"])
                    found_nm = True
                elif name == "XA":
                    if r["xa"]:
                        rest.append("XA:Z:%s" % r["xa"])
                    found_xa = True
                else:
                    rest.append(tag)
            if not found_nm:
                rest.append("NM:i:%d" % r["nm"])
            if not found_xa and r["xa"]:
                rest.append("XA:Z:%s" % r["xa"])

            out.write("\t".join(f[:11] + rest) + "\n")


def normalize_chrom(chrom):
    """Make 'chr1'/'1' and 'chrM'/'MT'/'M' compare equal, so a gene BED and a
    BAM/genome that disagree on chromosome-naming convention still match."""
    c = chrom.strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    if c.upper() in ("M", "MT"):
        c = "MT"
    return c.upper()


class IntervalIndex:
    """Generic in-memory interval lookup: chrom -> intervals sorted by start.

    Each BED row is chrom, start0, end, then `n_extra` arbitrary string
    fields, which overlapping() returns verbatim as a tuple.
    """

    def __init__(self, bed_path, n_extra, label, padding=0):
        self.by_chrom = {}
        n = 0
        with open(bed_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                cols = line.split("\t")
                chrom, start, end = cols[0], int(cols[1]), int(cols[2])
                if padding:
                    start = max(0, start - padding)
                    end = end + padding
                extra = tuple(cols[3:3 + n_extra])
                key = normalize_chrom(chrom)
                self.by_chrom.setdefault(key, []).append((start, end) + extra)
                n += 1
        for chrom in self.by_chrom:
            self.by_chrom[chrom].sort(key=lambda t: t[0])
        self._starts = {c: [iv[0] for iv in ivals] for c, ivals in self.by_chrom.items()}
        sys.stderr.write("%s: loaded %d intervals on %d chroms\n" % (label, n, len(self.by_chrom)))

    def overlapping(self, chrom, start_1based, end_1based):
        """start/end are 1-based inclusive (SAM-style); BED is 0-based half-open."""
        key = normalize_chrom(chrom)
        ivals = self.by_chrom.get(key)
        if not ivals:
            return []
        starts = self._starts[key]
        q_start0 = start_1based - 1
        # intervals starting after end_1based can't overlap; scan backwards from there
        hi = bisect.bisect_right(starts, end_1based)
        out = []
        for i in range(hi - 1, -1, -1):
            iv_start, iv_end, *extra = ivals[i]
            if iv_end <= q_start0:
                continue
            if iv_start < end_1based:
                out.append(tuple(extra))
        return out


def gene_list_str(overlaps):
    if not overlaps:
        return ""
    return ";".join("%s:%s" % (gid, gname) for gid, gname, _ in overlaps)


def gene_names_str(overlaps):
    if not overlaps:
        return ""
    return ";".join(sorted({gname for _, gname, _ in overlaps}))


def exon_list_str(exon_overlaps, has_gene_overlap, exon_index_available):
    """exon_overlaps: list of (gene_id, gene_name, transcript_id, exon_number, strand).

    "intronic" only means something if exon annotation was actually loaded --
    otherwise every gene-overlapping hit would look intronic just because
    there was nothing to check it against.
    """
    if not exon_index_available:
        return ""
    if exon_overlaps:
        pairs = sorted({(gname, enum) for _, gname, _, enum, _ in exon_overlaps})
        return ";".join("%s:exon%s" % (gname, enum) for gname, enum in pairs)
    return "intronic" if has_gene_overlap else ""


def declared_gene_in_overlaps(overlaps, gene_id, gene_name):
    return any(gid == gene_id or gname == gene_name for gid, gname, _ in overlaps)


def estimate_amplicon_bounds(anchor_hit, declared_size):
    """Estimate the amplicon's reference bounds from a single confidently-placed
    primer plus the declared amplicon size, assuming a normal PCR product:
    whichever primer sits on '+' is upstream/left of the amplicon and '-' is
    downstream/right -- true regardless of which primer (fwd or rev) it is."""
    if anchor_hit.strand == "+":
        start = anchor_hit.start
        end = start + declared_size - 1
    else:
        end = anchor_hit.end
        start = end - declared_size + 1
    return max(1, start), end


def pick_estimation_anchor(fwd_result, rev_result):
    """When fwd_primer/rev_primer failed to produce a normal amplicon (one
    unmapped, or mapped to different chromosomes -- e.g. the other one hit a
    repeat and BWA's tie-break landed miles away), decide whether either
    single primer is trustworthy enough to anchor a declared-size estimate.
    Requires the candidate to actually overlap the declared gene; returns
    (result, seq_type) or (None, None) if it's not safe to guess (neither
    mapped, or both did but neither/both overlap the gene -- no way to tell
    which one to trust)."""
    have_fwd, have_rev = fwd_result is not None, rev_result is not None
    if have_fwd and not have_rev:
        return (fwd_result, "fwd_primer") if fwd_result["gene_match"] else (None, None)
    if have_rev and not have_fwd:
        return (rev_result, "rev_primer") if rev_result["gene_match"] else (None, None)
    if have_fwd and have_rev:
        fwd_ok, rev_ok = fwd_result["gene_match"], rev_result["gene_match"]
        if fwd_ok and not rev_ok:
            return fwd_result, "fwd_primer"
        if rev_ok and not fwd_ok:
            return rev_result, "rev_primer"
    return None, None


def primer_orientation_ok(fwd_hit, rev_hit):
    """A normal PCR product needs fwd/rev on opposite strands, with whichever
    one is '+' sitting upstream (lower coordinate) of the '-' one -- true
    regardless of which primer happens to be the gene's sense strand."""
    if fwd_hit.chrom != rev_hit.chrom or fwd_hit.strand == rev_hit.strand:
        return False
    if fwd_hit.strand == "+":
        return fwd_hit.start < rev_hit.start
    return rev_hit.start < fwd_hit.start


def find_tied_alternate(current_hit, candidates, matches):
    """Among candidates tied with current_hit at the same edit distance (i.e.
    equally good, per BWA), return the first (sorted for determinism)
    satisfying `matches`, else None."""
    tied = sorted((h for h in candidates if h.nm == current_hit.nm), key=lambda h: (h.chrom, h.start))
    return next((h for h in tied if matches(h)), None)


def compute_status_bits(gene_match, overlaps, n_alt, n_best_hits, high_copy_threshold, reassigned_from):
    status_bits = []
    if gene_match:
        status_bits.append("OK")
    elif overlaps:
        status_bits.append("GENE_MISMATCH")
    else:
        status_bits.append("NO_GENE_OVERLAP")
    if n_alt > 0 or (n_best_hits or 1) > 1:
        status_bits.append("MULTI_MAPPED")
    if (n_best_hits or 0) > high_copy_threshold:
        # a huge tied-best-hit count (BWA's X0) is the signature of a sequence that
        # overlaps a common repeat (e.g. an Alu element) -- at that scale, bwa samse's
        # XA tag gets silently suppressed entirely (see the two-alignments note above),
        # so the gene-overlap/orientation tie-correction has no alternates left to search
        # and the reported position is essentially an arbitrary pick among thousands.
        status_bits.append("HIGH_COPY_NUMBER")
    if reassigned_from:
        status_bits.append("REASSIGNED_PRIMARY")
    return status_bits


def apply_reassignment(sample_id, seq_type, position_source, old_primary, match, candidates, rout):
    """Swap `match` in as the new primary, demoting the old primary into the
    alt list. Logs the swap to rout. Returns (new_primary, new_candidates,
    reassigned_from_str)."""
    reassigned_from = "%s:%d-%d(%s)" % (old_primary.chrom, old_primary.start, old_primary.end, old_primary.strand)
    demoted = Hit(
        old_primary.chrom, old_primary.start, old_primary.end, old_primary.strand,
        None, old_primary.nm, "XA", old_primary.cigar,
    )
    new_candidates = [h for h in candidates if h is not match] + [demoted]
    new_primary = Hit(
        match.chrom, match.start, match.end, match.strand,
        old_primary.mapq, match.nm, "primary", match.cigar,
    )
    rout.write("\t".join([
        sample_id, seq_type, position_source,
        old_primary.chrom, str(old_primary.start), str(old_primary.end), old_primary.strand,
        new_primary.chrom, str(new_primary.start), str(new_primary.end), new_primary.strand, str(new_primary.nm),
    ]) + "\n")
    return new_primary, new_candidates, reassigned_from


def extract_fasta_seq(ref_path, fai_path, chrom, start, end):
    """samtools faidx ref chrom:start-end -> sequence string (no header, no newlines)."""
    region = "%s:%d-%d" % (chrom, start, end)
    out = subprocess.run(
        ["samtools", "faidx", "--fai-idx", fai_path, ref_path, region],
        capture_output=True, text=True, check=True,
    )
    lines = out.stdout.splitlines()
    return "".join(lines[1:]) if lines else ""


def resolve_position(sample_id, seq_type, hits_exact, hits_mm):
    """Pick which alignment run supplies this read's reported position.

    The exact-match run is always preferred (its ties are naturally few, so
    they're never hidden the way a mismatch-tolerant search's XA tag can be).
    Falls back to the mismatch-tolerant run only if there's no exact hit at
    all (e.g. a real SNP/indel at the true target site).

    Returns (position_source, primary_hit, candidate_alt_hits) or (None, None, None).
    """
    qr_exact = hits_exact.get((sample_id, seq_type)) if hits_exact else None
    if qr_exact and qr_exact.primary:
        return "exact", qr_exact
    qr_mm = hits_mm.get((sample_id, seq_type))
    if qr_mm and qr_mm.primary:
        return "mismatch_tolerant", qr_mm
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--sam", required=True,
                     help="mismatch-tolerant alignment (BWA_N/MAX_XA) -- used for off_targets.tsv, and as a "
                          "position fallback for reads with no exact (0-mismatch) hit anywhere")
    ap.add_argument("--exact-sam", default=None,
                     help="exact-match-only alignment (bwa aln -n 0). Source of truth for position/gene-"
                          "overlap/amplicon-extraction when given; omit to use --sam for everything (legacy "
                          "single-alignment mode)")
    ap.add_argument("--ref", required=True, help="genome FASTA")
    ap.add_argument("--fai", default=None,
                     help="path to the .fai index (default: <ref>.fai, i.e. alongside the fasta)")
    ap.add_argument("--gene-bed", required=True)
    ap.add_argument("--gene-window", type=int, default=0,
                     help="pad each gene's start/end by this many bp before checking guide/primer overlap "
                          "-- e.g. primers are often designed just outside the gene body/UTR, which would "
                          "otherwise report a real on-target primer as NO_GENE_OVERLAP/GENE_MISMATCH. Also "
                          "widens the overlapping_genes/overlapping_gene_names columns and off-target gene "
                          "matching accordingly. Does not affect --exon-bed overlap.")
    ap.add_argument("--exon-bed", default=None,
                     help="optional exon BED (from gtf_to_genebed.py) for the overlapping_exons column")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-mismatch", type=int, default=2,
                     help="edit-distance threshold below which an extra hit is reported as a high-quality off-target")
    ap.add_argument("--high-copy-threshold", type=int, default=50,
                     help="flag a guide/primer as HIGH_COPY_NUMBER in the status column if its tied-best-hit "
                          "count (BWA's X0 tag, n_best_hits) exceeds this -- a strong signal it overlaps a "
                          "common repeat element (e.g. an Alu), where the true off-target list is likely "
                          "incomplete (bwa samse silently drops the XA tag once total hits are too numerous)")
    ap.add_argument("--reassign-primary", action="store_true",
                     help="also write a corrected SAM/BAM (<exact-sam>.reassigned.sam) reflecting whichever "
                          "gene-matching tied position was used in the report, for viewing in IGV. The report "
                          "itself always uses the corrected position regardless of this flag; this only "
                          "controls whether a corrected alignment file is also produced.")
    args = ap.parse_args()

    fai_path = args.fai or (args.ref + ".fai")
    if not os.path.exists(fai_path):
        sys.exit("ERROR: %s not found. Run setup_reference.sh first." % fai_path)

    reports_dir = os.path.join(args.outdir, "reports")
    amplicons_dir = os.path.join(args.outdir, "amplicons")
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(amplicons_dir, exist_ok=True)

    rows = dict(read_rows(args.tsv))
    hits_mm = parse_sam(args.sam)
    hits_exact = parse_sam(args.exact_sam) if args.exact_sam else None
    gene_index = IntervalIndex(  # gene_id, gene_name, strand
        args.gene_bed, n_extra=3, label="GeneIndex", padding=args.gene_window
    )
    exon_index = (
        IntervalIndex(args.exon_bed, n_extra=5, label="ExonIndex")  # gene_id, gene_name, transcript_id, exon_number, strand
        if args.exon_bed else None
    )

    all_hit_sources = list(hits_mm.values()) + (list(hits_exact.values()) if hits_exact else [])
    mapped_chroms = {normalize_chrom(qr.primary.chrom) for qr in all_hit_sources if qr.primary}
    if mapped_chroms and not (mapped_chroms & set(gene_index.by_chrom)):
        sys.stderr.write(
            "WARNING: none of the chromosomes in the alignments (%s) match any "
            "chromosome in --gene-bed (%s). Every hit will be reported as "
            "NO_GENE_OVERLAP even where IGV shows a real overlap. This usually "
            "means the gene BED and the genome FASTA use different chromosome "
            "naming (e.g. 'chr1' vs '1', or 'chrM' vs 'MT') and were not built "
            "from a matching pair -- rebuild GENE_BED from a GTF for the exact "
            "REF_FASTA you aligned against.\n"
            % (", ".join(sorted(mapped_chroms)[:5]), ", ".join(sorted(gene_index.by_chrom)[:5]))
        )
    if exon_index is not None and mapped_chroms and not (mapped_chroms & set(exon_index.by_chrom)):
        sys.stderr.write(
            "WARNING: none of the chromosomes in the alignments (%s) match any "
            "chromosome in --exon-bed (%s). Every hit will show 'intronic' even "
            "where IGV shows the guide/primer sitting inside an exon -- rebuild "
            "EXON_BED (and GENE_BED) from a GTF for the exact REF_FASTA you "
            "aligned against.\n"
            % (", ".join(sorted(mapped_chroms)[:5]), ", ".join(sorted(exon_index.by_chrom)[:5]))
        )

    validation_path = os.path.join(reports_dir, "validation_report.tsv")
    amplicon_report_path = os.path.join(reports_dir, "amplicon_report.tsv")
    offtarget_path = os.path.join(reports_dir, "off_targets.tsv")
    crispresso_path = os.path.join(reports_dir, "crispresso_amplicons.txt")
    reassigned_log_path = os.path.join(reports_dir, "reassigned_primaries.tsv")
    combined_fasta_path = os.path.join(amplicons_dir, "all_amplicons.fasta")

    n_ok, n_flagged, n_unmapped = 0, 0, 0
    n_amplicons_ok = 0
    n_estimated = 0
    n_crispresso = 0
    n_reassigned = 0
    reassignments = {}  # (sample_id, seq_type) -> dict of SAM fields, for write_reassigned_sam()

    with open(validation_path, "w") as vout, \
         open(amplicon_report_path, "w") as aout, \
         open(offtarget_path, "w") as oout, \
         open(combined_fasta_path, "w") as fout, \
         open(crispresso_path, "w") as cout, \
         open(reassigned_log_path, "w") as rout:

        vout.write("\t".join([
            "sample_id", "gene_name_declared", "gene_id_declared", "seq_type",
            "chrom", "start", "end", "strand", "mapq", "edit_distance",
            "n_best_hits", "n_suboptimal_hits", "n_alt_hits", "overlapping_genes",
            "overlapping_gene_names", "overlapping_exons", "position_source",
            "reassigned_from", "gene_declared_in_overlaps", "distance_from_amplicon_start",
            "status",
        ]) + "\n")

        rout.write("\t".join([
            "sample_id", "seq_type", "position_source", "old_chrom", "old_start", "old_end",
            "old_strand", "new_chrom", "new_start", "new_end", "new_strand", "edit_distance",
        ]) + "\n")

        aout.write("\t".join([
            "sample_id", "gene_name_declared", "gene_id_declared", "chrom",
            "amplicon_start", "amplicon_end", "amplicon_length",
            "declared_amplicon_size", "length_diff", "fwd_strand", "rev_strand",
            "orientation_ok", "fasta_file", "status",
        ]) + "\n")

        oout.write("\t".join([
            "sample_id", "gene_name_declared", "seq_type", "chrom", "start", "end",
            "strand", "edit_distance", "source", "overlapping_genes",
        ]) + "\n")

        for sample_id, row in rows.items():
            declared_gene_name = row["gene_name"]
            declared_gene_id = row["gene_id"]

            per_type_position = {}  # seq_type -> final Hit (or None)
            per_type_result = {}    # seq_type -> dict of everything needed to write its row, or None if unmapped

            # --- pass 1: resolve each guide/primer's position (incl. tie-correction) ---
            for seq_type in SEQ_TYPES:
                position_source, qr = resolve_position(sample_id, seq_type, hits_exact, hits_mm)
                per_type_position[seq_type] = None

                if qr is None:
                    n_unmapped += 1
                    per_type_result[seq_type] = None
                    continue

                p = qr.primary
                candidates = qr.alt
                overlaps = gene_index.overlapping(p.chrom, p.start, p.end)
                gene_match = declared_gene_in_overlaps(overlaps, declared_gene_id, declared_gene_name)

                # The reported position always uses the best gene-matching tie available
                # (cheap now: the exact-match run's ties are few, never hidden from XA).
                # --reassign-primary only controls whether a corrected BAM is also written.
                reassigned_from = ""
                if not gene_match and candidates:
                    match = find_tied_alternate(
                        p, candidates,
                        lambda h: declared_gene_in_overlaps(
                            gene_index.overlapping(h.chrom, h.start, h.end), declared_gene_id, declared_gene_name
                        ),
                    )
                    if match is not None:
                        p, candidates, reassigned_from = apply_reassignment(
                            sample_id, seq_type, position_source, p, match, candidates, rout
                        )
                        overlaps = gene_index.overlapping(p.chrom, p.start, p.end)
                        gene_match = True
                        n_reassigned += 1
                        if position_source == "exact":
                            reassignments[(sample_id, seq_type)] = {
                                "chrom": p.chrom, "pos": p.start, "strand": p.strand,
                                "cigar": p.cigar, "nm": p.nm, "xa": format_xa(candidates),
                            }

                per_type_position[seq_type] = p
                exon_overlaps = exon_index.overlapping(p.chrom, p.start, p.end) if exon_index else []
                n_alt = len(candidates)

                status_bits = compute_status_bits(
                    gene_match, overlaps, n_alt, qr.n_best_hits, args.high_copy_threshold, reassigned_from
                )
                if status_bits[0] == "OK":
                    n_ok += 1
                else:
                    n_flagged += 1

                per_type_result[seq_type] = {
                    "p": p, "qr": qr, "candidates": candidates, "overlaps": overlaps,
                    "gene_match": gene_match, "exon_overlaps": exon_overlaps, "n_alt": n_alt,
                    "status_bits": status_bits, "reassigned_from": reassigned_from,
                    "position_source": position_source,
                }

            # --- orientation fix: fwd/rev tied onto the same strand (or the wrong
            # order) is the same class of bad BWA tie-break as the gene-mismatch
            # case above, just caught by checking the *pair* instead of one read in
            # isolation. If either primer has an equally-good (same edit distance)
            # alternate that would put the pair back in correct orientation AND
            # still overlaps the declared gene, use it -- tried on fwd first, then
            # rev, and only one of the two is ever swapped.
            fwd_result, rev_result = per_type_result.get("fwd_primer"), per_type_result.get("rev_primer")
            if fwd_result is not None and rev_result is not None:
                fwd_p, rev_p = fwd_result["p"], rev_result["p"]
                if not primer_orientation_ok(fwd_p, rev_p):
                    fix = find_tied_alternate(
                        fwd_p, fwd_result["candidates"],
                        lambda h: primer_orientation_ok(h, rev_p) and declared_gene_in_overlaps(
                            gene_index.overlapping(h.chrom, h.start, h.end), declared_gene_id, declared_gene_name
                        ),
                    )
                    fixed_result, fixed_seq_type = fwd_result, "fwd_primer"
                    if fix is None:
                        fix = find_tied_alternate(
                            rev_p, rev_result["candidates"],
                            lambda h: primer_orientation_ok(fwd_p, h) and declared_gene_in_overlaps(
                                gene_index.overlapping(h.chrom, h.start, h.end), declared_gene_id, declared_gene_name
                            ),
                        )
                        fixed_result, fixed_seq_type = rev_result, "rev_primer"
                    if fix is not None:
                        old_p = fixed_result["p"]
                        new_p, new_candidates, reassigned_from = apply_reassignment(
                            sample_id, fixed_seq_type, fixed_result["position_source"],
                            old_p, fix, fixed_result["candidates"], rout,
                        )
                        overlaps = gene_index.overlapping(new_p.chrom, new_p.start, new_p.end)
                        exon_overlaps = exon_index.overlapping(new_p.chrom, new_p.start, new_p.end) if exon_index else []
                        gene_match = declared_gene_in_overlaps(overlaps, declared_gene_id, declared_gene_name)
                        n_alt = len(new_candidates)

                        old_category = fixed_result["status_bits"][0]
                        status_bits = compute_status_bits(
                            gene_match, overlaps, n_alt, fixed_result["qr"].n_best_hits,
                            args.high_copy_threshold, reassigned_from,
                        )
                        if old_category != status_bits[0]:
                            n_ok += 1 if status_bits[0] == "OK" else -1
                            n_flagged += 1 if status_bits[0] != "OK" else -1

                        fixed_result.update({
                            "p": new_p, "candidates": new_candidates, "overlaps": overlaps,
                            "gene_match": gene_match, "exon_overlaps": exon_overlaps, "n_alt": n_alt,
                            "status_bits": status_bits, "reassigned_from": reassigned_from,
                        })
                        per_type_position[fixed_seq_type] = new_p
                        n_reassigned += 1
                        if fixed_result["position_source"] == "exact":
                            reassignments[(sample_id, fixed_seq_type)] = {
                                "chrom": new_p.chrom, "pos": new_p.start, "strand": new_p.strand,
                                "cigar": new_p.cigar, "nm": new_p.nm, "xa": format_xa(new_candidates),
                            }

            # amplicon bounds, needed for distance_from_amplicon_start below. Ignores
            # strand/direction entirely -- just the leftmost (5') reference coordinate.
            amp_chrom = amp_start = amp_end = None
            _fwd, _rev = per_type_position["fwd_primer"], per_type_position["rev_primer"]
            if _fwd is not None and _rev is not None and _fwd.chrom == _rev.chrom:
                amp_chrom = _fwd.chrom
                amp_start = min(_fwd.start, _rev.start)
                amp_end = max(_fwd.end, _rev.end)

            # --- pass 2: write validation_report.tsv + off_targets.tsv rows ---
            for seq_type in SEQ_TYPES:
                result = per_type_result[seq_type]
                if result is None:
                    vout.write("\t".join([
                        sample_id, declared_gene_name, declared_gene_id, seq_type,
                        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "UNMAPPED",
                    ]) + "\n")
                    continue

                p, qr, candidates = result["p"], result["qr"], result["candidates"]
                overlaps, gene_match = result["overlaps"], result["gene_match"]
                exon_overlaps, n_alt = result["exon_overlaps"], result["n_alt"]
                status_bits, reassigned_from = result["status_bits"], result["reassigned_from"]
                position_source = result["position_source"]

                distance_from_amplicon_start = ""
                if amp_chrom is not None and p.chrom == amp_chrom:
                    midpoint = (p.start + p.end) // 2
                    distance_from_amplicon_start = str(midpoint - amp_start)

                vout.write("\t".join([
                    sample_id, declared_gene_name, declared_gene_id, seq_type,
                    p.chrom, str(p.start), str(p.end), p.strand,
                    str(p.mapq) if p.mapq is not None else "",
                    str(p.nm) if p.nm is not None else "",
                    str(qr.n_best_hits) if qr.n_best_hits is not None else "",
                    str(qr.n_suboptimal_hits) if qr.n_suboptimal_hits is not None else "",
                    str(n_alt), gene_list_str(overlaps), gene_names_str(overlaps),
                    exon_list_str(exon_overlaps, bool(overlaps), exon_index is not None),
                    position_source, reassigned_from, str(gene_match),
                    distance_from_amplicon_start, "+".join(status_bits),
                ]) + "\n")

                # off-target reporting always draws from the mismatch-tolerant run
                # (more sensitive), excluding whichever position was used above
                mm_qr = hits_mm.get((sample_id, seq_type))
                offtarget_hits = ([mm_qr.primary] + mm_qr.alt) if mm_qr and mm_qr.primary else ([p] + candidates)
                for h in offtarget_hits:
                    on_target_pos = (h.chrom == p.chrom and h.start == p.start)
                    h_overlaps = overlaps if on_target_pos else gene_index.overlapping(h.chrom, h.start, h.end)
                    on_target = on_target_pos or declared_gene_in_overlaps(h_overlaps, declared_gene_id, declared_gene_name)
                    if on_target:
                        continue
                    if h.nm is None or h.nm > args.max_mismatch:
                        continue
                    oout.write("\t".join([
                        sample_id, declared_gene_name, seq_type, h.chrom, str(h.start), str(h.end),
                        h.strand, str(h.nm), h.source, gene_list_str(h_overlaps),
                    ]) + "\n")

            # --- amplicon extraction ---
            fwd = per_type_position["fwd_primer"]
            rev = per_type_position["rev_primer"]
            declared_size = row.get("amplicon_size", "")
            fwd_strand_out = fwd.strand if fwd is not None else ""
            rev_strand_out = rev.strand if rev is not None else ""

            ext_chrom = ext_start = ext_end = None
            orientation_ok = None
            estimated_from = None

            if fwd is not None and rev is not None and amp_chrom is not None:
                ext_chrom, ext_start, ext_end = amp_chrom, amp_start, amp_end
                orientation_ok = primer_orientation_ok(fwd, rev)
            else:
                # one primer unmapped, or the pair landed on different chromosomes
                # (e.g. the other primer hit a repeat -- see HIGH_COPY_NUMBER above)
                # -- if exactly one primer is confidently on the declared gene, fall
                # back to estimating the amplicon from it plus the declared size,
                # rather than giving up on the sequence entirely.
                failure_status = "PRIMER_UNMAPPED" if (fwd is None or rev is None) else "PRIMERS_ON_DIFFERENT_CHROMS"
                declared_size_int = None
                if declared_size:
                    try:
                        declared_size_int = int(declared_size)
                    except ValueError:
                        pass
                if declared_size_int:
                    anchor_result, estimated_from = pick_estimation_anchor(
                        per_type_result.get("fwd_primer"), per_type_result.get("rev_primer")
                    )
                    if anchor_result is not None:
                        anchor_hit = anchor_result["p"]
                        ext_start, ext_end = estimate_amplicon_bounds(anchor_hit, declared_size_int)
                        ext_chrom = anchor_hit.chrom
                if ext_chrom is None:
                    aout.write("\t".join([
                        sample_id, declared_gene_name, declared_gene_id, "", "", "", "",
                        declared_size, "", fwd_strand_out, rev_strand_out, "", "", failure_status,
                    ]) + "\n")
                    continue

            ext_len = ext_end - ext_start + 1

            try:
                seq = extract_fasta_seq(args.ref, fai_path, ext_chrom, ext_start, ext_end)
            except subprocess.CalledProcessError as e:
                aout.write("\t".join([
                    sample_id, declared_gene_name, declared_gene_id, ext_chrom,
                    str(ext_start), str(ext_end), str(ext_len), declared_size, "",
                    fwd_strand_out, rev_strand_out, str(orientation_ok), "", "FAIDX_ERROR",
                ]) + "\n")
                sys.stderr.write("ERROR extracting %s: %s\n" % (sample_id, e.stderr))
                continue

            fasta_name = "%s.fa" % sample_id
            fasta_path = os.path.join(amplicons_dir, fasta_name)
            extra_hdr = (
                "estimated_from=%s" % estimated_from if estimated_from
                else "orientation_ok=%s" % orientation_ok
            )
            header = ">%s|%s|%s:%d-%d|len=%d|declared=%s|%s" % (
                sample_id, declared_gene_name, ext_chrom, ext_start, ext_end, ext_len,
                declared_size, extra_hdr,
            )
            with open(fasta_path, "w") as sfout:
                sfout.write(header + "\n" + seq + "\n")
            fout.write(header + "\n" + seq + "\n")

            # CRISPRessoPooled amplicon description file: amplicon_name, amplicon_seq,
            # guide_seq, expected_hdr_amplicon_seq, coding_seq (5 cols, no header, "NA"
            # for anything not applicable/available). guide_seq is pulled from the
            # reference at each guide's aligned coordinates (not the designed oligo),
            # so it's guaranteed to be an exact substring of amplicon_seq even where
            # bwa aln allowed a mismatch.
            guide_seqs = []
            for guide_type in ("guide1", "guide2", "guide3"):
                g = per_type_position.get(guide_type)
                if g is None:
                    continue
                if g.chrom != ext_chrom or not (ext_start <= g.start and g.end <= ext_end):
                    continue
                try:
                    guide_seqs.append(extract_fasta_seq(args.ref, fai_path, g.chrom, g.start, g.end))
                except subprocess.CalledProcessError:
                    pass
            cout.write("\t".join([
                sample_id, seq, ",".join(guide_seqs) if guide_seqs else "NA", "NA", "NA",
            ]) + "\n")
            n_crispresso += 1

            length_diff = ""
            if declared_size:
                try:
                    length_diff = str(ext_len - int(declared_size))
                except ValueError:
                    pass

            if estimated_from:
                status = "AMPLICON_ESTIMATED_FROM_%s" % estimated_from.upper()
                n_estimated += 1
                n_flagged += 1
            else:
                status = "OK" if orientation_ok else "UNEXPECTED_ORIENTATION"
                if status == "OK":
                    n_amplicons_ok += 1
                else:
                    n_flagged += 1

            aout.write("\t".join([
                sample_id, declared_gene_name, declared_gene_id, ext_chrom,
                str(ext_start), str(ext_end), str(ext_len), declared_size, length_diff,
                fwd_strand_out, rev_strand_out, str(orientation_ok) if orientation_ok is not None else "",
                fasta_name, status,
            ]) + "\n")

    sys.stderr.write(
        "build_report: %d guide/primer hits OK, %d flagged, %d unmapped; %d amplicons extracted OK "
        "(%d estimated from a single primer + declared size); %d rows written to CRISPResso amplicons "
        "file; %d positions corrected to a gene-matching tie\n"
        % (n_ok, n_flagged, n_unmapped, n_amplicons_ok, n_estimated, n_crispresso, n_reassigned)
    )
    if args.reassign_primary:
        if args.exact_sam:
            reassigned_sam_path = os.path.join(
                os.path.dirname(args.exact_sam),
                os.path.splitext(os.path.basename(args.exact_sam))[0] + ".reassigned.sam",
            )
            write_reassigned_sam(args.exact_sam, reassigned_sam_path, reassignments)
            sys.stderr.write("Reassigned SAM written to %s\n" % reassigned_sam_path)
        else:
            sys.stderr.write("WARNING: --reassign-primary given without --exact-sam; no corrected SAM written.\n")
    sys.stderr.write("Reports written to %s\n" % reports_dir)
    sys.stderr.write("Amplicon FASTAs written to %s\n" % amplicons_dir)


if __name__ == "__main__":
    main()
