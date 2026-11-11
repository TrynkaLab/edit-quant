#!/usr/bin/env python3
"""Shared packing for the amplicon viewer: CRISPResso per-amplicon tables -> compact arrays.

Used by pack_viewer_data.py (one task per sample, in the pipeline) and by
build_amplicon_viewer.py (standalone, any number of samples at once).

Per CRISPResso_on_<amplicon> folder it reads:
  Nucleotide_frequency_table.txt   reference sequence (header row) + A/C/G/T/N/- counts
  Modification_count_vectors.txt   insertion / deletion / substitution counts
  CRISPResso2_info.json            guide intervals, cut points, quantification window
All three survive --suppress_plots: CRISPRessoCORE writes them above the first in-loop
`if not args.suppress_plots` guard, which is why the viewer costs no CRISPResso runtime.

Counts become one byte per position per series (fraction of aligned reads scaled to
0-255, ~0.4% resolution), series-major. That is ~1.8 KB per amplicon gzipped, so a
30-sample run of 500 amplicons is tens of MB rather than the >100 MB the same data
would take as JSON numbers.
"""

import json
import os

# Byte layout per amplicon, matching SER in assets/amplicon_viewer.js
SERIES = ["A", "C", "G", "T", "N", "-", "ins", "sub"]

# A deletion block is called when at least this fraction of reads are deleted, and its
# ends must fall within BLOCK_CUT_SLOP of a cut site to count as inter-guide excision.
BLOCK_MIN_FRAC = 0.10
BLOCK_CUT_SLOP = 10
BLOCK_GAP_TOLERANCE = 5


def read_matrix(path):
    """Tab-separated table with a row label in column 0 -> (header_cols, {label: [float]})."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            out[parts[0]] = [float(x) for x in parts[1:]]
    return header[1:], out


def guide_strand(interval, cut, guide_len, qw_center):
    """CRISPResso places the cut at start+qw_center+len-1 for a forward match and at
    start-qw_center-1 for a reverse-complement one. With the SpCas9 default
    (qw_center=-3, 20mer) that is start+16 vs start+2, so which side of the interval the
    cut falls on identifies the strand without looking at the PAM. Pass the same
    quantification_window_center the run used, or Cas12a designs (qw_center=1) invert."""
    start, _ = interval
    fwd = start + qw_center + guide_len - 1
    rev = start - qw_center - 1
    return "+" if abs(cut - fwd) <= abs(cut - rev) else "-"


def deletion_block(frac, cut_positions):
    """Extent of the dominant deletion, and whether its ends sit on cut sites.

    In a multi-guide design the expected outcome is excision between two cuts, which
    shows as a contiguous block ending on cut sites. Taking the longest run rather than
    first-to-last matters: a single stray position over the threshold would otherwise
    stretch the block across the whole amplicon and flag it as unexplained.
    """
    peak = max(frac) if frac else 0.0
    if peak < BLOCK_MIN_FRAC:
        return None, True
    thresh = max(BLOCK_MIN_FRAC, 0.5 * peak)
    runs, start, end, gap = [], None, None, 0
    for i, v in enumerate(frac):
        if v >= thresh:
            if start is None:
                start = i
            end, gap = i, 0
        elif start is not None:
            gap += 1
            if gap > BLOCK_GAP_TOLERANCE:
                runs.append((start, end))
                start = None
    if start is not None:
        runs.append((start, end))
    if not runs:
        return None, True
    block = max(runs, key=lambda r: r[1] - r[0])
    if not cut_positions:
        return list(block), False
    near = lambda pos: min(abs(pos - c) for c in cut_positions)
    ok = near(block[0]) <= BLOCK_CUT_SLOP and near(block[1]) <= BLOCK_CUT_SLOP
    return list(block), ok


def load_amplicon(folder, sample, qw_center):
    """One CRISPResso_on_<amplicon> folder -> record dict, or None if unusable."""
    nuc = os.path.join(folder, "Nucleotide_frequency_table.txt")
    mod = os.path.join(folder, "Modification_count_vectors.txt")
    info_path = os.path.join(folder, "CRISPResso2_info.json")
    if not all(os.path.exists(p) for p in (nuc, mod, info_path)):
        return None

    refseq_cols, counts = read_matrix(nuc)
    _, mods = read_matrix(mod)
    length = len(refseq_cols)
    if length == 0:
        return None

    # Column 0 of every row sums to the aligned-read total for this amplicon.
    total = sum(counts[b][0] for b in "ACGTN-" if b in counts)
    if total <= 0:
        return None

    with open(info_path) as fh:
        info = json.load(fh)
    refs = info.get("results", {}).get("refs", {})
    if not refs:
        return None
    # Pooled runs quantify against a single "Reference" per amplicon folder.
    ref = refs.get("Reference") or next(iter(refs.values()))

    refseq = ref.get("sequence") or "".join(refseq_cols)
    intervals = ref.get("sgRNA_intervals") or []
    cuts = ref.get("sgRNA_cut_points") or []
    seqs = ref.get("sgRNA_sequences") or []
    guides = []
    for i, iv in enumerate(intervals):
        cut = cuts[i] if i < len(cuts) else iv[0]
        seq = seqs[i] if i < len(seqs) else ""
        guides.append({
            "seq": seq,
            "start": int(iv[0]),
            "end": int(iv[1]),
            "cut": int(cut),
            "strand": guide_strand(iv, cut, len(seq) or 20, qw_center),
        })

    qwin = ref.get("include_idxs") or []
    if isinstance(qwin, dict):          # numpy arrays serialise as {_type, value}
        qwin = qwin.get("value", [])
    qwin = [int(x) for x in qwin]

    series = []
    for name in SERIES:
        if name == "ins":
            series.append(mods.get("Insertions", [0.0] * length))
        elif name == "sub":
            series.append(mods.get("Substitutions", [0.0] * length))
        else:
            series.append(counts.get(name, [0.0] * length))

    packed = bytearray()
    for vec in series:
        for i in range(length):
            v = vec[i] / total if i < len(vec) else 0.0
            packed.append(max(0, min(255, round(v * 255))))

    deletions = counts.get("-", [0.0] * length)
    cut_positions = [g["cut"] for g in guides if 0 <= g["cut"] < length]
    cut_del = max((deletions[c] / total for c in cut_positions), default=0.0)
    block, block_ok = deletion_block([d / total for d in deletions], cut_positions)

    modified = 0.0
    quant = os.path.join(folder, "CRISPResso_quantification_of_editing_frequency.txt")
    if os.path.exists(quant):
        with open(quant) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            row = fh.readline().rstrip("\n").split("\t")
        if "Modified%" in head and len(row) == len(head):
            modified = float(row[head.index("Modified%")]) / 100.0

    name = os.path.basename(folder)
    name = name[len("CRISPResso_on_"):] if name.startswith("CRISPResso_on_") else name
    return {
        "sample": sample,
        "name": name,
        "len": length,
        "refseq": refseq[:length],
        "reads": int(round(total)),
        "modified": round(modified, 5),
        "cut_del": round(cut_del, 5),
        "block": block,
        "block_ok": block_ok,
        "guides": guides,
        "qwin": qwin,
        "_packed": bytes(packed),
    }


def sample_name(pooled_dir):
    base = os.path.basename(os.path.normpath(pooled_dir))
    prefix = "CRISPRessoPooled_on_"
    return base[len(prefix):] if base.startswith(prefix) else base


def pack_pooled_dir(pooled_dir, sample=None, qw_center=-3, min_reads=0):
    """One CRISPRessoPooled_on_<sample> dir -> (amplicon records, concatenated bytes).

    Each record carries 'off', its offset into the returned buffer.
    """
    sample = sample or sample_name(pooled_dir)
    amplicons, blob, skipped = [], bytearray(), 0
    folders = sorted(
        os.path.join(pooled_dir, d) for d in os.listdir(pooled_dir)
        if d.startswith("CRISPResso_on_") and os.path.isdir(os.path.join(pooled_dir, d))
    )
    for folder in folders:
        rec = load_amplicon(folder, sample, qw_center)
        if rec is None or rec["reads"] < min_reads:
            skipped += 1
            continue
        packed = rec.pop("_packed")
        rec["off"] = len(blob)
        blob += packed
        amplicons.append(rec)
    return amplicons, bytes(blob), skipped
