#!/usr/bin/env python3
"""Convert the amplicon design TSV into a FASTA of query sequences to align.

Guide sequences are given as RNA (U instead of T) and are converted to DNA.
Primers are taken from fwd_primer/rev_primer (not the adaptor-tagged
final_fwd_primer/final_rev_primer columns), since adaptor sequences are not
genomic and would not align.

Each FASTA record header encodes: sample_id|seq_type|gene_name|gene_id
so build_report.py can recover this metadata after alignment.

Usage: tsv_to_fasta.py input.tsv > queries.fasta
"""
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from pipeline_common import read_rows, rna_to_dna, dna_upper, SEQ_TYPES  # noqa: E402

VALID_BASES = set("ACGTN")


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: tsv_to_fasta.py input.tsv > queries.fasta\n")
        sys.exit(1)
    tsv_path = sys.argv[1]

    n_records = 0
    n_rows = 0
    for sid, row in read_rows(tsv_path):
        n_rows += 1
        for seq_type in SEQ_TYPES:
            raw = row.get(seq_type, "")
            if not raw:
                sys.stderr.write(
                    "WARNING line %s: sample %s has no sequence for %s, skipping\n"
                    % (row["_lineno"], sid, seq_type)
                )
                continue
            seq = rna_to_dna(raw) if seq_type.startswith("guide") else dna_upper(raw)
            bad = set(seq) - VALID_BASES
            if bad:
                sys.stderr.write(
                    "WARNING line %s: sample %s %s has unexpected characters %s\n"
                    % (row["_lineno"], sid, seq_type, sorted(bad))
                )
            header = "%s|%s|%s|%s" % (sid, seq_type, row["gene_name"], row["gene_id"])
            print(">%s" % header)
            print(seq)
            n_records += 1

    sys.stderr.write("tsv_to_fasta: %d rows -> %d query sequences\n" % (n_rows, n_records))


if __name__ == "__main__":
    main()
