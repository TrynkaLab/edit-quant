"""Shared helpers for the amplicon mapping pipeline.

Column names in the input TSV are matched after stripping all non-alphanumeric
characters and lowercasing, so "guide_.1", "guide_1" and "Guide1" are all
treated the same way. This keeps tsv_to_fasta.py and build_report.py in sync
without duplicating the parsing logic.
"""
import csv
import re

# logical field name -> normalized header key it should match
FIELD_KEYS = {
    "project_name": "projectname",
    "plate_name": "platename",
    "well_position": "wellposition",
    "gene_name": "genename",
    "gene_id": "geneid",
    "species": "species",
    "guide1": "guide1",
    "guide2": "guide2",
    "guide3": "guide3",
    "amplicon_size": "ampliconsize",
    "fwd_primer": "fwdprimer",
    "rev_primer": "revprimer",
}

# the five sequences we align per row, in output order
SEQ_TYPES = ["guide1", "guide2", "guide3", "fwd_primer", "rev_primer"]


def _normalize(header):
    return re.sub(r"[^a-z0-9]", "", header.lower())


def build_column_map(fieldnames):
    """Map logical field name -> actual column name in this TSV's header."""
    norm_to_actual = {_normalize(h): h for h in fieldnames}
    colmap = {}
    missing = []
    for logical, key in FIELD_KEYS.items():
        if key in norm_to_actual:
            colmap[logical] = norm_to_actual[key]
        else:
            missing.append(logical)
    if missing:
        raise KeyError(
            "Could not find columns for: %s\nAvailable headers: %s"
            % (", ".join(missing), ", ".join(fieldnames))
        )
    return colmap


def read_rows(tsv_path):
    """Yield (sample_id, row_dict) for each data row, row_dict keyed by logical field name."""
    with open(tsv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header row: %s" % tsv_path)
        colmap = build_column_map(reader.fieldnames)
        for lineno, raw_row in enumerate(reader, start=2):
            row = {}
            for logical, actual in colmap.items():
                val = raw_row.get(actual)
                row[logical] = val.strip() if val is not None else ""
            if not any(row.values()):
                continue  # skip blank lines
            sid = sample_id(row)
            row["_sample_id"] = sid
            row["_lineno"] = lineno
            yield sid, row


def sanitize(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s.strip())


def sample_id(row):
    return "%s_%s_%s" % (
        sanitize(row["project_name"]),
        sanitize(row["plate_name"]),
        sanitize(row["well_position"]),
    )


def rna_to_dna(seq):
    return seq.strip().upper().replace("U", "T")


def dna_upper(seq):
    return seq.strip().upper()
