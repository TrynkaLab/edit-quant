#!/usr/bin/env python3
"""Match each fastq-manifest sample to the amplicons designed for its project.

Splits crispresso_amplicons.txt by the amplicon manifest's project_name,
writing amplicons_by_project/<project_name>.amplicons.txt.

CRISPRessoPooled requires unique amplicon sequences, and wells sharing an
amplicon within a project (technical replicates) can't be told apart in a
sample's pooled reads anyway. So each group of identical sequences is merged
into one amplicon, named after all its wells (e.g. IL2RA_G4-D8-B9_<plate>_<project>),
with the union of the wells' guides. reports/amplicon_replicates.tsv maps every
well's amplicon to the merged CRISPResso amplicon name.

Then writes reports/sample_check.tsv with one row per fastq-manifest sample;
main.nf runs only the samples with status OK, and warns about and skips the rest:

  OK                               project has at least one extracted amplicon
  NO_PROJECT_IN_AMPLICON_MANIFEST  project_name has no rows in the amplicon manifest
  NO_AMPLICONS_EXTRACTED           project is in the manifest, but none of its
                                   amplicons could be extracted from the reference
"""
import argparse
import csv
import os
import sys

from pipeline_common import read_rows, sanitize


def unique(values):
    """Unique values, in first-seen order."""
    return list(dict.fromkeys(values))


def merge_replicates(group, manifest_rows):
    """Merge crispresso amplicon rows sharing one amplicon_seq into a single row.

    The name keeps the sample_id layout (gene_well_plate_project), joining the
    differing values of each part with '-'.
    """
    if len(group) == 1:
        return group[0], False
    rows = [manifest_rows[fields[0]] for fields in group]
    name = "_".join(
        "-".join(unique(sanitize(r[col]) for r in rows))
        for col in ("gene_name", "well_position", "plate_name", "project_name")
    )
    guide_sets = [tuple(g for g in fields[2].split(",") if g != "NA") for fields in group]
    guides = unique(g for gs in guide_sets for g in gs)
    merged = [name, group[0][1], ",".join(guides) if guides else "NA"] + group[0][3:]
    return merged, len(set(guide_sets)) > 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amplicon-manifest", required=True)
    ap.add_argument("--crispresso-amplicons", required=True)
    ap.add_argument("--fastq-manifest", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    # amplicon name -> manifest row, from the design manifest
    manifest_rows = dict(read_rows(args.amplicon_manifest))
    designed_projects = {row["project_name"] for row in manifest_rows.values()}

    # project_name -> {amplicon_seq: [rows]}, grouping technical replicates
    amplicons = {}
    with open(args.crispresso_amplicons) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if not fields[0]:
                continue
            project = manifest_rows[fields[0]]["project_name"]
            amplicons.setdefault(project, {}).setdefault(fields[1], []).append(fields)

    reports_dir = os.path.join(args.outdir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    amplicon_dir = os.path.join(args.outdir, "amplicons_by_project")
    os.makedirs(amplicon_dir, exist_ok=True)
    amplicon_files = {}
    with open(os.path.join(reports_dir, "amplicon_replicates.tsv"), "w", newline="") as rep:
        rep_writer = csv.writer(rep, delimiter="\t", lineterminator="\n")
        rep_writer.writerow(["project_name", "amplicon_name", "crispresso_amplicon_name",
                             "n_replicates", "guides_differ"])
        for project, by_seq in amplicons.items():
            name = "%s.amplicons.txt" % sanitize(project)
            if name in (n for n, _ in amplicon_files.values()):
                sys.exit("ERROR: two project_name values map to the same file name %s" % name)
            with open(os.path.join(amplicon_dir, name), "w") as out:
                for group in by_seq.values():
                    merged, guides_differ = merge_replicates(group, manifest_rows)
                    if guides_differ:
                        print("WARNING: replicates %s share an amplicon but have different guides - "
                              "quantifying against all of them" % [f[0] for f in group], file=sys.stderr)
                    out.write("\t".join(merged) + "\n")
                    for fields in group:
                        rep_writer.writerow([project, fields[0], merged[0], len(group), guides_differ])
            amplicon_files[project] = (name, len(by_seq))

    with open(args.fastq_manifest, newline="", encoding="utf-8-sig") as fh, \
            open(os.path.join(reports_dir, "sample_check.tsv"), "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "project_name", "donor", "fastq_1", "fastq_2",
                         "n_amplicons", "amplicons_file", "status"])
        for row in csv.DictReader(fh, delimiter="\t"):
            row = {k: (v or "").strip() for k, v in row.items()}
            if not any(row.values()):
                continue
            project = row["project_name"]
            if project not in designed_projects:
                status = "NO_PROJECT_IN_AMPLICON_MANIFEST"
            elif project not in amplicon_files:
                status = "NO_AMPLICONS_EXTRACTED"
            else:
                status = "OK"
            name, n = amplicon_files.get(project, ("NA", 0))
            writer.writerow([row["sample"], project, row.get("donor") or "none",
                             row["fastq_1"], row["fastq_2"], n, name, status])


if __name__ == "__main__":
    main()
