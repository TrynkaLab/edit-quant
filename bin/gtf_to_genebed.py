#!/usr/bin/env python3
"""Build gene and exon BED files from a GTF.

genes.bed: chrom, start0, end, gene_id, gene_name, strand
exons.bed: chrom, start0, end, gene_id, gene_name, transcript_id, exon_number, strand

Used by build_report.py for the gene/exon-overlap checks. Only "gene" and
"exon" feature rows are used. Handles plain-text or gzipped GTF.

Usage: gtf_to_genebed.py genes.gtf[.gz] genes.bed exons.bed
"""
import sys
import re
import gzip

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: gtf_to_genebed.py genes.gtf[.gz] genes.bed exons.bed\n")
        sys.exit(1)
    gtf_path, genes_out_path, exons_out_path = sys.argv[1:4]

    n_genes = 0
    n_exons = 0
    with opener(gtf_path) as fh, open(genes_out_path, "w") as gout, open(exons_out_path, "w") as eout:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("gene", "exon"):
                continue
            chrom, feature, start, end, strand = fields[0], fields[2], fields[3], fields[4], fields[6]
            attrs = dict(ATTR_RE.findall(fields[8]))
            gene_id = attrs.get("gene_id", "NA")
            gene_name = attrs.get("gene_name", gene_id)
            bed_start = int(start) - 1  # GTF is 1-based inclusive -> BED 0-based half-open

            if feature == "gene":
                gout.write("%s\t%d\t%s\t%s\t%s\t%s\n" % (chrom, bed_start, end, gene_id, gene_name, strand))
                n_genes += 1
            else:
                transcript_id = attrs.get("transcript_id", "NA")
                exon_number = attrs.get("exon_number", "NA")
                eout.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    chrom, bed_start, end, gene_id, gene_name, transcript_id, exon_number, strand
                ))
                n_exons += 1

    sys.stderr.write("gtf_to_genebed: wrote %d gene intervals -> %s\n" % (n_genes, genes_out_path))
    sys.stderr.write("gtf_to_genebed: wrote %d exon intervals -> %s\n" % (n_exons, exons_out_path))


if __name__ == "__main__":
    main()
