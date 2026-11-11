#!/usr/bin/env python3
"""Pack one sample's CRISPRessoPooled output into a compact file for the report viewer.

Runs once per sample so the per-amplicon tables (tens of files per amplicon, hundreds of
amplicons per sample) are read in parallel and never staged into the single report task.
The output is one small gzipped JSON per sample, named after the sample so that
collecting them into RENDER_REPORT cannot collide -- Nextflow stages inputs by basename
and does not rename duplicates.

    pack_viewer_data.py --pooled CRISPRessoPooled_on_S1 --sample S1 --out S1.viewer.json.gz
"""

import argparse
import base64
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from viewer_common import SERIES, pack_pooled_dir, sample_name  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pooled", required=True, help="CRISPRessoPooled_on_<sample> directory")
    ap.add_argument("--sample", help="sample name (default: parsed from the directory name)")
    ap.add_argument("--out", required=True, help="output .viewer.json.gz")
    ap.add_argument("--min-reads", type=int, default=0, help="skip amplicons below this aligned-read count")
    ap.add_argument("--quantification-window-center", type=int, default=-3,
                    help="must match the run's setting; only affects guide strand labels")
    args = ap.parse_args()

    sample = args.sample or sample_name(args.pooled)
    if not os.path.isdir(args.pooled):
        # A sample whose CRISPResso task produced nothing still needs an output file,
        # or the collect() in the workflow would silently drop it.
        amplicons, blob, skipped = [], b"", 0
        print(f"{sample}: no CRISPRessoPooled directory, writing empty payload", file=sys.stderr)
    else:
        amplicons, blob, skipped = pack_pooled_dir(
            args.pooled, sample=sample,
            qw_center=args.quantification_window_center,
            min_reads=args.min_reads,
        )

    payload = {
        "sample": sample,
        "series": SERIES,
        "amplicons": amplicons,
        "blob": base64.b64encode(gzip.compress(blob, 9)).decode("ascii"),
    }
    with gzip.open(args.out, "wt", encoding="utf-8") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    print(
        f"{sample}: {len(amplicons)} amplicons, {skipped} skipped, "
        f"{len(blob) / 1024:.0f} KB packed -> {os.path.getsize(args.out) / 1024:.0f} KB file",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
