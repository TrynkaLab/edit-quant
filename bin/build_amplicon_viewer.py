#!/usr/bin/env python3
"""Build a standalone, self-contained amplicon viewer HTML from CRISPRessoPooled output.

Development/ad-hoc companion to the viewer tab in the pipeline report: the pipeline
route is PACK_VIEWER_DATA (one task per sample) feeding RENDER_REPORT, while this packs
any number of pooled directories into one page in a single pass, with no Nextflow.

    build_amplicon_viewer.py -o viewer.html CRISPRessoPooled_on_*/

The page is assembled from the same assets the report uses -- amplicon_viewer.css,
amplicon_viewer.body.html and amplicon_viewer.js -- so the two cannot drift.
"""

import argparse
import base64
import gzip
import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from viewer_common import SERIES, pack_pooled_dir, sample_name  # noqa: E402

ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets")


def read_asset(name):
    with open(os.path.join(ASSETS, name), encoding="utf-8") as fh:
        return fh.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pooled_dirs", nargs="+", help="CRISPRessoPooled_on_<sample> directories")
    ap.add_argument("-o", "--out", default="amplicon_viewer.html")
    ap.add_argument("--min-reads", type=int, default=0, help="skip amplicons below this aligned-read count")
    ap.add_argument("--quantification-window-center", type=int, default=-3,
                    help="must match the run's setting; only affects guide strand labels")
    ap.add_argument("--template", default=os.path.join(ASSETS, "amplicon_viewer.html"))
    args = ap.parse_args()

    amplicons, blob, skipped = [], bytearray(), 0
    for pooled in args.pooled_dirs:
        if not os.path.isdir(pooled):
            print(f"warning: not a directory, skipping: {pooled}", file=sys.stderr)
            continue
        sample = sample_name(pooled)
        recs, buf, n_skipped = pack_pooled_dir(
            pooled, sample=sample,
            qw_center=args.quantification_window_center,
            min_reads=args.min_reads,
        )
        base = len(blob)
        for rec in recs:
            rec["off"] += base          # rebase onto the combined buffer
        blob += buf
        amplicons += recs
        skipped += n_skipped
        print(f"{sample}: {len(recs)} amplicons", file=sys.stderr)

    if not amplicons:
        sys.exit("no amplicons with usable tables found")

    gz = gzip.compress(bytes(blob), 9)
    payload = {
        "meta": {
            "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "series": SERIES,
        },
        "blob": base64.b64encode(gz).decode("ascii"),
        "amplicons": amplicons,
    }

    html = read_asset(os.path.basename(args.template)) if os.path.dirname(args.template) == ASSETS \
        else open(args.template, encoding="utf-8").read()
    # json.dumps does not escape '<', so close the tag defensively for inline <script>
    blob_json = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    for marker, value in (
        ("/*__VIEWER_CSS__*/", read_asset("amplicon_viewer.css")),
        ("/*__VIEWER_BODY__*/", read_asset("amplicon_viewer.body.html")),
        ("/*__VIEWER_JS__*/", read_asset("amplicon_viewer.js")),
        ("/*__VIEWER_DATA__*/", blob_json),
    ):
        if marker not in html:
            sys.exit(f"template {args.template} has no {marker} placeholder")
        html = html.replace(marker, value)

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(
        f"\n{len(amplicons)} amplicons, {skipped} skipped\n"
        f"packed {len(blob) / 1024:.0f} KB -> gzip {len(gz) / 1024:.0f} KB\n"
        f"wrote {args.out} ({os.path.getsize(args.out) / 1e6:.1f} MB)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
