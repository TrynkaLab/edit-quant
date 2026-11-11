#!/usr/bin/env python3
"""Build the self-contained edit-quant summary report HTML page, aggregated over all samples.

Layout follows the tglow-pipeline / tglow-core QC report (tglow.qc.render): one
Jinja2 template (assets/report.html.j2) with Plotly figures embedded as divs and
plotly.js injected once, so the report is a single file that works offline.

  Tab 1 Summary        - run-level counts, edited and WT samples split, read depth histogram
  Tab 2 Design QC      - amplicon extraction and guide orientation/PAM checks, i.e. whether
                         CRISPResso was handed a sane amplicon and guide in the first place
  Tab 3 QC             - per-sample success/success_wt/fail/not-tested counts and the full pair table
  Tab 4 Plates         - per design plate heatmaps of edit % for each sample
  Tab 5 Amplicon/donor - amplicon x donor heatmap of edit %

A "pair" is one sample x one (replicate-merged) amplicon in that sample's project panel:

  success     edited sample, reads_aligned >= min_reads and Modified% >= edit_threshold
  success_wt  WT sample, reads_aligned >= min_reads and Modified% <= wt_max_edit (no editing)
  fail        reads_aligned >= min_reads and either an edited sample with
              Modified% < edit_threshold, or a WT sample with Modified% > wt_max_edit
  not_tested  no reads assigned (no_reads), reads assigned but none aligned to the
              amplicon (not_aligned), fewer than min_reads aligned (low_depth), or
              the sample's CRISPResso run produced no output (crispresso_failed)

Samples whose donor (or sample name, if donor is 'none') matches wt_pattern are
counted separately, so they don't contribute to the edited-sample numbers.
Also writes the pair table as TSV (the machine-readable companion to the report).
"""

import argparse
import base64
import glob
import gzip
import json
import logging
import math
import os
import re
import string
from datetime import datetime

import jinja2
import pandas as pd
import plotly.graph_objects as go
import plotly.offline

from pipeline_common import read_rows

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

STATUS_ORDER = ["success", "success_wt", "fail", "not_tested"]
STATUS_COLORS = {"success": "#0e9e86", "success_wt": "#3d85c6", "fail": "#e07a5f", "not_tested": "#9e9e9e"}

# Design-QC statuses (build_report.py) that mean "usable as designed". Everything
# else is shown as flagged, so a status added upstream shows up rather than being
# silently counted as a pass.
AMPLICON_PASS = {"OK"}
GUIDE_PASS = {"OK"}

# rows x cols for each supported plate format (as tglow.qc.plate_layout)
PLATE_FORMATS = {24: (4, 6), 96: (8, 12), 384: (16, 24), 1536: (32, 48)}
HEATMAP_WIDTH = 520
HEATMAP_HEIGHT = 350


def fig_to_div(fig):
    """Render a Plotly figure to a <div> HTML fragment (plotly.js itself is injected once, see build_report)."""
    if fig is None:
        return None
    return fig.to_html(include_plotlyjs=False, full_html=False, config={"displayModeBar": False})


def slugify_sample(sample):
    """Mirror of lib/utils.nf slugifySample - the name CRISPResso writes its output under."""
    return re.sub(r"_{2,}", "_", re.sub(r"[\s'*\"/\\\[\]:;|,<>?]", "_", sample))


def fmt_pct(value):
    return "n/a" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.1f}%"


# ---------------------------------------------------------------------------------- #
#                                      Loading
# ---------------------------------------------------------------------------------- #
def load_samples(sample_check_path, wt_pattern):
    samples = pd.read_csv(sample_check_path, sep="\t", dtype=str, keep_default_na=False)
    wt_re = re.compile(wt_pattern)
    samples["label"] = [d if d != "none" else s for s, d in zip(samples["sample"], samples["donor"])]
    samples["is_wt"] = [bool(wt_re.search(label)) for label in samples["label"]]
    samples["slug"] = samples["sample"].map(slugify_sample)
    return samples


def load_amplicons(amplicon_manifest_path, replicates_path, amplicon_report_path):
    """One row per designed well: project, plate, well, gene, merged CRISPResso amplicon name (or NA)."""
    wells = pd.DataFrame([
        {"amplicon_name": sid, "project_name": row["project_name"], "plate_name": row["plate_name"],
         "well": row["well_position"], "gene_name": row["gene_name"]}
        for sid, row in read_rows(amplicon_manifest_path)
    ])
    replicates = pd.read_csv(replicates_path, sep="\t", dtype=str)
    extraction = pd.read_csv(amplicon_report_path, sep="\t", dtype=str)[["sample_id", "status"]]
    extraction = extraction.rename(columns={"sample_id": "amplicon_name", "status": "extraction_status"})

    wells = wells.merge(replicates[["amplicon_name", "crispresso_amplicon_name"]], on="amplicon_name", how="left")
    wells = wells.merge(extraction, on="amplicon_name", how="left")
    wells["extracted"] = wells["crispresso_amplicon_name"].notna()
    return wells


def load_design(amplicon_report_path, guide_qc_path, wells):
    """The two BUILD_REPORT design-QC tables, each annotated with the well they
    came from (project/plate/well), which the TSVs identify only by sample_id."""
    well_cols = ["amplicon_name", "project_name", "plate_name", "well"]
    keys = wells[well_cols].drop_duplicates("amplicon_name")

    amplicons = pd.read_csv(amplicon_report_path, sep="\t", dtype=str, keep_default_na=False)
    guides = pd.read_csv(guide_qc_path, sep="\t", dtype=str, keep_default_na=False)
    for table in (amplicons, guides):
        table["status"] = table["status"].replace("", "UNKNOWN")
    return (
        amplicons.merge(keys, left_on="sample_id", right_on="amplicon_name", how="left"),
        guides.merge(keys, left_on="sample_id", right_on="amplicon_name", how="left"),
    )


def load_quantification(crispresso_dir):
    """dict sample slug -> DataFrame of SAMPLES_QUANTIFICATION_SUMMARY.txt (Name, Modified%, Reads_aligned)."""
    quant = {}
    for path in glob.glob(os.path.join(crispresso_dir, "CRISPRessoPooled_on_*", "SAMPLES_QUANTIFICATION_SUMMARY.txt")):
        slug = os.path.basename(os.path.dirname(path))[len("CRISPRessoPooled_on_"):]
        df = pd.read_csv(path, sep="\t")
        quant[slug] = df.set_index("Name")[["Modified%", "Reads_total", "Reads_aligned"]]
    return quant


def load_viewer_data(viewer_dir):
    """Merge the per-sample PACK_VIEWER_DATA payloads into one viewer payload.

    Each sample was packed independently, so its amplicon offsets are relative to its
    own buffer; concatenating the buffers means rebasing them. Re-gzipping the whole
    thing once compresses better than keeping 30 separately-compressed blobs.
    """
    paths = sorted(glob.glob(os.path.join(viewer_dir, "*.viewer.json.gz")))
    if not paths:
        return None
    amplicons, blob, series = [], bytearray(), None
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            payload = json.load(fh)
        buf = gzip.decompress(base64.b64decode(payload["blob"]))
        base = len(blob)
        for rec in payload["amplicons"]:
            rec["off"] += base
            amplicons.append(rec)
        blob += buf
        series = series or payload.get("series")
    if not amplicons:
        return None
    log.info("amplicon viewer: %d amplicons from %d sample(s), %d KB packed",
             len(amplicons), len(paths), len(blob) / 1024)
    return {
        "meta": {"generated": datetime.now().strftime("%Y-%m-%d %H:%M"), "series": series},
        "blob": base64.b64encode(gzip.compress(bytes(blob), 9)).decode("ascii"),
        "amplicons": amplicons,
    }


def build_viewer_tab(viewer_dir, assets_dir):
    """Viewer tab context: the merged payload plus the shared css/markup/js assets.

    The same three assets drive bin/build_amplicon_viewer.py's standalone page, so the
    report tab and the standalone build cannot drift.
    """
    if not viewer_dir:
        return {"available": False}
    data = load_viewer_data(viewer_dir)
    if data is None:
        return {"available": False}

    def asset(name):
        path = os.path.join(assets_dir, name)
        if not os.path.exists(path):
            # Under Nextflow the template is staged into the task directory, so its
            # own directory is NOT the assets directory -- the viewer assets have to
            # be staged too and pointed at with --viewer_assets.
            raise SystemExit(
                f"viewer asset '{name}' not found in {assets_dir!r}. Pass --viewer_assets "
                f"with the directory holding amplicon_viewer.css/.body.html/.js, or run "
                f"with --viewer_dir omitted to build the report without the amplicon tab."
            )
        with open(path, encoding="utf-8") as fh:
            return fh.read()

    return {
        "available": True,
        "n_amplicons": len(data["amplicons"]),
        "n_samples": len({a["sample"] for a in data["amplicons"]}),
        "css": asset("amplicon_viewer.css"),
        "body": asset("amplicon_viewer.body.html"),
        "js": asset("amplicon_viewer.js"),
        # json.dumps leaves '<' alone, which would end the inline <script> early
        "json": json.dumps(data, separators=(",", ":")).replace("</", "<\\/"),
    }


def build_pairs(samples, wells, quant, edit_threshold, wt_max_edit, min_reads):
    """One row per (sample, merged amplicon in its project panel) with status/reason."""
    panel = (wells[wells["extracted"]]
             .groupby(["project_name", "crispresso_amplicon_name"], as_index=False)
             .agg(gene_name=("gene_name", "first"), plate_name=("plate_name", "first"),
                  wells=("well", lambda w: ",".join(w)), n_replicates=("well", "size")))

    rows = []
    for _, sample in samples[samples["status"] == "OK"].iterrows():
        sample_quant = quant.get(sample["slug"])
        for _, amp in panel[panel["project_name"] == sample["project_name"]].iterrows():
            reads_total = reads_aligned = modified = float("nan")
            if sample_quant is None:
                status, reason = "not_tested", "crispresso_failed"
            elif amp["crispresso_amplicon_name"] not in sample_quant.index:
                status, reason = "not_tested", "missing_from_output"
            else:
                q = sample_quant.loc[amp["crispresso_amplicon_name"]]
                reads_total = float(q["Reads_total"]) if pd.notna(q["Reads_total"]) else 0.0
                reads_aligned = float(q["Reads_aligned"]) if pd.notna(q["Reads_aligned"]) else 0.0
                modified = float(q["Modified%"]) if pd.notna(q["Modified%"]) else float("nan")
                if reads_total == 0:
                    status, reason = "not_tested", "no_reads"
                elif reads_aligned == 0:
                    status, reason = "not_tested", "not_aligned"
                elif reads_aligned < min_reads:
                    status, reason = "not_tested", "low_depth"
                elif sample["is_wt"]:
                    # WT should be unedited: no editing is a success in its own right
                    status, reason = ("success_wt", "") if modified <= wt_max_edit else ("fail", "")
                elif modified >= edit_threshold:
                    status, reason = "success", ""
                else:
                    status, reason = "fail", ""
            rows.append({
                "sample": sample["sample"], "project_name": sample["project_name"], "donor": sample["donor"],
                "label": sample["label"], "is_wt": sample["is_wt"],
                "amplicon": amp["crispresso_amplicon_name"], "gene_name": amp["gene_name"],
                "plate_name": amp["plate_name"], "wells": amp["wells"], "n_replicates": amp["n_replicates"],
                "reads_total": reads_total, "reads_aligned": reads_aligned, "modified_pct": modified,
                "status": status, "reason": reason,
            })
    return pd.DataFrame(rows, columns=[
        "sample", "project_name", "donor", "label", "is_wt", "amplicon", "gene_name", "plate_name", "wells",
        "n_replicates", "reads_total", "reads_aligned", "modified_pct", "status", "reason",
    ])


# ---------------------------------------------------------------------------------- #
#                                  Tab 1: Summary
# ---------------------------------------------------------------------------------- #
def group_stats(pairs, samples, wells, is_wt):
    """Counts for either the edited (is_wt=False) or WT (is_wt=True) samples."""
    group_samples = samples[samples["is_wt"] == is_wt]
    run = group_samples[group_samples["status"] == "OK"]
    group_pairs = pairs[pairs["is_wt"] == is_wt]
    tested = group_pairs[group_pairs["status"] != "not_tested"]
    n_success = int((group_pairs["status"] == "success").sum())
    n_success_wt = int((group_pairs["status"] == "success_wt").sum())
    n_fail = int((group_pairs["status"] == "fail").sum())

    return {
        "n_samples": len(group_samples),
        "n_run": len(run),
        "n_skipped": int((group_samples["status"] != "OK").sum()),
        "n_failed": int(group_pairs.loc[group_pairs["reason"] == "crispresso_failed", "sample"].nunique()),
        "n_pairs": len(group_pairs),
        "n_tested": len(tested),
        "n_not_tested": len(group_pairs) - len(tested),
        "n_success": n_success,
        "pct_success": fmt_pct(100 * n_success / len(tested) if len(tested) else None),
        "n_success_wt": n_success_wt,
        "pct_success_wt": fmt_pct(100 * n_success_wt / len(tested) if len(tested) else None),
        "n_fail": n_fail,
        "pct_fail": fmt_pct(100 * n_fail / len(tested) if len(tested) else None),
        "mean_edit": fmt_pct(tested["modified_pct"].mean() if len(tested) else None),
        "median_edit": fmt_pct(tested["modified_pct"].median() if len(tested) else None),
    }


def build_read_depth_histogram(pairs, min_reads):
    depth = pairs[pairs["reads_aligned"].notna()]
    fig = go.Figure()
    for is_wt, name in ((False, "Edited"), (True, "WT")):
        values = depth.loc[depth["is_wt"] == is_wt, "reads_aligned"]
        if len(values):
            # log10(reads + 1) so amplicons with no reads still show up at 0
            fig.add_trace(go.Histogram(x=(values + 1).map(math.log10), name=name, nbinsx=50, opacity=0.75))
    fig.add_vline(x=math.log10(min_reads + 1), line_dash="dash", annotation_text=f"min_reads = {min_reads}")
    fig.update_layout(
        title="Aligned reads per sample x amplicon pair", barmode="overlay",
        xaxis_title="log10(aligned reads + 1)", yaxis_title="Pairs",
    )
    return fig


def build_summary_tab(pairs, samples, wells, min_reads):
    run_projects = set(samples.loc[samples["status"] == "OK", "project_name"])
    run_wells = wells[wells["project_name"].isin(run_projects)]
    unique_amplicons = run_wells.loc[run_wells["extracted"], "crispresso_amplicon_name"].unique()

    # An amplicon "failed to amplify" if it has CRISPResso output in at least one
    # edited sample, but no edited sample has enough aligned reads for it
    edited_pairs = pairs[(~pairs["is_wt"]) & (~pairs["reason"].isin(["crispresso_failed", "missing_from_output"]))]
    edited_tested = edited_pairs[edited_pairs["status"] != "not_tested"]

    amplicon_stats = {
        "n_designed": len(run_wells),
        "n_unique": len(unique_amplicons),
        "n_not_extracted": int((~run_wells["extracted"]).sum()),
        "n_not_amplified": len(set(edited_pairs["amplicon"]) - set(edited_tested["amplicon"])),
    }

    return {
        "edited": group_stats(pairs, samples, wells, is_wt=False),
        "wt": group_stats(pairs, samples, wells, is_wt=True),
        "amplicons": amplicon_stats,
        "n_input_samples": len(samples),
        "read_depth_html": fig_to_div(build_read_depth_histogram(pairs, min_reads)),
    }


# ---------------------------------------------------------------------------------- #
#                                  Tab 2: Design QC
# ---------------------------------------------------------------------------------- #
def status_counts(table, pass_statuses):
    """[{status, n, pass}] ordered by count, passing statuses first."""
    counts = table["status"].value_counts()
    return sorted(
        ({"status": status, "n": int(n), "pass": status in pass_statuses} for status, n in counts.items()),
        key=lambda c: (not c["pass"], -c["n"], c["status"]),
    )


def build_design_tab(amplicons, guides, pam, pam_side, qw_center):
    """Did each design well give CRISPResso a usable amplicon, and a guide it can
    actually find in that amplicon and orient correctly? Both are upstream of any
    editing number: a guide CRISPResso can't place leaves the quantification
    window spanning the whole amplicon, and one placed the wrong way round puts
    it at the far end of the protospacer, so the edit % below is measured in the
    wrong place rather than being obviously absent."""
    n_amplicons_pass = int(amplicons["status"].isin(AMPLICON_PASS).sum())
    n_supplied = int((amplicons.get("amplicon_source", pd.Series(dtype=str)) == "manifest").sum())
    n_guides_pass = int(guides["status"].isin(GUIDE_PASS).sum())
    n_guides_passed_on = int((guides["passed_to_crispresso"] == "yes").sum())

    amplicon_cols = ["sample_id", "project_name", "plate_name", "well", "gene_name_declared",
                     "chrom", "amplicon_start", "amplicon_end", "amplicon_length",
                     "declared_amplicon_size", "length_diff", "orientation_ok",
                     "amplicon_source", "reference_check", "primers_in_supplied",
                     "reference_length", "status"]
    guide_cols = ["sample_id", "project_name", "plate_name", "well", "gene_name_declared",
                  "seq_type", "guide_seq", "strand_in_amplicon", "match_start", "n_fwd_matches",
                  "n_rc_matches", "pam_seq", "cut_index", "passed_to_crispresso", "status"]

    return {
        "pam": pam,
        "pam_side": "3'" if pam_side == "3prime" else "5'",
        "qw_center": qw_center,
        "amplicons": {
            "n_total": len(amplicons),
            "n_pass": n_amplicons_pass,
            "n_flagged": len(amplicons) - n_amplicons_pass,
            "n_supplied": n_supplied,
            "counts": status_counts(amplicons, AMPLICON_PASS),
            "pass_statuses": sorted(AMPLICON_PASS),
            "rows": amplicons.reindex(columns=amplicon_cols).fillna("").to_dict("records"),
        },
        "guides": {
            "n_total": len(guides),
            "n_pass": n_guides_pass,
            "n_flagged": len(guides) - n_guides_pass,
            "n_passed_on": n_guides_passed_on,
            "n_dropped": len(guides) - n_guides_passed_on,
            "counts": status_counts(guides, GUIDE_PASS),
            "pass_statuses": sorted(GUIDE_PASS),
            "rows": guides.reindex(columns=guide_cols).fillna("").to_dict("records"),
        },
    }


# ---------------------------------------------------------------------------------- #
#                                     Tab 3: QC
# ---------------------------------------------------------------------------------- #
def build_status_barplot(pairs, samples):
    order = samples.loc[samples["status"] == "OK", "sample"].tolist()
    counts = pairs.groupby(["sample", "status"]).size().unstack(fill_value=0).reindex(order, fill_value=0)
    fig = go.Figure()
    for status in STATUS_ORDER:
        fig.add_trace(go.Bar(
            x=counts.index, y=counts[status] if status in counts else [0] * len(counts),
            name=status, marker_color=STATUS_COLORS[status],
        ))
    fig.update_layout(title="Pairs per sample", barmode="stack", xaxis_title="Sample", yaxis_title="Amplicons")
    return fig


def build_qc_tab(pairs, samples):
    per_sample = []
    for _, sample in samples.iterrows():
        sample_pairs = pairs[pairs["sample"] == sample["sample"]]
        tested = sample_pairs[sample_pairs["status"] != "not_tested"]
        counts = sample_pairs["status"].value_counts()
        if sample["status"] != "OK":
            state = sample["status"]
        elif (sample_pairs["reason"] == "crispresso_failed").any():
            state = "CRISPRESSO_FAILED"
        else:
            state = "OK"
        per_sample.append({
            "sample": sample["sample"], "project_name": sample["project_name"], "donor": sample["donor"],
            "is_wt": sample["is_wt"], "state": state,
            "success": int(counts.get("success", 0)), "success_wt": int(counts.get("success_wt", 0)),
            "fail": int(counts.get("fail", 0)),
            "not_tested": int(counts.get("not_tested", 0)),
            "mean_edit": fmt_pct(tested["modified_pct"].mean() if len(tested) else None),
        })

    table = pairs.assign(
        reads_aligned=pairs["reads_aligned"].map(lambda v: "" if math.isnan(v) else f"{v:.0f}"),
        modified_pct=pairs["modified_pct"].map(lambda v: "" if math.isnan(v) else f"{v:.1f}"),
    ).to_dict("records")

    return {
        "status_barplot_html": fig_to_div(build_status_barplot(pairs, samples)),
        "per_sample": per_sample,
        "pairs": table,
        "statuses": STATUS_ORDER,
    }


# ---------------------------------------------------------------------------------- #
#                               Tab 4: Plate heatmaps
# ---------------------------------------------------------------------------------- #
def well_to_index(well):
    match = re.match(r"^([A-Za-z]+)0*([0-9]+)$", well.strip())
    if not match:
        return None
    row = 0
    for ch in match.group(1).upper():
        row = row * 26 + (string.ascii_uppercase.index(ch) + 1)
    return row - 1, int(match.group(2)) - 1


def infer_plate_format(plate_wells, override="auto"):
    if override not in (None, "auto"):
        return PLATE_FORMATS[int(override)]
    indices = [i for i in map(well_to_index, plate_wells) if i is not None]
    max_row = max(r for r, _ in indices) + 1
    max_col = max(c for _, c in indices) + 1
    for fmt in sorted(PLATE_FORMATS):
        n_rows, n_cols = PLATE_FORMATS[fmt]
        if max_row <= n_rows and max_col <= n_cols:
            return n_rows, n_cols
    return max_row, max_col


def style_heatmap_axes(fig):
    """As tglow.qc.plate_layout.style_heatmap_axes - boxed, every label, fixed size."""
    axis_style = dict(
        showgrid=False, zeroline=False, showline=True, linewidth=1,
        linecolor="rgba(136, 136, 136, 0.5)", mirror=True, ticks="",
        tickmode="linear", dtick=1, tickfont=dict(size=8),
    )
    fig.update_xaxes(**axis_style)
    fig.update_yaxes(**axis_style)
    fig.update_layout(autosize=False, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT)
    return fig


def build_plate_heatmap(plate_wells, sample_pairs, plate_format, title):
    """Edit % per well; designed wells that weren't tested are grey, empty wells blank."""
    n_rows, n_cols = plate_format
    row_labels = [string.ascii_uppercase[i % 26] if i < 26 else
                  f"{string.ascii_uppercase[i // 26 - 1]}{string.ascii_uppercase[i % 26]}" for i in range(n_rows)]
    col_labels = [str(i + 1) for i in range(n_cols)]

    values = [[None] * n_cols for _ in range(n_rows)]
    designed = [[None] * n_cols for _ in range(n_rows)]
    hover = [[""] * n_cols for _ in range(n_rows)]

    by_amplicon = sample_pairs.set_index("amplicon")
    for _, well in plate_wells.iterrows():
        index = well_to_index(well["well"])
        if index is None or index[0] >= n_rows or index[1] >= n_cols:
            continue
        r, c = index
        designed[r][c] = 1
        text = f"{well['well']} {well['gene_name']}"
        if not well["extracted"]:
            hover[r][c] = f"{text}<br>not extracted ({well['extraction_status']})"
            continue
        pair = by_amplicon.loc[well["crispresso_amplicon_name"]] if well["crispresso_amplicon_name"] in by_amplicon.index else None
        if pair is None:
            hover[r][c] = f"{text}<br>no result"
            continue
        reads = "n/a" if math.isnan(pair["reads_aligned"]) else f"{pair['reads_aligned']:.0f}"
        if pair["status"] == "not_tested":
            hover[r][c] = f"{text}<br>not tested ({pair['reason']})<br>reads: {reads}"
        else:
            values[r][c] = pair["modified_pct"]
            hover[r][c] = f"{text}<br>edit: {pair['modified_pct']:.1f}% ({pair['status']})<br>reads: {reads}"
        if pair["n_replicates"] > 1:
            hover[r][c] += f"<br>replicate wells: {pair['wells']}"

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=designed, x=col_labels, y=row_labels, colorscale=[[0, "#bdbdbd"], [1, "#bdbdbd"]],
        showscale=False, text=hover, hovertemplate="%{text}<extra></extra>",
    ))
    fig.add_trace(go.Heatmap(
        z=values, x=col_labels, y=row_labels, colorscale="Viridis", zmin=0, zmax=100,
        colorbar=dict(title="edit %"), text=hover, hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(title=title, xaxis_title="Column", yaxis_title="Row")
    return style_heatmap_axes(fig)


def build_plate_tab(pairs, samples, wells, plate_format_override):
    by_project = {}
    run = samples[samples["status"] == "OK"]
    for project, project_samples in run.groupby("project_name", sort=True):
        project_wells = wells[wells["project_name"] == project]
        plates = []
        for plate, plate_wells in project_wells.groupby("plate_name", sort=True):
            plate_format = infer_plate_format(plate_wells["well"], plate_format_override)
            figures = []
            # edited donors first, WT last
            for _, sample in project_samples.sort_values(["is_wt", "label"]).iterrows():
                sample_pairs = pairs[pairs["sample"] == sample["sample"]]
                title = f"{sample['label']} - plate {plate}" + (" (WT)" if sample["is_wt"] else "")
                figures.append(fig_to_div(build_plate_heatmap(plate_wells, sample_pairs, plate_format, title)))
            plates.append({"plate": plate, "figures": figures})
        by_project[project] = plates
    return {"available": bool(by_project), "by_project": by_project}


# ---------------------------------------------------------------------------------- #
#                           Tab 5: Amplicon x donor heatmap
# ---------------------------------------------------------------------------------- #
def build_donor_heatmap(project_pairs, project):
    tested = project_pairs.assign(
        value=project_pairs["modified_pct"].where(project_pairs["status"] != "not_tested")
    )
    columns = (tested[["label", "is_wt"]].drop_duplicates().sort_values(["is_wt", "label"])["label"].tolist())
    rows = tested.sort_values(["gene_name", "amplicon"])["amplicon"].unique().tolist()

    grid = tested.pivot_table(index="amplicon", columns="label", values="value", aggfunc="first", dropna=False)
    grid = grid.reindex(index=rows, columns=columns)

    info = tested.set_index(["amplicon", "label"])
    hover = []
    for amplicon in rows:
        line = []
        for label in columns:
            if (amplicon, label) not in info.index:
                line.append(f"{amplicon}<br>{label}: no result")
                continue
            pair = info.loc[(amplicon, label)]
            reads = "n/a" if math.isnan(pair["reads_aligned"]) else f"{pair['reads_aligned']:.0f}"
            if pair["status"] == "not_tested":
                line.append(f"{amplicon}<br>{label}: not tested ({pair['reason']})<br>reads: {reads}")
            else:
                line.append(f"{amplicon}<br>{label}: {pair['modified_pct']:.1f}% ({pair['status']})<br>reads: {reads}")
        hover.append(line)

    fig = go.Figure(data=go.Heatmap(
        z=grid.values, x=columns, y=rows, colorscale="Viridis", zmin=0, zmax=100,
        colorbar=dict(title="edit %"), text=hover, hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=8), tickmode="linear", dtick=1)
    fig.update_xaxes(side="top")
    fig.update_layout(
        title=f"Edit % per amplicon and donor - {project}", plot_bgcolor="#bdbdbd",
        height=max(400, 14 * len(rows) + 160), width=max(500, 80 * len(columns) + 380),
        margin=dict(l=320),
    )
    return fig


def build_donor_tab(pairs):
    figures = {project: fig_to_div(build_donor_heatmap(project_pairs, project))
               for project, project_pairs in pairs.groupby("project_name", sort=True)}
    return {"available": bool(figures), "by_project": figures}


# ---------------------------------------------------------------------------------- #
#                                      Report
# ---------------------------------------------------------------------------------- #
def build_report(args):
    samples = load_samples(args.sample_check, args.wt_pattern)
    wells = load_amplicons(args.amplicon_manifest, args.amplicon_replicates, args.amplicon_report)
    design_amplicons, design_guides = load_design(args.amplicon_report, args.guide_qc, wells)
    quant = load_quantification(args.crispresso_dir)
    pairs = build_pairs(samples, wells, quant, args.edit_threshold, args.wt_max_edit, args.min_reads)

    pairs.drop(columns=["label"]).to_csv(args.pairs_output, sep="\t", index=False, na_rep="NA")
    log.info(f"Wrote pair table to {args.pairs_output}")

    context = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "pipeline_version": args.pipeline_version,
        "params": {
            "edit_threshold": args.edit_threshold, "wt_max_edit": args.wt_max_edit,
            "min_reads": args.min_reads, "wt_pattern": args.wt_pattern,
        },
        "summary": build_summary_tab(pairs, samples, wells, args.min_reads),
        "design": build_design_tab(design_amplicons, design_guides, args.pam, args.pam_side, args.quantification_window_center),
        "qc": build_qc_tab(pairs, samples),
        "plates": build_plate_tab(pairs, samples, wells, args.plate_format),
        "donors": build_donor_tab(pairs),
        "viewer": build_viewer_tab(
            args.viewer_dir,
            args.viewer_assets or os.path.dirname(os.path.abspath(args.template)),
        ),
    }

    with open(args.template) as fh:
        template_source = fh.read()
    env = jinja2.Environment(autoescape=False, trim_blocks=True, lstrip_blocks=True)
    html = env.from_string(template_source).render(plotlyjs=plotly.offline.get_plotlyjs(), **context)

    with open(args.output, "w") as fh:
        fh.write(html)
    log.info(f"Wrote report to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the edit-quant summary report HTML page.")
    parser.add_argument("--crispresso_dir", required=True, help="Dir holding the CRISPRessoPooled_on_<sample> result dirs")
    parser.add_argument("--viewer_dir", help="Dir of per-sample *.viewer.json.gz from PACK_VIEWER_DATA; "
                                             "the amplicon tab is omitted when empty or absent")
    parser.add_argument("--viewer_assets", help="Dir holding amplicon_viewer.css/.body.html/.js "
                                                "(default: alongside --template, which only holds under "
                                                "a direct run -- Nextflow stages the template on its own)")
    parser.add_argument("--sample_check", required=True, help="reports/sample_check.tsv from CHECK_SAMPLES")
    parser.add_argument("--amplicon_replicates", required=True, help="reports/amplicon_replicates.tsv from CHECK_SAMPLES")
    parser.add_argument("--amplicon_report", required=True, help="reports/amplicon_report.tsv from BUILD_REPORT")
    parser.add_argument("--guide_qc", required=True, help="reports/guide_qc.tsv from BUILD_REPORT")
    parser.add_argument("--amplicon_manifest", required=True, help="Amplicon design manifest")
    parser.add_argument("--template", required=True, help="Jinja2 report template (assets/report.html.j2)")
    parser.add_argument("--edit_threshold", type=float, default=50, help="Modified%% at or above which a pair counts as a successful edit")
    parser.add_argument("--wt_max_edit", type=float, default=10, help="Modified%% at or below which a tested WT pair counts as success_wt (no editing)")
    parser.add_argument("--min_reads", type=int, default=100, help="Min aligned reads for a pair to count as tested")
    parser.add_argument("--wt_pattern", default="(^|_)WT$", help="Regex matched against donor (or sample if donor is 'none') to flag WT samples")
    parser.add_argument("--pam", default="NGG", help="PAM used for the guide orientation check (shown in the Design QC tab)")
    parser.add_argument("--pam_side", default="3prime", choices=["3prime", "5prime"],
                        help="Which end of the protospacer the PAM sits at (shown in the Design QC tab)")
    parser.add_argument("--quantification_window_center", type=int, default=-3,
                        help="CRISPResso -wc used for the reported cut index (shown in the Design QC tab)")
    parser.add_argument("--plate_format", default="auto", help="Plate format for heatmaps: auto, 24, 96, 384 or 1536")
    parser.add_argument("--pipeline_version", default="unknown", help="Shown in the report header")
    parser.add_argument("--output", required=True, help="Output HTML path")
    parser.add_argument("--pairs_output", required=True, help="Output TSV path for the pair table")
    build_report(parser.parse_args())
