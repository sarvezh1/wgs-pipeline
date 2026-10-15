#!/usr/bin/env python3
"""Render repository-authored, genome-scale Plotly tracks from pipeline outputs.

The regional renderer remains VCF2Circos-derived.  This renderer only consumes
standardized/native derivatives already produced by the pipeline and does not
perform variant, SV, CNV, or coverage calling.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import Counter
from pathlib import Path


PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
COLORS = {
    "ideogram_a": "#475569", "ideogram_b": "#94A3B8",
    "snv": "#356FA8", "indel": "#D97706",
    "sample": "#0F766E", "tumor": "#B91C1C", "normal": "#2563EB",
    "gain": "#C2410C", "loss": "#1D4ED8", "neutral": "#9CA3AF",
    "DEL": "#DC2626", "DUP": "#7C3AED", "INV": "#EA580C", "BND": "#0891B2",
}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_fai(path: Path):
    contigs = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split("\t")
        if len(fields) < 2 or not fields[1].isdigit() or int(fields[1]) < 1:
            raise ValueError(f"malformed reference FAI at line {number}")
        contigs.append((fields[0], int(fields[1])))
    if not contigs:
        raise ValueError("reference FAI is empty")
    lengths = dict(contigs)
    if any(name in lengths for name in PRIMARY):
        included = [(name, lengths[name]) for name in PRIMARY if name in lengths]
    else:
        bare = [name.removeprefix("chr") for name in PRIMARY]
        included = [(name, lengths[name]) for name in bare if name in lengths]
    if not included:
        raise ValueError("reference FAI contains no canonical human chromosomes")
    included_names = {name for name, _ in included}
    excluded = [(name, length) for name, length in contigs if name not in included_names]
    return included, excluded


def info_map(text: str):
    result = {}
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        result[key] = value if separator else ""
    return result


def variant_kind(ref: str, alt: str):
    alleles = alt.split(",")
    return "SNV" if len(ref) == 1 and all(len(value) == 1 and value not in ".*" for value in alleles) else "INDEL"


def read_small_variants(path: Path, lengths: dict[str, int], bin_size: int):
    counts = Counter()
    stats = Counter()
    saw_header = False
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith("#CHROM\t"):
                saw_header = True
                continue
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 7 or not fields[1].isdigit():
                raise ValueError(f"malformed small-variant VCF record at line {number}")
            stats["input_records"] += 1
            chrom, pos, ref, alt, filt = fields[0], int(fields[1]), fields[3], fields[4], fields[6]
            if filt != "PASS":
                stats["non_pass_records"] += 1
                continue
            if chrom not in lengths:
                stats["excluded_contig_records"] += 1
                continue
            if pos < 1 or pos > lengths[chrom]:
                stats["invalid_coordinate_records"] += 1
                continue
            kind = variant_kind(ref, alt)
            counts[(chrom, (pos - 1) // bin_size, kind)] += 1
            stats[f"pass_{kind.lower()}"] += 1
    if not saw_header:
        raise ValueError("small-variant VCF is missing its #CHROM header")
    return counts, stats


def parse_labeled_path(text: str):
    label, separator, raw_path = text.partition("=")
    if not separator or not label or not raw_path:
        raise ValueError(f"expected LABEL=PATH, received {text!r}")
    return label, Path(raw_path)


def read_coverage(specs, lengths):
    records, stats = [], Counter()
    for spec in specs:
        label, path = parse_labeled_path(spec)
        with open_text(path) as handle:
            for number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 4:
                    raise ValueError(f"malformed mosdepth regions row at {path}:{number}")
                try:
                    start, end, depth = int(fields[1]), int(fields[2]), float(fields[-1])
                except ValueError as error:
                    raise ValueError(f"malformed mosdepth coordinates/depth at {path}:{number}") from error
                chrom = fields[0]
                if chrom not in lengths:
                    stats["excluded_contig_records"] += 1
                elif start < 0 or end <= start or end > lengths[chrom] or not math.isfinite(depth) or depth < 0:
                    stats["invalid_coordinate_records"] += 1
                else:
                    records.append((label, chrom, start, end, depth))
                    stats["records"] += 1
    return records, stats


def read_cnv(path: Path | None, cnv_format: str, lengths):
    if path is None:
        return [], Counter()
    records, stats, columns = [], Counter(), None
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("##") or raw.startswith("@"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if cnv_format == "germline_vcf":
                if raw.startswith("#CHROM\t"):
                    continue
                if len(fields) < 8 or not fields[1].isdigit():
                    raise ValueError(f"malformed germline CNV VCF at line {number}")
                info = info_map(fields[7])
                chrom, start, end = fields[0], int(fields[1]) - 1, int(info.get("END", fields[1]))
                kind = {"DUP": "gain", "DEL": "loss"}.get(info.get("SVTYPE", ""))
                if kind is None:
                    stats["undefined_class_records"] += 1
                    continue
            else:
                if "CONTIG" in fields and "START" in fields and "END" in fields and "CALL" in fields:
                    columns = {name: index for index, name in enumerate(fields)}
                    continue
                if columns is None:
                    raise ValueError("somatic called SEG is missing a CONTIG/START/END/CALL header")
                try:
                    chrom = fields[columns["CONTIG"]]
                    start = int(fields[columns["START"]]) - 1
                    end = int(fields[columns["END"]])
                    kind = {"+": "gain", "-": "loss", "0": "neutral"}.get(fields[columns["CALL"]].strip())
                except (IndexError, ValueError) as error:
                    raise ValueError(f"malformed somatic CNV SEG at line {number}") from error
                if kind is None:
                    stats["undefined_class_records"] += 1
                    continue
            if chrom not in lengths:
                stats["excluded_contig_records"] += 1
            elif start < 0 or end <= start or end > lengths[chrom]:
                stats["invalid_coordinate_records"] += 1
            else:
                records.append((chrom, start, end, kind))
                stats[kind] += 1
    return records, stats


def sv_kind(alt: str, info):
    kind = info.get("SVTYPE", "").upper()
    if not kind and ("[" in alt or "]" in alt):
        kind = "BND"
    return kind


def read_sv(path: Path | None, lengths):
    if path is None:
        return [], Counter()
    candidates, stats = [], Counter()
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[1].isdigit():
                raise ValueError(f"malformed Manta VCF at line {number}")
            chrom, pos, record_id, alt, filt = fields[0], int(fields[1]), fields[2], fields[4], fields[6]
            if filt != "PASS":
                stats["non_pass_records"] += 1
                continue
            info = info_map(fields[7])
            kind = sv_kind(alt, info)
            if kind not in {"BND", "DEL", "DUP", "INV"}:
                stats["unsupported_class_records"] += 1
                continue
            if kind == "BND":
                match = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
                if not match:
                    stats["unpaired_bnd_records"] += 1
                    continue
                other_chrom, other_pos = match.group(1), int(match.group(2))
            else:
                other_chrom, raw_end = chrom, info.get("END", "")
                if not raw_end.isdigit():
                    stats["missing_end_records"] += 1
                    continue
                other_pos = int(raw_end)
            if chrom not in lengths or other_chrom not in lengths:
                stats["excluded_contig_records"] += 1
                continue
            if pos < 1 or pos > lengths[chrom] or other_pos < 1 or other_pos > lengths[other_chrom]:
                stats["invalid_coordinate_records"] += 1
                continue
            mate = info.get("MATEID", "")
            key = tuple(sorted((record_id, mate))) if kind == "BND" and mate else (kind, chrom, pos, other_chrom, other_pos)
            candidates.append((key, kind, chrom, pos, other_chrom, other_pos))
    links, seen = [], set()
    for key, kind, chrom, pos, other_chrom, other_pos in candidates:
        if key in seen:
            continue
        seen.add(key)
        links.append((kind, chrom, pos, other_chrom, other_pos))
        stats[kind] += 1
    return links, stats


def genome_geometry(contigs):
    gap = sum(length for _, length in contigs) * 0.0025
    total = sum(length for _, length in contigs) + gap * len(contigs)
    geometry, cursor = {}, 0.0
    for name, length in contigs:
        start = 2 * math.pi * cursor / total
        end = 2 * math.pi * (cursor + length) / total
        geometry[name] = (start, end, length)
        cursor += length + gap
    return geometry


def angle(geometry, chrom, position):
    start, end, length = geometry[chrom]
    return start + (end - start) * position / length


def xy(theta, radius):
    return math.sin(theta) * radius, math.cos(theta) * radius


def arc_points(geometry, chrom, start, end, radius, count=24):
    a, b = angle(geometry, chrom, start), angle(geometry, chrom, end)
    return tuple(zip(*(xy(a + (b - a) * i / (count - 1), radius) for i in range(count))))


def rgba(hex_color, alpha):
    return f"rgba({int(hex_color[1:3], 16)},{int(hex_color[3:5], 16)},{int(hex_color[5:7], 16)},{alpha})"


def radial_band(geometry, chrom, start, end, inner_radius, outer_radius, count=4):
    outer_x, outer_y = arc_points(geometry, chrom, start, end, outer_radius, count)
    inner_x, inner_y = arc_points(geometry, chrom, end, start, inner_radius, count)
    return list(outer_x) + list(inner_x), list(outer_y) + list(inner_y)


def density_track(figure, contigs, bins, geometry, bin_size, kind, inner_radius, thickness, color):
    from plotly import graph_objects as go

    maximum = max([value for (chrom, _bin, record_kind), value in bins.items() if record_kind == kind] or [1])
    baseline_x, baseline_y = [], []
    polygon_x, polygon_y, polygon_text = [], [], []
    for chrom, length in contigs:
        base_x, base_y = arc_points(geometry, chrom, 0, length, inner_radius, 40)
        baseline_x.extend(list(base_x) + [None]); baseline_y.extend(list(base_y) + [None])
        for start in range(0, length, bin_size):
            end = min(start + bin_size, length)
            count = bins.get((chrom, start // bin_size, kind), 0)
            if count <= 0:
                continue
            outer_radius = inner_radius + thickness * math.sqrt(count / maximum)
            xs, ys = radial_band(geometry, chrom, start, end, inner_radius, outer_radius)
            polygon_x.extend(xs + [xs[0], None]); polygon_y.extend(ys + [ys[0], None])
            polygon_text.extend([f"<b>{kind} density</b><br>{chrom}:{start + 1:,}–{end:,}<br>{count:,} PASS alleles"] * (len(xs) + 1) + [None])
    figure.add_trace(go.Scatter(x=baseline_x, y=baseline_y, mode="lines", hoverinfo="skip", showlegend=False,
                                line={"color": rgba(color, .32), "width": 1}))
    figure.add_trace(go.Scatter(x=polygon_x, y=polygon_y, text=polygon_text, hoverinfo="text", mode="lines",
                                fill="toself", fillcolor=rgba(color, .82), name="SNV density" if kind == "SNV" else "INDEL density",
                                line={"color": color, "width": .35}, legendrank=1 if kind == "SNV" else 2))


def coverage_track(figure, contigs, rows, geometry, label, inner_radius, thickness, color):
    from plotly import graph_objects as go

    selected = [row for row in rows if row[0] == label]
    maximum = max([row[4] for row in selected] or [1]) or 1
    baseline_x, baseline_y = [], []
    polygon_x, polygon_y, polygon_text = [], [], []
    for chrom, length in contigs:
        base_x, base_y = arc_points(geometry, chrom, 0, length, inner_radius, 40)
        baseline_x.extend(list(base_x) + [None]); baseline_y.extend(list(base_y) + [None])
    for _, chrom, start, end, depth in selected:
        outer_radius = inner_radius + thickness * min(depth / maximum, 1)
        xs, ys = radial_band(geometry, chrom, start, end, inner_radius, outer_radius, 5)
        polygon_x.extend(xs + [xs[0], None]); polygon_y.extend(ys + [ys[0], None])
        polygon_text.extend([f"<b>{label} depth</b><br>{chrom}:{start + 1:,}–{end:,}<br>mean {depth:g}×"] * (len(xs) + 1) + [None])
    role = label.split(":", 1)[0].lower()
    display = "Tumor depth" if role == "tumor" else "Normal depth" if role == "normal" else "Depth"
    figure.add_trace(go.Scatter(x=baseline_x, y=baseline_y, mode="lines", hoverinfo="skip", showlegend=False,
                                line={"color": rgba(color, .30), "width": 1}))
    figure.add_trace(go.Scatter(x=polygon_x, y=polygon_y, text=polygon_text, hoverinfo="text", mode="lines",
                                fill="toself", fillcolor=rgba(color, .60), name=display, line={"color": color, "width": .65},
                                legendrank=3 if role == "tumor" else 4))


def build_figure(contigs, bins, coverage, cnv, links, bin_size, sample, caller, analysis, showcase=False):
    from plotly import graph_objects as go

    geometry = genome_geometry(contigs)
    figure = go.Figure()
    for index, (chrom, length) in enumerate(contigs):
        xs, ys = arc_points(geometry, chrom, 0, length, 1.0, 80)
        figure.add_trace(go.Scatter(x=xs, y=ys, mode="lines", hovertemplate=f"<b>{chrom}</b><br>{length:,} bp<extra></extra>",
                                    line={"color": COLORS[f'ideogram_{"a" if index % 2 == 0 else "b"}'], "width": 16}, showlegend=False))
        label_x, label_y = xy((geometry[chrom][0] + geometry[chrom][1]) / 2, 1.12)
        figure.add_annotation(x=label_x, y=label_y, text=chrom, showarrow=False, font={"size": 12, "color": "#1F2937"})

    density_track(figure, contigs, bins, geometry, bin_size, "SNV", .86, .085, COLORS["snv"])
    density_track(figure, contigs, bins, geometry, bin_size, "INDEL", .755, .085, COLORS["indel"])

    for cov_index, label in enumerate(dict.fromkeys(row[0] for row in coverage)):
        role = label.split(":", 1)[0].lower()
        coverage_track(figure, contigs, coverage, geometry, label, .65 - .105 * cov_index, .085,
                       COLORS.get(role, COLORS["sample"]))

    shown_cnv = set()
    for chrom, start, end, kind in cnv:
        xs, ys = arc_points(geometry, chrom, start, end, .43, 18)
        figure.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=f"CNV {kind}", legendgroup=f"cnv-{kind}", showlegend=kind not in shown_cnv,
                                    hovertemplate=f"<b>CNV {kind}</b><br>{chrom}:{start + 1:,}–{end:,}<extra></extra>",
                                    line={"width": 12, "color": COLORS[kind]}, legendrank={"gain": 5, "loss": 6, "neutral": 7}.get(kind, 7)))
        shown_cnv.add(kind)

    shown_sv = set()
    for kind, chrom, pos, other_chrom, other_pos in links:
        x1, y1 = xy(angle(geometry, chrom, pos), .34)
        x2, y2 = xy(angle(geometry, other_chrom, other_pos), .34)
        xs, ys = [], []
        for index in range(31):
            fraction = index / 30
            scale = (1 - fraction) ** 2
            reverse = fraction ** 2
            xs.append(scale * x1 + reverse * x2)
            ys.append(scale * y1 + reverse * y2)
        scope = "inter-chromosomal" if chrom != other_chrom else "intra-chromosomal"
        figure.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=kind, legendgroup=f"sv-{kind}", showlegend=kind not in shown_sv,
                                    hovertemplate=f"<b>{kind} · {scope}</b><br>{chrom}:{pos:,} ↔ {other_chrom}:{other_pos:,}<extra></extra>",
                                    opacity=.48, line={"width": 1.6, "color": COLORS[kind]}, legendrank={"DEL": 8, "DUP": 9, "INV": 10, "BND": 11}.get(kind, 11)))
        shown_sv.add(kind)

    if showcase:
        title_text = "<b>Genome-mode renderer demonstration</b><br><span style='font-size:16px;color:#475569'>Assembly: GRCh38 · Visualization mode: genome · Data status: visualization demonstration</span>"
        centre_text = "<b>Visualization-only demonstration</b><br><span style='font-size:11px'>Synthetic renderer test data · not a production result</span>"
    else:
        title_text = f"<b>{sample} {analysis} genome overview</b><br><span style='font-size:15px;color:#475569'>GRCh38 · {caller} · {bin_size:,} bp density bins</span>"
        centre_text = "<b>Genome-scale<br>variant context</b>"
    figure.update_layout(
        title={"text": title_text, "x": .5, "xanchor": "center"},
        width=1500, height=1100, paper_bgcolor="white", plot_bgcolor="white", font={"family": "Arial", "color": "#1F2937"},
        xaxis={"visible": False, "range": [-1.25, 1.25], "scaleanchor": "y"}, yaxis={"visible": False, "range": [-1.25, 1.25]},
        margin={"l": 50, "r": 300, "t": 125, "b": 55}, hovermode="closest",
        legend={"title": {"text": "<b>Genome tracks</b>"}, "font": {"size": 14}, "x": 1.01, "xanchor": "left", "y": .96, "bgcolor": "rgba(255,255,255,.95)", "bordercolor": "#CBD5E1", "borderwidth": 1},
        annotations=[*figure.layout.annotations, {"x": 0, "y": 0, "text": centre_text, "showarrow": False, "align": "center", "bgcolor": "rgba(255,255,255,.78)", "bordercolor": "rgba(203,213,225,.8)", "borderwidth": 1, "borderpad": 4, "font": {"size": 14, "color": "#475569"}}],
    )
    return figure


def write_outputs(args, contigs, excluded, bins, variant_stats, coverage, coverage_stats, cnv, cnv_stats, links, sv_stats):
    lengths = dict(contigs)
    with args.track_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["track", "series", "chrom", "start_0based", "end_0based_exclusive", "value", "unit"])
        for chrom, length in contigs:
            for start in range(0, length, args.bin_size):
                end = min(start + args.bin_size, length)
                index = start // args.bin_size
                writer.writerow(["small_variant_density", "PASS_SNV", chrom, start, end, bins[(chrom, index, "SNV")], "alleles_per_bin"])
                writer.writerow(["small_variant_density", "PASS_INDEL", chrom, start, end, bins[(chrom, index, "INDEL")], "alleles_per_bin"])
        for label, chrom, start, end, depth in coverage:
            writer.writerow(["coverage", label, chrom, start, end, f"{depth:g}", "mean_depth"])
        for chrom, start, end, kind in cnv:
            writer.writerow(["cnv", kind, chrom, start, end, 1, "segment"])
        for kind, chrom, pos, other_chrom, other_pos in links:
            writer.writerow(["sv_link", kind, chrom, pos - 1, pos, f"{other_chrom}:{other_pos}", "breakpoint_link"])

    rows = [
        ("renderer", "repository-authored Plotly genome renderer"), ("renderer_scope", "visualization only"),
        ("regional_implementation", "VCF2Circos 1.2.0 plus repository refinement; not used for this genome artifact"),
        ("sample_or_patient", args.sample), ("analysis_type", args.analysis_type), ("caller", args.caller),
        ("data_status", "visualization demonstration" if args.showcase else "pipeline-derived visualization"),
        ("reference_fai", args.fai.name), ("bin_size_bp", args.bin_size),
        ("included_contigs", ",".join(name for name, _ in contigs)),
        ("excluded_contigs", ",".join(name for name, _ in excluded) or "NONE"),
        ("small_variant_source", args.vcf.name), ("small_variant_selection", "FILTER=PASS"),
        ("pass_snv_alleles", variant_stats["pass_snv"]), ("pass_indel_or_other_alleles", variant_stats["pass_indel"]),
        ("small_variant_non_pass_omitted", variant_stats["non_pass_records"]),
        ("small_variant_excluded_contig_records", variant_stats["excluded_contig_records"]),
        ("small_variant_invalid_coordinate_records", variant_stats["invalid_coordinate_records"]),
        ("coverage_sources", ",".join(Path(value.partition("=")[2]).name for value in args.coverage) or "NOT_PROVIDED"),
        ("coverage_records", coverage_stats["records"]), ("cnv_source", args.cnv.name if args.cnv else "NOT_PROVIDED"),
        ("cnv_format", args.cnv_format if args.cnv else "NA"), ("cnv_segments", len(cnv)),
        ("sv_source", args.sv.name if args.sv else "NOT_PROVIDED"), ("sv_pass_links", len(links)),
        ("sv_link_classes", ",".join(f"{kind}:{sv_stats[kind]}" for kind in ("DEL", "DUP", "INV", "BND"))),
        ("invalid_optional_track_records", coverage_stats["invalid_coordinate_records"] + cnv_stats["invalid_coordinate_records"] + sv_stats["invalid_coordinate_records"]),
    ]
    with args.provenance.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--fai", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--analysis-type", choices=["germline", "somatic"], required=True)
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    parser.add_argument("--coverage", action="append", default=[], metavar="ROLE:SAMPLE=PATH")
    parser.add_argument("--cnv", type=Path)
    parser.add_argument("--cnv-format", choices=["germline_vcf", "somatic_seg"])
    parser.add_argument("--sv", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--track-summary", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--no-render", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--showcase", action="store_true", help="label output as a synthetic visualization demonstration")
    args = parser.parse_args(argv)
    if args.bin_size < 1:
        parser.error("--bin-size must be a positive integer")
    if bool(args.cnv) != bool(args.cnv_format):
        parser.error("--cnv and --cnv-format must be supplied together")
    return args


def main(argv=None):
    args = parse_args(argv)
    contigs, excluded = read_fai(args.fai)
    lengths = dict(contigs)
    bins, variant_stats = read_small_variants(args.vcf, lengths, args.bin_size)
    coverage, coverage_stats = read_coverage(args.coverage, lengths)
    cnv, cnv_stats = read_cnv(args.cnv, args.cnv_format, lengths)
    links, sv_stats = read_sv(args.sv, lengths)
    write_outputs(args, contigs, excluded, bins, variant_stats, coverage, coverage_stats, cnv, cnv_stats, links, sv_stats)
    if not args.no_render:
        figure = build_figure(contigs, bins, coverage, cnv, links, args.bin_size, args.sample, args.caller, args.analysis_type, args.showcase)
        figure.write_html(args.output, include_plotlyjs=True, full_html=True, config={"responsive": True, "displaylogo": False, "scrollZoom": True})
        html = args.output.read_text(encoding="utf-8")
        if any(value in html for value in ("/tmp/", "/work/", "file://", "<script src=")):
            raise ValueError("generated HTML contains a non-portable path or external Plotly reference")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"Genome visualization error: {error}", file=sys.stderr)
        raise SystemExit(2)
