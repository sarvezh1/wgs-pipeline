#!/usr/bin/env python3
"""Render a portable, descriptive WGS analysis report from published TSV contracts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


NA = {"", "NA", "N/A", "na", "null", "None", "-"}
CALLER_VERSIONS = {
    "haplotypecaller": "GATK 4.6.2.0", "deepvariant": "1.10.0",
    "mutect2": "GATK 4.6.2.0", "strelka2": "2.9.10", "deepsomatic": "1.10.0",
}
SCHEMAS = {
    "variant_summary.tsv": {"caller", "analysis_type", "sample_or_patient", "total_records", "pass_records", "snvs", "indels"},
    "benchmark_summary.tsv": {"sample_or_patient", "analysis_type", "caller", "variant_class", "fp", "fn", "precision", "recall", "f1", "engine"},
    "benchmark_stratified_summary.tsv": {"caller", "stratum", "variant_class", "precision", "recall", "f1"},
    "concordance_summary.tsv": {"analysis_type", "variant_class"},
    "pairwise_jaccard.tsv": {"caller_a", "caller_b", "variant_class", "jaccard"},
    "genotype_concordance.tsv": {"variant_class", "genotype_concordant", "genotype_discordant", "missing_or_nocall", "genotype_concordance_rate"},
    "coverage_summary.tsv": {"sample_id", "patient_id", "role", "mean_depth", "region_bases", "engine", "engine_version"},
    "coverage_threshold_summary.tsv": {"sample_id", "threshold", "breadth_fraction"},
    "sv_summary.tsv": {"sample_or_patient", "total_sv", "pass_sv", "engine", "engine_version", "output_semantics"},
    "sv_class_summary.tsv": {"sv_class", "total", "pass"},
    "cnv_summary.tsv": {"sample_or_patient", "segments", "gain_segments", "loss_segments", "engine", "engine_version"},
    "performance_summary.tsv": {"process", "stage", "status", "cached", "realtime_seconds", "peak_rss_bytes"},
    "stage_summary.tsv": {"stage", "executed_tasks", "cached_tasks", "summed_task_wall_seconds", "summed_task_cpu_hours", "max_peak_rss_bytes"},
    "workflow_summary.tsv": {"workflow_elapsed_seconds", "executed_successful_tasks", "cached_tasks", "analysis_type", "caller_selection"},
    "caller_performance.tsv": {"caller", "summed_task_wall_seconds", "summed_task_cpu_hours", "max_peak_rss_bytes"},
}


class ReportError(ValueError):
    pass


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def available(value: object) -> bool:
    return str(value).strip() not in NA


def safe_display(value: object) -> str:
    """Remove machine-specific paths while retaining useful logical filenames."""
    text = str(value).strip()
    if not available(text):
        return "NA"
    text = re.sub(r"file://\S+", "[local file]", text)
    if text.startswith(("/", "~/")) or re.match(r"^[A-Za-z]:[\\/]", text):
        return Path(text).name or "[local path]"
    text = re.sub(r"(?<!\w)/(?:tmp|mnt|home|Users)/\S+", "[local path]", text)
    return text


def number(value: object, digits: int = 3) -> str:
    if not available(value):
        return "NA"
    try:
        parsed = float(str(value))
    except ValueError:
        return safe_display(value)
    if math.isnan(parsed) or math.isinf(parsed):
        return "NA"
    return f"{parsed:.{digits}f}"


def integer(value: object) -> str:
    if not available(value):
        return "NA"
    try:
        return f"{int(float(str(value))):,}"
    except ValueError:
        return safe_display(value)


def duration(value: object) -> str:
    if not available(value):
        return "NA"
    try:
        seconds = float(str(value))
    except ValueError:
        return safe_display(value)
    if seconds < 60:
        return f"{seconds:.1f} s"
    if seconds < 3600:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.2f} h"


def memory(value: object) -> str:
    if not available(value):
        return "NA"
    try:
        amount = float(str(value))
    except ValueError:
        return safe_display(value)
    return f"{amount / 2**30:.2f} GiB" if amount >= 2**30 else f"{amount / 2**20:.1f} MiB"


def bool_status(enabled: bool) -> str:
    return "Enabled" if enabled else "Not run"


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if not fields or any(not field for field in fields) or len(fields) != len(set(fields)):
                raise ReportError(f"Malformed summary input {path.name}: missing or duplicate header")
            missing = (required or set()) - set(fields)
            if missing:
                raise ReportError(f"Malformed summary input {path.name}: missing column(s): {', '.join(sorted(missing))}")
            rows = list(reader)
            if any(None in row for row in rows):
                raise ReportError(f"Malformed summary input {path.name}: inconsistent column count")
            return rows
    except UnicodeDecodeError as error:
        raise ReportError(f"Malformed summary input {path.name}: not UTF-8 text") from error


def discover(root: Path, suffix: str) -> list[Path]:
    paths = sorted(path for path in root.rglob(f"*{suffix}") if "report/figures" not in path.as_posix())
    if suffix in {"concordance_summary.tsv", "genotype_concordance.tsv"}:
        paths = [path for path in paths if ".stratified_" not in path.name]
    return paths


def load_contract(root: Path, suffix: str) -> tuple[list[dict[str, str]], list[Path]]:
    paths = discover(root, suffix)
    rows: list[dict[str, str]] = []
    required = SCHEMAS.get(suffix)
    for path in paths:
        rows.extend(read_tsv(path, required))
    return rows, paths


def local_href(report_dir: Path, target: Path) -> str:
    try:
        rel = target.resolve().relative_to(report_dir.parent.resolve()).as_posix()
        rel = "../" + rel
    except ValueError:
        return ""
    if rel.startswith("/") or "file:" in rel or "\\" in rel:
        return ""
    return rel


def configured_link(raw: str, label: str) -> str:
    if not raw:
        return ""
    posix = PurePosixPath(raw)
    parts = posix.parts
    allowed_parent = len(parts) >= 3 and parts[0] == ".." and parts[1] in {"multiqc", "vcf2circos"} and ".." not in parts[1:]
    if posix.is_absolute() or "\\" in raw or raw.startswith("file:") or (".." in parts and not allowed_parent):
        raise ReportError(f"Unsafe local link for {label}: use a relative path below the output directory")
    return f'<a class="button" href="{esc(raw)}">Open {esc(label)}</a>'


def table(headers: list[tuple[str, str]], rows: list[dict[str, object]], empty: str = "Not available") -> str:
    if not rows:
        return f'<p class="empty">{esc(empty)}</p>'
    head = "".join(f'<th scope="col">{esc(label)}</th>' for _, label in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{esc(safe_display(row.get(key, 'NA')))}</td>" for key, _ in headers) + "</tr>")
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def metric_card(label: str, value: str, note: str = "") -> str:
    return f'<div class="metric"><span>{esc(label)}</span><strong>{esc(value)}</strong>{f"<small>{esc(note)}</small>" if note else ""}</div>'


def svg_bars(title: str, labels: list[str], values: list[float], color: str, unit: str = "") -> str:
    if not labels or not values or max(values, default=0) <= 0:
        return ""
    width, left, right, row_h = 820, 180, 70, 34
    height = 54 + row_h * len(labels)
    span = width - left - right
    maximum = max(values)
    items = [f'<svg role="img" aria-label="{esc(title)}" viewBox="0 0 {width} {height}">', f'<title>{esc(title)}</title>']
    for index, (label, value) in enumerate(zip(labels, values)):
        y = 34 + index * row_h
        bar = max(1, span * value / maximum)
        items.extend([
            f'<text x="{left - 10}" y="{y + 15}" text-anchor="end">{esc(label)}</text>',
            f'<rect x="{left}" y="{y}" width="{bar:.2f}" height="20" rx="3" fill="{color}"/>',
            f'<text x="{min(left + bar + 8, width - 4)}" y="{y + 15}">{esc(number(value, 3))}{esc(unit)}</text>',
        ])
    items.append("</svg>")
    return "".join(items)


def svg_heatmap(rows: list[dict[str, str]]) -> str:
    usable = [row for row in rows if available(row.get("f1", ""))]
    callers = sorted({row["caller"] for row in usable})
    strata = sorted({row["stratum"] for row in usable})
    if not callers or not strata:
        return ""
    cell_w, cell_h, left, top = 130, 40, 190, 62
    width, height = left + cell_w * len(callers) + 20, top + cell_h * len(strata) + 25
    lookup = {(row["stratum"], row["caller"]): row["f1"] for row in usable}
    out = [f'<svg role="img" aria-label="Stratified F1 heatmap" viewBox="0 0 {width} {height}"><title>Stratified F1 heatmap</title>']
    for col, caller in enumerate(callers):
        out.append(f'<text x="{left + col * cell_w + cell_w / 2}" y="28" text-anchor="middle">{esc(caller)}</text>')
    for row_i, stratum in enumerate(strata):
        y = top + row_i * cell_h
        out.append(f'<text x="{left - 10}" y="{y + 25}" text-anchor="end">{esc(stratum)}</text>')
        for col, caller in enumerate(callers):
            raw = lookup.get((stratum, caller), "NA")
            value = float(raw) if available(raw) else None
            shade = "#edf0f4" if value is None else f"hsl(266 42% {92 - 43 * max(0, min(1, value)):.1f}%)"
            x = left + col * cell_w
            out.append(f'<rect x="{x}" y="{y}" width="{cell_w - 4}" height="{cell_h - 4}" rx="3" fill="{shade}"/>')
            out.append(f'<text x="{x + (cell_w - 4) / 2}" y="{y + 24}" text-anchor="middle" class="heat-value">{esc(number(raw, 3))}</text>')
    out.append("</svg>")
    return "".join(out)


def write_figure(directory: Path, name: str, svg: str) -> None:
    if svg:
        directory.mkdir(parents=True, exist_ok=True)
        standalone = svg.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" ', 1)
        (directory / name).write_text(standalone + "\n", encoding="utf-8")


def fastp_rows(root: Path) -> list[dict[str, str]]:
    rows = []
    for path in sorted(root.rglob("*_fastp.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            before, after = data["summary"]["before_filtering"], data["summary"]["after_filtering"]
            rows.append({
                "sample": path.name.removesuffix("_fastp.json"),
                "reads": integer(before.get("total_reads")),
                "retained": number(100 * after.get("total_reads", 0) / before.get("total_reads", 1), 2) + "%",
                "q30": number(100 * after.get("q30_rate", 0), 2) + "%",
            })
        except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as error:
            raise ReportError(f"Malformed fastp summary input {path.name}: {error}") from error
    return rows


def fastqc_status(root: Path) -> list[dict[str, str]]:
    rows = []
    for path in sorted(root.rglob("*_fastqc.zip")):
        try:
            with zipfile.ZipFile(path) as archive:
                member = next(name for name in archive.namelist() if name.endswith("/summary.txt"))
                statuses = [line.split("\t", 1)[0] for line in archive.read(member).decode("utf-8").splitlines() if line]
            overall = "FAIL" if "FAIL" in statuses else "WARN" if "WARN" in statuses else "PASS"
            rows.append({"read": path.name.removesuffix("_fastqc.zip"), "status": overall})
        except (OSError, zipfile.BadZipFile, StopIteration, UnicodeDecodeError) as error:
            raise ReportError(f"Malformed FastQC archive {path.name}: {error}") from error
    return rows


def alignment_qc(root: Path) -> list[dict[str, str]]:
    paths = sorted(root.rglob("multiqc_general_stats.txt"))
    if not paths:
        return []
    rows = read_tsv(paths[0])
    output = []
    for row in rows:
        sample = row.get("Sample", "")
        if sample not in {"tumor", "normal", "germline"}:
            continue
        output.append({
            "sample": sample,
            "mapped": f"{number(row.get('samtools_stats-reads_mapped_percent'), 3)}%",
            "paired": f"{number(row.get('samtools_stats-reads_properly_paired_percent'), 3)}%",
            "duplicates": f"{number(row.get('samtools_markdup-duplicate_fraction'), 3)}%",
            "error": number(row.get("samtools_stats-error_rate"), 5),
        })
    return output


def all_rows(root: Path, suffixes: list[str]) -> tuple[dict[str, list[dict[str, str]]], list[Path]]:
    contracts, paths = {}, []
    for suffix in suffixes:
        contracts[suffix], found = load_contract(root, suffix)
        paths.extend(found)
    return contracts, paths


def render(args: argparse.Namespace) -> str:
    root, output = args.input_root.resolve(), args.output.resolve()
    if not root.is_dir():
        raise ReportError(f"Input root does not exist: {args.input_root}")
    suffixes = list(SCHEMAS) + [
        "concordance_summary.tsv", "stratified_concordance_summary.tsv", "stratified_jaccard.tsv",
        "comparison_provenance.tsv", "benchmark_provenance.tsv", "benchmark_stratification_provenance.tsv",
        "performance_provenance.tsv", "resource_manifest.tsv", "gatk_path_manifest.tsv", "wgs_modality_provenance.tsv",
    ]
    contracts, source_paths = all_rows(root, list(dict.fromkeys(suffixes)))
    variants = contracts["variant_summary.tsv"]
    benchmarks = contracts["benchmark_summary.tsv"]
    stratified = contracts["benchmark_stratified_summary.tsv"]
    concordance = contracts["concordance_summary.tsv"]
    jaccard = contracts["pairwise_jaccard.tsv"]
    genotypes = contracts["genotype_concordance.tsv"]
    stratified_concordance = contracts["stratified_concordance_summary.tsv"]
    stratified_jaccard = contracts["stratified_jaccard.tsv"]
    coverage = contracts["coverage_summary.tsv"]
    breadth = contracts["coverage_threshold_summary.tsv"]
    sv = contracts["sv_summary.tsv"]
    sv_classes = contracts["sv_class_summary.tsv"]
    cnv = contracts["cnv_summary.tsv"]
    stages = contracts["stage_summary.tsv"]
    workflows = contracts["workflow_summary.tsv"]
    caller_perf = contracts["caller_performance.tsv"]
    perf_prov = contracts["performance_provenance.tsv"]
    benchmark_provenance = contracts.get("benchmark_provenance.tsv", [])

    samples = []
    if args.samplesheet:
        sample_rows = read_tsv(args.samplesheet.resolve(), {"sample_id", "patient_id", "role"})
        samples = [safe_display(row["sample_id"]) for row in sample_rows]
        patients = [safe_display(row["patient_id"]) for row in sample_rows]
        report_id = args.report_id or (patients[0] if patients else samples[0])
        roles = {row["role"]: safe_display(row["sample_id"]) for row in sample_rows}
    else:
        report_id = args.report_id or next((safe_display(row.get("sample_or_patient")) for row in variants), "Analysis")
        roles = {}
        samples = sorted({safe_display(row.get("sample_id")) for row in coverage if available(row.get("sample_id"))})
    callers = sorted({row.get("caller", "") for row in variants if row.get("caller")})
    if not callers:
        callers = [item.strip() for item in args.callers.split(",") if item.strip()]

    workflow = workflows[0] if workflows else {}
    provenance = perf_prov[0] if perf_prov else {}
    analysis = args.analysis_type
    fixture = args.resource_mode == "fixture"
    benchmark_enabled = args.run_benchmark
    comparison_enabled = args.caller == "compare"
    title = args.title or f"{report_id} {analysis.title()} WGS Analysis Report"
    subtitle = f"{args.assembly} | {'Tumor-normal' if analysis == 'somatic' else 'Germline'} | {'Technical fixture demonstration' if fixture else 'Integrated analysis summary'}"

    overall_bench = [row for row in benchmarks if row.get("variant_class", "").upper() in {"OVERALL", "ALL"}]
    mean_cov = next((row.get("mean_depth") for row in coverage if row.get("role") in {"tumor", "germline"}), coverage[0].get("mean_depth") if coverage else "NA")
    small_total = sum(int(row.get("total_records", "0")) for row in variants if row.get("total_records", "").isdigit()) if variants else None
    best_label_free_f1 = max((float(row["f1"]) for row in overall_bench if available(row.get("f1"))), default=None)
    shared = next((row.get("count") or row.get("allele_count") for row in concordance if (row.get("comparison_group") or row.get("agreement_category")) in {"all_callers", "all_three_shared", "shared"} and available(row.get("count") or row.get("allele_count"))), "NA")
    total_sv = sv[0].get("total_sv", "NA") if sv else ("Not run" if not args.run_sv else "Not available")
    cnv_segments = cnv[0].get("segments", "NA") if cnv else ("Not run" if not args.run_cnv else "Not available")
    elapsed = workflow.get("workflow_elapsed_seconds", "NA")
    peak = max((float(row["max_peak_rss_bytes"]) for row in stages if available(row.get("max_peak_rss_bytes"))), default=None)
    cards = [
        metric_card("Samples analyzed", str(len(samples)) if samples else "Not available", ", ".join(samples)),
        metric_card("Mean coverage", f"{number(mean_cov, 2)}×" if available(mean_cov) else ("Not run" if not args.run_coverage_qc else "Not available")),
        metric_card("Normalized small variants", integer(small_total) if small_total is not None else "Not available", "sum across caller outputs"),
        metric_card("Benchmark F1", number(best_label_free_f1, 3) if best_label_free_f1 is not None else ("Not run" if not benchmark_enabled else "Not available"), "highest displayed value; no ranking implied" if best_label_free_f1 is not None and len(overall_bench) > 1 else "overall"),
        metric_card("Shared caller calls", integer(shared) if available(shared) else ("Not run" if not comparison_enabled else "Not available")),
        metric_card("Final SVs", integer(total_sv) if available(total_sv) else total_sv),
        metric_card("CNV segments", integer(cnv_segments) if available(cnv_segments) else cnv_segments),
        metric_card("Workflow elapsed", duration(elapsed) if available(elapsed) else "Not available"),
        metric_card("Peak memory", memory(peak) if peak is not None else "Not available"),
    ]

    overview_rows = [
        {"field": "Analysis", "value": analysis.title()}, {"field": "Report identifier", "value": report_id},
        {"field": "Reference assembly", "value": args.assembly}, {"field": "Resource mode", "value": args.resource_mode},
        {"field": "Selected caller(s)", "value": ", ".join(callers) or args.caller},
        {"field": "BQSR", "value": bool_status(args.run_bqsr)}, {"field": "Benchmark", "value": bool_status(benchmark_enabled)},
        {"field": "Comparison", "value": bool_status(comparison_enabled)}, {"field": "VEP", "value": bool_status(args.run_vep)},
        {"field": "SV", "value": bool_status(args.run_sv)}, {"field": "CNV", "value": bool_status(args.run_cnv)},
        {"field": "Coverage QC", "value": bool_status(args.run_coverage_qc)},
        {"field": "Pipeline / Nextflow", "value": f"{args.pipeline_version} / {args.nextflow_version}"},
    ]
    if analysis == "somatic":
        overview_rows.insert(2, {"field": "Tumor / normal", "value": f"{roles.get('tumor', 'Not available')} / {roles.get('normal', 'Not available')}"})

    coverage_rows = [{"sample": row.get("sample_id", "NA"), "role": row.get("role", "NA"), "mean": f"{number(row.get('mean_depth'), 2)}×", "median": row.get("median_depth", "NA"), "bases": integer(row.get("region_bases"))} for row in coverage]
    breadth_rows = [{"sample": row.get("sample_id", "NA"), "threshold": f"≥{row.get('threshold', 'NA')}×", "breadth": number(row.get("breadth_fraction"), 4)} for row in breadth]
    coverage_svg = svg_bars("Mean coverage by sample", [f"{r.get('sample_id')} ({r.get('role')})" for r in coverage], [float(r["mean_depth"]) for r in coverage if available(r.get("mean_depth"))], "#3478b8", "×") if coverage and all(available(r.get("mean_depth")) for r in coverage) else ""

    variant_rows = [{
        "caller": row.get("caller", "NA"), "version": CALLER_VERSIONS.get(row.get("caller", ""), "Not available"),
        "records": integer(row.get("total_records")), "snv": integer(row.get("snvs")), "indel": integer(row.get("indels")),
        "pass": integer(row.get("pass_records")), "samples": row.get("sample_or_patient", "NA"),
        "benchmark": "Available" if any(b.get("caller") == row.get("caller") for b in benchmarks) else ("Not run" if not benchmark_enabled else "Not available"),
    } for row in variants]
    bench_rows = [{"caller": r.get("caller", "NA"), "class": r.get("variant_class", "NA"), "tp": integer(r.get("query_tp", r.get("truth_tp"))), "fp": integer(r.get("fp")), "fn": integer(r.get("fn")), "precision": number(r.get("precision"), 4), "recall": number(r.get("recall"), 4), "f1": number(r.get("f1"), 4)} for r in benchmarks]
    bench_svg = svg_bars("Overall benchmark F1 by caller", [r.get("caller", "NA") for r in overall_bench], [float(r["f1"]) for r in overall_bench if available(r.get("f1"))], "#7654a8") if overall_bench and all(available(r.get("f1")) for r in overall_bench) else ""

    jac_rows = [{"pair": f"{r.get('caller_a', 'NA')} / {r.get('caller_b', 'NA')}", "class": r.get("variant_class", "NA"), "jaccard": number(r.get("jaccard"), 4), "intersection": integer(r.get("intersection_count", r.get("intersection"))), "union": integer(r.get("union_count", r.get("union")))} for r in jaccard]
    genotype_rows = [{"class": r.get("variant_class", "NA"), "concordant": integer(r.get("genotype_concordant")), "discordant": integer(r.get("genotype_discordant")), "missing": integer(r.get("missing_or_nocall")), "rate": number(r.get("genotype_concordance_rate"), 4)} for r in genotypes]
    concord_rows = []
    for r in concordance:
        category = r.get("comparison_group") or r.get("agreement_category") or r.get("category") or r.get("set") or "summary"
        count = r.get("count") or r.get("allele_count") or r.get("variant_count") or "NA"
        concord_rows.append({"class": r.get("variant_class", "NA"), "category": category, "count": integer(count), "truth": number(r.get("truth_match_fraction"), 4)})
    overlap_source = [r for r in concordance if r.get("variant_class", "").upper() in {"ALL", "OVERALL"} and r.get("comparison_group") != "union" and available(r.get("count"))]
    overlap_svg = svg_bars("Caller overlap categories", [r.get("comparison_group", "category") for r in overlap_source], [float(r["count"]) for r in overlap_source], "#267f80") if any(float(r["count"]) > 0 for r in overlap_source) else ""

    strat_rows = [{"caller": r.get("caller", "NA"), "stratum": r.get("stratum", "NA"), "class": r.get("variant_class", "NA"), "precision": number(r.get("precision"), 4), "recall": number(r.get("recall"), 4), "f1": number(r.get("f1"), 4)} for r in stratified]
    heat_rows = [r for r in stratified if r.get("variant_class", "").upper() in {"ALL", "OVERALL"}]
    heat_svg = svg_heatmap(heat_rows)
    strat_concord_rows = [{"stratum": r.get("stratum", "NA"), "class": r.get("variant_class", "NA"), "category": r.get("comparison_group", "NA"), "count": integer(r.get("count")), "truth": number(r.get("truth_match_fraction"), 4)} for r in stratified_concordance]
    strat_jac_rows = [{"stratum": r.get("stratum", "NA"), "pair": f"{r.get('caller_a', 'NA')} / {r.get('caller_b', 'NA')}", "class": r.get("variant_class", "NA"), "jaccard": number(r.get("jaccard"), 4)} for r in stratified_jaccard]

    sv_rows = [{"caller": f"{r.get('engine', 'Manta')} {r.get('engine_version', '')}".strip(), "semantics": r.get("output_semantics", "final"), "total": integer(r.get("total_sv")), "pass": integer(r.get("pass_sv"))} for r in sv]
    sv_class_rows = [{"class": r.get("sv_class", "NA"), "total": integer(r.get("total")), "pass": integer(r.get("pass"))} for r in sv_classes]
    sv_svg = svg_bars("Final structural variants by class", [r.get("sv_class", "NA") for r in sv_classes], [float(r.get("total", 0)) for r in sv_classes], "#b06b32") if any(float(r.get("total", 0)) > 0 for r in sv_classes) else ""
    cnv_rows = [{"engine": f"{r.get('engine', 'NA')} {r.get('engine_version', '')}".strip(), "segments": integer(r.get("segments")), "gains": integer(r.get("gain_segments")), "losses": integer(r.get("loss_segments")), "neutral": integer(r.get("neutral_segments")), "altered": integer(r.get("total_altered_bases")), "largest": integer(r.get("largest_event_bp"))} for r in cnv]

    stage_rows = [{"stage": r.get("stage", "NA"), "fresh": integer(r.get("executed_tasks")), "cached": integer(r.get("cached_tasks")), "runtime": duration(r.get("summed_task_wall_seconds")), "cpu": number(r.get("summed_task_cpu_hours"), 3), "memory": memory(r.get("max_peak_rss_bytes"))} for r in stages]
    runtime_svg = svg_bars("Summed fresh-task runtime by stage", [r.get("stage", "NA") for r in stages if available(r.get("summed_task_wall_seconds"))], [float(r["summed_task_wall_seconds"]) for r in stages if available(r.get("summed_task_wall_seconds"))], "#5367a5", " s")
    compute_rows = []
    bench_by_caller = {r.get("caller"): r for r in overall_bench}
    for row in caller_perf:
        metric = bench_by_caller.get(row.get("caller"), {})
        compute_rows.append({"caller": row.get("caller", "NA"), "precision": number(metric.get("precision"), 4), "recall": number(metric.get("recall"), 4), "f1": number(metric.get("f1"), 4), "runtime": duration(row.get("summed_task_wall_seconds")), "cpu": number(row.get("summed_task_cpu_hours"), 3), "memory": memory(row.get("max_peak_rss_bytes"))})

    provenance_rows = [
        {"field": "Pipeline version", "value": args.pipeline_version}, {"field": "Nextflow version", "value": args.nextflow_version},
        {"field": "Analysis / resources", "value": f"{analysis} / {args.resource_mode}"}, {"field": "Reference assembly", "value": args.assembly},
        {"field": "Callers", "value": "; ".join(f"{c} {CALLER_VERSIONS.get(c, '')}".strip() for c in callers) or "Not available"},
        {"field": "Benchmark", "value": f"{args.benchmark_name} {args.benchmark_version}".strip() if benchmark_enabled else "Not run"},
        {"field": "Benchmark engine", "value": safe_display(overall_bench[0].get("engine")) if overall_bench else ("Not run" if not benchmark_enabled else "Not available")},
        {"field": "Benchmark ploidy handling", "value": (
            "Male-aware (hap.py)" if any(r.get("sex_aware_benchmarking", "").lower() == "true" and r.get("benchmark_sex_resolved", "").lower() == "male" for r in benchmark_provenance)
            else "Female-aware (hap.py)" if any(r.get("sex_aware_benchmarking", "").lower() == "true" and r.get("benchmark_sex_resolved", "").lower() == "female" for r in benchmark_provenance)
            else "Sex-aware handling not enabled" if benchmark_provenance
            else ("Not run" if not benchmark_enabled else "Not available")
        )},
        {"field": "VEP", "value": "Ensembl VEP 116" if args.run_vep else "Not run"},
        {"field": "SV / CNV", "value": f"{'Manta 1.6.0' if args.run_sv else 'Not run'} / {'GATK CNV 4.6.2.0' if args.run_cnv else 'Not run'}"},
        {"field": "Containers", "value": safe_display(provenance.get("container_identities", "Not available"))},
        {"field": "Metric scope", "value": safe_display(provenance.get("metric_scope", "TECHNICAL_VALIDATION_METRICS" if fixture else "PRODUCTION_PERFORMANCE_METRICS"))},
    ]
    resource_rows = contracts.get("resource_manifest.tsv", [])
    for resource in resource_rows:
        logical_name = resource.get("resource_name") or resource.get("name")
        if not logical_name:
            continue
        descriptor = "; ".join(part for part in [
            safe_display(resource.get("version", "")),
            Path(resource.get("path_or_filename", "")).name if resource.get("path_or_filename") else "",
            f"sha256 {resource.get('checksum')}" if available(resource.get("checksum", "")) else "",
        ] if available(part))
        provenance_rows.append({"field": f"Resource: {logical_name}", "value": descriptor or "Declared"})

    export_rows = []
    seen = set()
    for path in sorted(source_paths):
        if path in seen or not path.is_file():
            continue
        seen.add(path)
        href = local_href(output.parent, path)
        export_rows.append({"name": path.name, "type": path.parent.name, "link": href or "Retained in output tree"})

    figure_dir = output.parent / "figures"
    figures = {"coverage_summary.svg": coverage_svg, "benchmark_f1.svg": bench_svg, "caller_overlap.svg": overlap_svg, "stratified_f1_heatmap.svg": heat_svg, "sv_class_counts.svg": sv_svg, "stage_runtime.svg": runtime_svg}
    for name, svg in figures.items():
        write_figure(figure_dir, name, svg)

    link_buttons = " ".join(filter(None, [configured_link(args.multiqc_link, "MultiQC report"), configured_link(args.vcf2circos_link, "VCF2Circos")])) or '<span class="muted">No companion visualization generated</span>'
    fixture_note = '<div class="notice"><strong>Technical validation report.</strong> This compact fixture demonstrates workflow and reporting behavior. Metrics are not full-WGS performance claims, biological truth claims, clinical findings, or caller rankings.</div>' if fixture else ""
    benchmark_note = "hap.py/vcfeval metrics are genotype/haplotype aware." if analysis == "germline" and any("hap" in r.get("engine", "").lower() for r in benchmarks) else "Somatic metrics use normalized PASS-allele comparison." if analysis == "somatic" and benchmark_enabled else ""
    cnv_empty = "Not run" if not args.run_cnv else ("Production CNV architecture configured; compact fixture not biologically suitable for CNV validation." if fixture else "Not available")
    sv_empty = "Not run" if not args.run_sv else "Not available"
    comparison_empty = "Not run" if not comparison_enabled else "Not available"
    benchmark_empty = "Not run" if not benchmark_enabled else "Not available"
    generated = args.generated_at or datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    style = """
:root{--ink:#17212b;--muted:#5c6976;--line:#dce3e8;--paper:#fff;--wash:#f5f7f8;--blue:#3478b8;--green:#2f7d5b;--orange:#bd642f;--purple:#7654a8;--teal:#267f80;--amber:#a7672d;--indigo:#5367a5}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:#eef2f4;color:var(--ink);font:15px/1.5 Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif}.page{max-width:1260px;margin:0 auto;background:var(--paper);box-shadow:0 0 32px #20304018}header{padding:44px 54px 34px;border-top:7px solid #354552;background:#f9fafb}h1{font:700 clamp(28px,4vw,46px)/1.1 Georgia,serif;margin:0 0 10px;color:#1d2b36}.subtitle{font-size:17px;color:var(--muted)}.meta{margin-top:17px;color:var(--muted);font-size:13px}.notice{margin:0 54px;padding:14px 18px;border-left:4px solid var(--amber);background:#fff7e9}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;padding:28px 54px}.metric{min-height:105px;border:1px solid var(--line);border-radius:7px;padding:14px 16px;background:#fff}.metric span,.metric small{display:block;color:var(--muted)}.metric strong{display:block;font:700 25px/1.2 Georgia,serif;margin:8px 0 3px}.metric small{font-size:11px}nav{position:sticky;top:0;z-index:2;padding:10px 54px;background:#fffefef2;border-block:1px solid var(--line);backdrop-filter:blur(7px);overflow:auto;white-space:nowrap}nav a{color:#344858;text-decoration:none;margin-right:20px;font-size:13px}main{padding:0 54px 48px}section{padding:35px 0;border-bottom:1px solid var(--line);scroll-margin-top:55px}section:last-child{border:0}h2{font:700 25px/1.2 Georgia,serif;margin:0 0 17px;border-left:5px solid var(--section,var(--blue));padding-left:12px}h3{font-size:16px;margin:24px 0 10px}.lead,.muted,.empty{color:var(--muted)}.empty{padding:16px;border:1px dashed #c9d1d7;border-radius:6px;background:var(--wash)}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(350px,1fr));gap:20px}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:6px}table{border-collapse:collapse;width:100%;font-size:13px}th{background:#f1f4f6;text-align:left;color:#334554;font-weight:650}th,td{padding:9px 11px;border-bottom:1px solid var(--line);vertical-align:top}tbody tr:last-child td{border-bottom:0}tbody tr:hover{background:#fafcfd}.figure{margin:18px 0;padding:12px;border:1px solid var(--line);border-radius:6px}.figure svg{display:block;width:100%;height:auto;max-height:480px}.figure svg text{font:12px Inter,system-ui,sans-serif;fill:#273744}.heat-value{font-weight:650}.button{display:inline-block;border:1px solid #667986;border-radius:5px;padding:8px 12px;color:#273b49;text-decoration:none;background:#fff;margin:4px 8px 4px 0}.button:hover{background:#f1f5f6}details{border:1px solid var(--line);border-radius:6px;padding:10px 13px;margin-top:12px}summary{cursor:pointer;font-weight:650}.qc{--section:var(--blue)}.variants{--section:var(--green)}.benchmark{--section:var(--purple)}.comparison{--section:var(--teal)}.svcnv{--section:var(--amber)}.performance{--section:var(--indigo)}.provenance{--section:#53616c}footer{padding:25px 54px;background:#f4f6f7;color:var(--muted);font-size:12px}@media(max-width:700px){header,.metrics,main,nav,footer{padding-left:22px;padding-right:22px}.notice{margin:0 22px}.grid{grid-template-columns:1fr}}@media print{body{background:#fff}.page{box-shadow:none;max-width:none}nav{display:none}section{break-inside:avoid}.button{border:0;padding:0}.table-wrap{overflow:visible}}
"""
    nav = "".join(f'<a href="#{key}">{label}</a>' for key, label in [("overview","Overview"),("qc","QC"),("variants","Small variants"),("benchmark","Benchmark"),("comparison","Concordance"),("stratification","Stratification"),("svcnv","SV/CNV"),("performance","Performance"),("provenance","Provenance")])
    html_doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(title)}</title><style>{style}</style></head><body><div class="page">
<header><h1>{esc(title)}</h1><div class="subtitle">{esc(subtitle)}</div><div class="meta">Generated {esc(generated)} · descriptive reporting layer · no workflow-routing effect</div></header>{fixture_note}<div class="metrics">{''.join(cards)}</div><nav aria-label="Report sections">{nav}</nav><main>
<section id="overview" class="provenance"><h2>Run overview</h2><div class="grid"><div>{table([('field','Field'),('value','Value')],overview_rows)}</div><div><h3>Companion reports</h3><p>{link_buttons}</p><p class="lead">MultiQC retains detailed native QC views; VCF2Circos retains its interactive genomic visualization. Neither is embedded or required to render this report.</p></div></div></section>
<section id="qc" class="qc"><h2>QC and coverage</h2><div class="grid"><div><h3>Coverage</h3>{table([('sample','Sample'),('role','Role'),('mean','Mean depth'),('median','Median depth'),('bases','Region bases')],coverage_rows,'Not run' if not args.run_coverage_qc else 'Not available')}</div><div class="figure">{coverage_svg or '<p class="empty">No meaningful coverage plot available</p>'}</div></div><details><summary>Coverage breadth</summary>{table([('sample','Sample'),('threshold','Threshold'),('breadth','Breadth fraction')],breadth_rows,'Not run' if not args.run_coverage_qc else 'Not available')}</details><div class="grid"><div><h3>Alignment and duplicate QC</h3>{table([('sample','Sample'),('mapped','Mapped'),('paired','Properly paired'),('duplicates','Duplicate rate'),('error','Error rate')],alignment_qc(root),'See the linked MultiQC report')}</div><div><h3>fastp</h3>{table([('sample','Sample'),('reads','Input reads'),('retained','Reads retained'),('q30','Post-filter Q30')],fastp_rows(root))}</div><div><h3>FastQC archive status</h3>{table([('read','Read'),('status','Worst module status')],fastqc_status(root))}</div></div></section>
<section id="variants" class="variants"><h2>Small-variant summary</h2><p class="lead">Counts describe standardized normalized caller outputs. Raw record semantics can differ between callers; these values are labeled rather than treated as interchangeable biological quantities.</p>{table([('caller','Caller'),('version','Version'),('records','Normalized records'),('snv','SNVs'),('indel','Indels'),('pass','PASS'),('samples','Sample/patient'),('benchmark','Benchmark')],variant_rows)}</section>
<section id="benchmark" class="benchmark"><h2>Benchmark accuracy</h2><p class="lead">{esc(benchmark_note)} Fixture values are technical validation results and are not full-WGS accuracy claims.</p><div class="figure">{bench_svg or f'<p class="empty">{esc(benchmark_empty)}</p>'}</div>{table([('caller','Caller'),('class','Class'),('tp','TP'),('fp','FP'),('fn','FN'),('precision','Precision'),('recall','Recall'),('f1','F1')],bench_rows,benchmark_empty)}</section>
<section id="comparison" class="comparison"><h2>Caller concordance and disagreement</h2><p class="lead">Shared and caller-specific are neutral set-membership descriptions—not consensus, validation, confidence, or truth labels.</p><div class="grid"><div>{table([('class','Class'),('category','Agreement category'),('count','Alleles'),('truth','Truth-match fraction')],concord_rows,comparison_empty)}</div><div class="figure">{overlap_svg or f'<p class="empty">{esc(comparison_empty)}</p>'}</div></div><h3>Pairwise allele Jaccard</h3>{table([('pair','Caller pair'),('class','Class'),('jaccard','Jaccard'),('intersection','Intersection'),('union','Union')],jac_rows,comparison_empty)}<h3>Germline genotype concordance</h3>{table([('class','Class'),('concordant','Concordant'),('discordant','Discordant'),('missing','Missing/no-call'),('rate','Callable concordance')],genotype_rows,'Not applicable' if analysis == 'somatic' else comparison_empty)}</section>
<section id="stratification" class="benchmark"><h2>Stratified benchmarking and concordance</h2><p class="lead">Strata may overlap and are non-additive. {'Fixture strata are coordinate-window technical validation strata.' if fixture else 'Configured authoritative stratum names are preserved.'}</p><div class="figure">{heat_svg or '<p class="empty">Not run or not available</p>'}</div>{table([('caller','Caller'),('stratum','Stratum'),('class','Class'),('precision','Precision'),('recall','Recall'),('f1','F1')],strat_rows,'Not run or not available')}<details><summary>Stratified concordance details</summary>{table([('stratum','Stratum'),('class','Class'),('category','Category'),('count','Alleles'),('truth','Truth-match fraction')],strat_concord_rows,'Not available')}<h3>Jaccard by stratum</h3>{table([('stratum','Stratum'),('pair','Caller pair'),('class','Class'),('jaccard','Jaccard')],strat_jac_rows,'Not available')}</details></section>
<section id="svcnv" class="svcnv"><h2>Structural and copy-number variants</h2><div class="grid"><div><h3>Structural variants</h3>{table([('caller','Caller'),('semantics','Output semantics'),('total','Final SVs'),('pass','PASS')],sv_rows,sv_empty)}{table([('class','Class'),('total','Total'),('pass','PASS')],sv_class_rows,sv_empty)}</div><div class="figure">{sv_svg or '<p class="empty">No non-empty SV class plot generated</p>'}</div></div><h3>Copy-number variants</h3>{table([('engine','Engine'),('segments','Segments'),('gains','Gains'),('losses','Losses'),('neutral','Neutral'),('altered','Altered bases'),('largest','Largest event (bp)')],cnv_rows,cnv_empty)}</section>
<section id="performance" class="performance"><h2>Computational performance</h2><p class="lead">Workflow elapsed time is run-level wall time; summed task runtime is task compute time and may exceed elapsed time when tasks overlap. Cached tasks are reported separately and excluded from fresh-compute aggregates. Fixture measurements are technical validation only.</p><div class="grid"><div>{table([('stage','Stage'),('fresh','Fresh tasks'),('cached','Cached'),('runtime','Summed runtime'),('cpu','CPU-hours'),('memory','Peak RSS')],stage_rows,'Not available')}</div><div class="figure">{runtime_svg or '<p class="empty">Not available</p>'}</div></div><h3>Accuracy and compute evidence table</h3><p class="lead">No combined score, ordering, winner, or recommendation is calculated.</p>{table([('caller','Caller'),('precision','Precision'),('recall','Recall'),('f1','F1'),('runtime','Runtime'),('cpu','CPU-hours'),('memory','Peak RAM')],compute_rows,'Not available')}</section>
<section id="provenance" class="provenance"><h2>Provenance and reproducibility</h2>{table([('field','Field'),('value','Value')],provenance_rows)}<h3>Machine-readable exports</h3><p class="lead">The HTML summarizes but does not replace standardized outputs. Links below remain relative and portable within the results directory.</p>{table([('name','File'),('type','Output group'),('link','Relative location')],export_rows,'No machine-readable exports discovered')}</section>
</main><footer>Generated by wgs-pipeline integrated reporting · static offline HTML · source values remain authoritative in independent TSV outputs</footer></div></body></html>"""
    forbidden = ["file://", "/mnt/d/", "/tmp/", str(root), str(output.parent.resolve())]
    hits = [token for token in forbidden if token and token in html_doc]
    if hits:
        raise ReportError(f"Portability validation failed; report contains machine-specific token(s): {', '.join(hits)}")
    if re.search(r'<(?:script|link)[^>]+(?:src|href)=["\']https?://', html_doc, flags=re.I):
        raise ReportError("Portability validation failed; report has a remote rendering dependency")
    return html_doc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samplesheet", type=Path)
    parser.add_argument("--report-id", default="")
    parser.add_argument("--title", default="")
    parser.add_argument("--analysis-type", choices=["germline", "somatic"], required=True)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--callers", default="")
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument("--resource-mode", choices=["fixture", "production"], required=True)
    parser.add_argument("--pipeline-version", default="NA")
    parser.add_argument("--nextflow-version", default="NA")
    parser.add_argument("--benchmark-name", default="")
    parser.add_argument("--benchmark-version", default="")
    parser.add_argument("--multiqc-link", default="")
    parser.add_argument("--vcf2circos-link", default="")
    parser.add_argument("--generated-at", default="")
    for flag in ["run-bqsr", "run-benchmark", "run-vep", "run-sv", "run-cnv", "run-coverage-qc"]:
        parser.add_argument(f"--{flag}", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    document = render(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(document, encoding="utf-8")
    digest = hashlib.sha256(document.encode()).hexdigest()
    print(f"Integrated report: {args.output} ({len(document.encode()):,} bytes; sha256 {digest})")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ReportError, OSError) as error:
        print(f"Integrated report error: {error}", file=sys.stderr)
        raise SystemExit(2)
