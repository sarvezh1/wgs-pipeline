#!/usr/bin/env python3
"""Create publication-ready performance tables from a Nextflow trace file."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import re
import sys
import time
from collections import defaultdict
from pathlib import Path


NA = "NA"
REQUIRED_COLUMNS = {
    "hash",
    "name",
    "tag",
    "status",
    "attempt",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "read_bytes",
    "write_bytes",
    "memory",
    "cpus",
    "container",
}
ANALYSES = {"germline", "somatic"}
CALLERS = {
    "germline": {"haplotypecaller", "deepvariant", "compare"},
    "somatic": {"mutect2", "strelka2", "deepsomatic", "compare"},
}

PROCESS_STAGE = {
    "INPUT_CHECK": "input/QC",
    "FASTQC_RAW": "input/QC",
    "FASTQC": "input/QC",
    "FASTP_TRIM": "input/QC",
    "FASTP": "input/QC",
    "FASTQC_TRIMMED": "input/QC",
    "FIXTURE_REFERENCE_PREP": "alignment",
    "FIXTURE_BWA_INDEX": "alignment",
    "BWA_MEM2_ALIGN": "alignment",
    "SORT_BAM": "alignment post-processing",
    "SAMTOOLS_SORT_INDEX": "alignment post-processing",
    "MARKDUP": "alignment post-processing",
    "SAMTOOLS_MARKDUP_QC": "alignment post-processing",
    "BQSR_SUPPORT_PREP": "alignment post-processing",
    "GATK_BQSR": "BQSR",
    "SAMTOOLS_BQSR_INDEX": "BQSR",
    "GATK_HAPLOTYPECALLER": "germline calling",
    "GATK_GENOTYPEGVCFS": "germline calling",
    "GATK_HAPLOTYPECALLER_GVCF": "germline calling",
    "GATK_GENOTYPE_GVCF": "germline calling",
    "DEEPVARIANT": "germline calling",
    "MUTECT2": "somatic calling",
    "MUTECT2_PRODUCTION": "somatic calling",
    "LEARN_READ_ORIENTATION_MODEL": "somatic calling",
    "GET_PILEUP_SUMMARIES": "somatic calling",
    "CALCULATE_CONTAMINATION": "somatic calling",
    "FILTERMUTECTCALLS": "somatic calling",
    "FILTERMUTECTCALLS_PRODUCTION": "somatic calling",
    "STRELKA2_REGIONS": "somatic calling",
    "STRELKA2": "somatic calling",
    "MERGE_STRELKA2_VCFS": "somatic calling",
    "DEEPSOMATIC": "somatic calling",
    "VCF_PROCESSING": "variant processing",
    "VEP_ANNOTATE": "annotation",
    "VEP": "annotation",
    "VCF2CIRCOS_PREP": "visualization",
    "VCF2CIRCOS": "visualization",
    "BENCHMARK_GERMLINE_EXACT": "benchmarking",
    "EXACT_ALLELE_BENCHMARK": "benchmarking",
    "BENCHMARK_GERMLINE_HAPPY": "benchmarking",
    "HAPPY_BENCHMARK": "benchmarking",
    "BENCHMARK_SOMATIC": "benchmarking",
    "STRATIFICATION_PREP": "benchmarking",
    "BENCHMARK_GERMLINE_STRATIFIED": "benchmarking",
    "BENCHMARK_SOMATIC_STRATIFIED": "benchmarking",
    "STRATIFIED_EXACT_BENCHMARK": "benchmarking",
    "COMPARE_GERMLINE": "caller comparison",
    "CALLSET_COMPARISON": "caller comparison",
    "COMPARE_SOMATIC": "caller comparison",
    "COMPARE_STRATIFIED": "caller comparison",
    "STRATIFIED_GERMLINE_COMPARISON": "caller comparison",
    "STRATIFIED_SOMATIC_COMPARISON": "caller comparison",
    "MULTIQC": "reporting",
    "ANALYSIS_REPORT": "reporting",
    "MANTA_REGIONS": "structural variation",
    "MANTA_GERMLINE": "structural variation",
    "MANTA_SOMATIC": "structural variation",
    "SV_SUMMARY": "structural variation",
    "MOSDEPTH": "coverage QC",
    "COVERAGE_SUMMARY": "coverage QC",
    "CNV_COLLECT_READ_COUNTS": "copy-number analysis",
    "GERMLINE_CONTIG_PLOIDY_CASE": "copy-number analysis",
    "GERMLINE_CNV_CALLER_CASE": "copy-number analysis",
    "POSTPROCESS_GERMLINE_CNV": "copy-number analysis",
    "SOMATIC_CNV_DENOISE": "copy-number analysis",
    "CNV_COLLECT_ALLELIC_COUNTS": "copy-number analysis",
    "SOMATIC_CNV_MODEL_SEGMENTS": "copy-number analysis",
    "SOMATIC_CNV_CALL_SEGMENTS": "copy-number analysis",
    "CNV_SUMMARY": "copy-number analysis",
}

PROCESS_CALLER = {
    "GATK_HAPLOTYPECALLER": "haplotypecaller",
    "GATK_GENOTYPEGVCFS": "haplotypecaller",
    "DEEPVARIANT": "deepvariant",
    "MUTECT2": "mutect2",
    "MUTECT2_PRODUCTION": "mutect2",
    "LEARN_READ_ORIENTATION_MODEL": "mutect2",
    "GET_PILEUP_SUMMARIES": "mutect2",
    "CALCULATE_CONTAMINATION": "mutect2",
    "FILTERMUTECTCALLS": "mutect2",
    "FILTERMUTECTCALLS_PRODUCTION": "mutect2",
    "STRELKA2_REGIONS": "strelka2",
    "STRELKA2": "strelka2",
    "MERGE_STRELKA2_VCFS": "strelka2",
    "DEEPSOMATIC": "deepsomatic",
}

DURATION_FACTORS = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}
BYTE_FACTORS = {
    "b": 1,
    "kb": 1024,
    "mb": 1024**2,
    "gb": 1024**3,
    "tb": 1024**4,
}


def fail(message: str) -> None:
    raise ValueError(message)


def is_na(value: str | None) -> bool:
    return value is None or value.strip() in {"", "-", "NA", "N/A", "null"}


def parse_nonnegative_number(value: str, field: str) -> float | None:
    if is_na(value):
        return None
    text = value.strip().replace(",", "")
    try:
        parsed = float(text)
    except ValueError as exc:
        fail(f"Invalid {field} value: {value!r}")
        raise exc
    if parsed < 0:
        fail(f"Impossible negative {field} value: {value!r}")
    return parsed


def parse_duration(value: str, field: str = "duration") -> float | None:
    if is_na(value):
        return None
    text = value.strip().lower()
    if text == "0":
        return 0.0
    if text.startswith("-"):
        fail(f"Impossible negative {field} value: {value!r}")
    matches = list(re.finditer(r"([0-9]+(?:\.[0-9]+)?)\s*(ms|s|m|h|d)", text))
    compact = re.sub(r"\s+", "", text)
    consumed = "".join(match.group(0).replace(" ", "") for match in matches)
    if not matches or consumed != compact:
        fail(f"Invalid {field} value: {value!r}")
    total = sum(float(match.group(1)) * DURATION_FACTORS[match.group(2)] for match in matches)
    if total < 0:
        fail(f"Impossible negative {field} value: {value!r}")
    return total


def parse_bytes(value: str, field: str) -> int | None:
    if is_na(value):
        return None
    text = value.strip().lower().replace(" ", "").replace(",", "")
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)(b|kb|mb|gb|tb)?", text)
    if not match:
        if text.startswith("-"):
            fail(f"Impossible negative {field} value: {value!r}")
        fail(f"Invalid {field} value: {value!r}")
    number = float(match.group(1))
    return int(round(number * BYTE_FACTORS[match.group(2) or "b"]))


def parse_cpu(value: str) -> float | None:
    if is_na(value):
        return None
    parsed = parse_nonnegative_number(value.strip().removesuffix("%"), "CPU utilization")
    return parsed


def fmt_number(value: float | int | None, digits: int = 6) -> str:
    if value is None:
        return NA
    if isinstance(value, int):
        return str(value)
    rendered = f"{value:.{digits}f}".rstrip("0").rstrip(".")
    return rendered if rendered else "0"


def process_name(task_name: str) -> str:
    base = task_name.split(" (", 1)[0]
    return base.rsplit(":", 1)[-1]


def infer_caller(process: str, tag: str) -> str:
    if process in PROCESS_CALLER:
        return PROCESS_CALLER[process]
    lower = tag.lower()
    for caller in ("haplotypecaller", "deepvariant", "mutect2", "strelka2", "deepsomatic"):
        if re.search(rf"(^|[:_\-]){caller}($|[:_\-])", lower):
            return caller
    return NA


def infer_sample(tag: str, caller: str) -> str:
    if is_na(tag):
        return NA
    token = tag.split(":", 1)[0].strip()
    if not token or token.lower() in {
        "input validation",
        "fixture_reference",
        "fixture_bwa_index",
        "fixture reference preparation",
        "fixture bwa-mem2 indexing",
        "strelka2 call regions",
        "multiqc",
    }:
        return NA
    if token.lower() == caller:
        return NA
    return token


def container_version(identity: str) -> str:
    if is_na(identity):
        return NA
    if "@" in identity:
        return identity.rsplit("@", 1)[1]
    final = identity.rsplit("/", 1)[-1]
    return final.rsplit(":", 1)[1] if ":" in final else NA


def read_trace(path: Path) -> list[dict[str, str]]:
    for _ in range(25):
        if path.exists() and path.stat().st_size > 0:
            break
        time.sleep(0.2)
    if not path.exists():
        fail(f"Nextflow trace file is missing: {path}")
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                fail(f"Nextflow trace has no header: {path}")
            missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames))
            if missing:
                fail("Nextflow trace is missing required column(s): " + ", ".join(missing))
            rows = list(reader)
    except UnicodeDecodeError as exc:
        fail(f"Nextflow trace is not valid UTF-8: {path}")
        raise exc
    if not rows:
        fail(f"Nextflow trace contains no task records: {path}")
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows: list[dict], key: str) -> list[dict[str, str]]:
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row[key] != NA:
            groups[row[key]].append(row)
    output = []
    for label in sorted(groups):
        members = groups[label]
        fresh = [row for row in members if row["status"] == "COMPLETED"]
        cached = sum(row["status"] == "CACHED" for row in members)
        failed = sum(row["status"] in {"FAILED", "ABORTED"} for row in members)
        wall = [row["realtime_seconds_raw"] for row in fresh if row["realtime_seconds_raw"] is not None]
        cpu_hours = [row["cpu_hours_raw"] for row in fresh if row["cpu_hours_raw"] is not None]
        rss = [row["peak_rss_bytes_raw"] for row in fresh if row["peak_rss_bytes_raw"] is not None]
        requested = [row["requested_memory_bytes_raw"] for row in fresh if row["requested_memory_bytes_raw"] is not None]
        output.append(
            {
                key: label,
                "executed_tasks": str(len(fresh)),
                "cached_tasks": str(cached),
                "failed_tasks": str(failed),
                "summed_task_wall_seconds": fmt_number(sum(wall)) if wall else NA,
                "summed_task_cpu_hours": fmt_number(sum(cpu_hours)) if cpu_hours else NA,
                "max_peak_rss_bytes": fmt_number(max(rss)) if rss else NA,
                "max_requested_memory_bytes": fmt_number(max(requested)) if requested else NA,
                "max_task_duration_seconds": fmt_number(max(wall)) if wall else NA,
            }
        )
    return output


def component_summary(rows: list[dict]) -> list[dict[str, str]]:
    preprocessing = {"input/QC", "alignment", "alignment post-processing", "BQSR"}
    caller = {"germline calling", "somatic calling"}
    sections = {
        "shared_preprocessing": preprocessing,
        "caller_execution": caller,
        "downstream_analysis": set(PROCESS_STAGE.values()) - preprocessing - caller,
    }
    staged = []
    for section, stages in sections.items():
        for row in rows:
            if row["stage"] in stages:
                copy = dict(row)
                copy["component"] = section
                staged.append(copy)
    return aggregate(staged, "component")


def available_memory() -> int | None:
    path = Path("/proc/meminfo")
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MemTotal:"):
            return int(line.split()[1]) * 1024
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--analysis", required=True, choices=sorted(ANALYSES))
    parser.add_argument("--caller", required=True)
    parser.add_argument("--sample", default=NA)
    parser.add_argument("--resource-mode", required=True, choices=("fixture", "production"))
    parser.add_argument("--reference-assembly", default=NA)
    parser.add_argument("--bqsr", required=True, choices=("true", "false"))
    parser.add_argument("--benchmark", required=True, choices=("true", "false"))
    parser.add_argument("--benchmark-mode", default=NA)
    parser.add_argument("--benchmark-name", default=NA)
    parser.add_argument("--benchmark-version", default=NA)
    parser.add_argument("--comparison", required=True, choices=("true", "false"))
    parser.add_argument("--stratification", required=True, choices=("true", "false"))
    parser.add_argument("--vep", default="false", choices=("true", "false"))
    parser.add_argument("--visualization", default="false", choices=("true", "false"))
    parser.add_argument("--multiqc", default="false", choices=("true", "false"))
    parser.add_argument("--sv", default="false", choices=("true", "false"))
    parser.add_argument("--cnv", default="false", choices=("true", "false"))
    parser.add_argument("--coverage-qc", default="false", choices=("true", "false"))
    parser.add_argument("--ploidy-configured", default="false", choices=("true", "false"))
    parser.add_argument("--executor", default=NA)
    parser.add_argument("--container-engine", default=NA)
    parser.add_argument("--nextflow-version", default=NA)
    parser.add_argument("--pipeline-version", default=NA)
    parser.add_argument("--workflow-start", required=True)
    parser.add_argument("--workflow-completion", required=True)
    parser.add_argument("--workflow-elapsed-seconds", required=True)
    parser.add_argument("--workflow-success", required=True, choices=("true", "false"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    caller = args.caller.lower()
    if caller not in CALLERS[args.analysis]:
        fail(
            f"Unsupported caller metadata {args.caller!r} for analysis {args.analysis!r}; "
            f"expected one of {', '.join(sorted(CALLERS[args.analysis]))}"
        )
    workflow_elapsed = parse_nonnegative_number(args.workflow_elapsed_seconds, "workflow elapsed time")
    source_scope = "TECHNICAL_VALIDATION_METRICS" if args.resource_mode == "fixture" else "RUN_PERFORMANCE_METRICS"
    raw_rows = read_trace(args.trace)
    parsed_rows = []
    containers = set()
    for raw in raw_rows:
        process = process_name(raw["name"])
        stage = PROCESS_STAGE.get(process, "other")
        status = raw["status"].strip().upper()
        realtime = parse_duration(raw["realtime"], "realtime")
        cpu = parse_cpu(raw["%cpu"])
        requested_cpus = parse_nonnegative_number(raw["cpus"], "requested CPUs")
        requested_memory = parse_bytes(raw["memory"], "requested memory")
        peak_rss = parse_bytes(raw["peak_rss"], "peak RSS")
        peak_vmem = parse_bytes(raw["peak_vmem"], "peak virtual memory")
        read_bytes = parse_bytes(raw["read_bytes"], "read I/O")
        write_bytes = parse_bytes(raw["write_bytes"], "write I/O")
        attempt = parse_nonnegative_number(raw["attempt"], "attempt")
        if attempt is not None and not attempt.is_integer():
            fail(f"Task attempt must be an integer: {raw['attempt']!r}")
        if attempt is not None and attempt < 1:
            fail(f"Task attempt must be at least 1: {raw['attempt']!r}")
        if requested_cpus is not None and requested_cpus <= 0:
            fail(f"Requested CPUs must be positive: {raw['cpus']!r}")
        if requested_memory is not None and requested_memory <= 0:
            fail(f"Requested memory must be positive: {raw['memory']!r}")
        task_cpu_hours = None if realtime is None or cpu is None else realtime * cpu / 100.0 / 3600.0
        inferred_caller = infer_caller(process, raw["tag"])
        identity = raw["container"].strip() if not is_na(raw["container"]) else NA
        if identity != NA:
            containers.add(identity)
        parsed_rows.append(
            {
                "process": process,
                "task": raw["name"],
                "sample_or_patient": infer_sample(raw["tag"], inferred_caller),
                "analysis_type": args.analysis,
                "caller": inferred_caller,
                "stage": stage,
                "attempt": fmt_number(int(attempt)) if attempt is not None else NA,
                "status": status,
                "cached": "true" if status == "CACHED" else "false",
                "realtime_seconds": fmt_number(realtime),
                "cpu_utilization_percent": fmt_number(cpu),
                "requested_cpus": fmt_number(requested_cpus),
                "requested_memory_bytes": fmt_number(requested_memory),
                "peak_rss_bytes": fmt_number(peak_rss),
                "peak_vmem_bytes": fmt_number(peak_vmem),
                "read_bytes": fmt_number(read_bytes),
                "write_bytes": fmt_number(write_bytes),
                "task_compute_cpu_hours": fmt_number(task_cpu_hours),
                "workdir_id": raw["hash"].strip() if not is_na(raw["hash"]) else NA,
                "container_identity": identity,
                "container_version_or_digest": container_version(identity),
                "metric_scope": source_scope,
                "realtime_seconds_raw": realtime,
                "cpu_hours_raw": task_cpu_hours,
                "peak_rss_bytes_raw": peak_rss,
                "requested_memory_bytes_raw": requested_memory,
            }
        )

    public_fields = [
        "process", "task", "sample_or_patient", "analysis_type", "caller", "stage", "attempt",
        "status", "cached", "realtime_seconds", "cpu_utilization_percent", "requested_cpus",
        "requested_memory_bytes", "peak_rss_bytes", "peak_vmem_bytes", "read_bytes", "write_bytes",
        "task_compute_cpu_hours", "workdir_id", "container_identity", "container_version_or_digest",
        "metric_scope",
    ]
    public_rows = [{field: row[field] for field in public_fields} for row in parsed_rows]
    write_tsv(args.output_dir / "performance_summary.tsv", public_fields, public_rows)

    aggregate_fields = [
        "stage", "executed_tasks", "cached_tasks", "failed_tasks", "summed_task_wall_seconds",
        "summed_task_cpu_hours", "max_peak_rss_bytes", "max_requested_memory_bytes",
        "max_task_duration_seconds",
    ]
    write_tsv(args.output_dir / "stage_summary.tsv", aggregate_fields, aggregate(parsed_rows, "stage"))
    caller_fields = ["caller"] + aggregate_fields[1:]
    caller_stage_rows = [
        row for row in parsed_rows if row["stage"] in {"germline calling", "somatic calling"}
    ]
    write_tsv(args.output_dir / "caller_performance.tsv", caller_fields, aggregate(caller_stage_rows, "caller"))
    component_fields = ["component"] + aggregate_fields[1:]
    write_tsv(
        args.output_dir / "compute_component_summary.tsv",
        component_fields,
        component_summary(parsed_rows),
    )

    executed = sum(row["status"] == "COMPLETED" for row in parsed_rows)
    cached = sum(row["status"] == "CACHED" for row in parsed_rows)
    failed = sum(row["status"] in {"FAILED", "ABORTED"} for row in parsed_rows)
    workflow_fields = [
        "workflow_start", "workflow_completion", "workflow_elapsed_seconds", "workflow_success",
        "executed_successful_tasks", "cached_tasks", "failed_tasks", "analysis_type", "caller_selection",
        "bqsr_enabled", "resource_mode", "benchmark_enabled", "comparison_enabled",
        "stratification_enabled", "sv_enabled", "cnv_enabled", "coverage_qc_enabled", "ploidy_configured",
        "vep_enabled", "visualization_enabled", "multiqc_enabled", "executor",
        "nextflow_version", "operating_system", "cpu_architecture", "available_cpus",
        "available_memory_bytes", "container_engine", "metric_scope",
    ]
    workflow_row = {
        "workflow_start": args.workflow_start,
        "workflow_completion": args.workflow_completion,
        "workflow_elapsed_seconds": fmt_number(workflow_elapsed),
        "workflow_success": args.workflow_success,
        "executed_successful_tasks": str(executed),
        "cached_tasks": str(cached),
        "failed_tasks": str(failed),
        "analysis_type": args.analysis,
        "caller_selection": caller,
        "bqsr_enabled": args.bqsr,
        "resource_mode": args.resource_mode,
        "benchmark_enabled": args.benchmark,
        "comparison_enabled": args.comparison,
        "stratification_enabled": args.stratification,
        "sv_enabled": args.sv,
        "cnv_enabled": args.cnv,
        "coverage_qc_enabled": args.coverage_qc,
        "ploidy_configured": args.ploidy_configured,
        "vep_enabled": args.vep,
        "visualization_enabled": args.visualization,
        "multiqc_enabled": args.multiqc,
        "executor": args.executor,
        "nextflow_version": args.nextflow_version,
        "operating_system": platform.system() or NA,
        "cpu_architecture": platform.machine() or NA,
        "available_cpus": str(os.cpu_count()) if os.cpu_count() is not None else NA,
        "available_memory_bytes": fmt_number(available_memory()),
        "container_engine": args.container_engine,
        "metric_scope": source_scope,
    }
    write_tsv(args.output_dir / "workflow_summary.tsv", workflow_fields, [workflow_row])

    provenance_fields = [
        "pipeline_version", "nextflow_version", "executor", "container_engine", "container_identities",
        "cpu_architecture", "available_cpus", "available_memory_bytes", "resource_mode",
        "reference_assembly", "analysis_type", "selected_callers", "bqsr_enabled", "benchmark_enabled",
        "benchmark_mode", "benchmark_name", "benchmark_version", "comparison_enabled",
        "stratification_enabled", "sv_enabled", "cnv_enabled", "coverage_qc_enabled", "ploidy_configured",
        "vep_enabled", "visualization_enabled", "multiqc_enabled",
        "trace_metric_semantics", "metric_scope",
    ]
    selected = (
        "haplotypecaller,deepvariant"
        if args.analysis == "germline" and caller == "compare"
        else "mutect2,strelka2,deepsomatic"
        if args.analysis == "somatic" and caller == "compare"
        else caller
    )
    provenance_row = {
        "pipeline_version": args.pipeline_version,
        "nextflow_version": args.nextflow_version,
        "executor": args.executor,
        "container_engine": args.container_engine,
        "container_identities": ";".join(sorted(containers)) if containers else NA,
        "cpu_architecture": workflow_row["cpu_architecture"],
        "available_cpus": workflow_row["available_cpus"],
        "available_memory_bytes": workflow_row["available_memory_bytes"],
        "resource_mode": args.resource_mode,
        "reference_assembly": args.reference_assembly,
        "analysis_type": args.analysis,
        "selected_callers": selected,
        "bqsr_enabled": args.bqsr,
        "benchmark_enabled": args.benchmark,
        "benchmark_mode": args.benchmark_mode if args.benchmark == "true" else NA,
        "benchmark_name": args.benchmark_name if args.benchmark == "true" else NA,
        "benchmark_version": args.benchmark_version if args.benchmark == "true" else NA,
        "comparison_enabled": args.comparison,
        "stratification_enabled": args.stratification,
        "sv_enabled": args.sv,
        "cnv_enabled": args.cnv,
        "coverage_qc_enabled": args.coverage_qc,
        "ploidy_configured": args.ploidy_configured,
        "vep_enabled": args.vep,
        "visualization_enabled": args.visualization,
        "multiqc_enabled": args.multiqc,
        "trace_metric_semantics": (
            "Task metrics are parsed from the native Nextflow trace; displayed memory units use 1024-byte "
            "multipliers; CPU-hours equal realtime multiplied by trace CPU utilization; cached tasks are "
            "reported but excluded from fresh-compute aggregates"
        ),
        "metric_scope": source_scope,
    }
    write_tsv(args.output_dir / "performance_provenance.tsv", provenance_fields, [provenance_row])
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    except OSError as exc:
        print(f"ERROR: Performance report filesystem operation failed: {exc}", file=sys.stderr)
        sys.exit(2)
