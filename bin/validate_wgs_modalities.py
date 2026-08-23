#!/usr/bin/env python3
"""Validate Block 7 modality configuration and emit public-safe provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
from pathlib import Path


def fail(message: str) -> None:
    raise ValueError(f"WGS modality configuration error: {message}")


def boolean(value: str, name: str) -> bool:
    text = value.strip().lower()
    if text not in {"true", "false"}:
        fail(f"{name} must be true or false, observed {value!r}")
    return text == "true"


def thresholds(value: str) -> list[int]:
    if not value.strip():
        fail("coverage_thresholds must be a comma-delimited list of positive integers")
    fields = value.split(",")
    if any(not re.fullmatch(r"[1-9][0-9]*", item.strip()) for item in fields):
        fail("coverage_thresholds must contain unique positive integers, for example 10,20,30")
    parsed = [int(item) for item in fields]
    if len(parsed) != len(set(parsed)):
        fail("coverage_thresholds contains duplicate values")
    if any(item > 100000 for item in parsed):
        fail("coverage_thresholds contains an implausibly large value (>100000)")
    return sorted(set([1, *parsed]))


def require_file(value: str, label: str) -> Path:
    if not value:
        fail(f"missing required {label}")
    path = Path(value)
    if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
        fail(f"{label} is missing, empty, or unreadable: {value}")
    return path


def require_directory(value: str, label: str) -> Path:
    if not value:
        fail(f"missing required {label}")
    path = Path(value)
    if not path.is_dir() or not os.access(path, os.R_OK) or not any(path.iterdir()):
        fail(f"{label} directory is missing, empty, or unreadable: {value}")
    return path


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def read_ploidy(path: Path | None, samples: list[dict[str, str]]) -> dict[str, str]:
    sample_ids = {row["sample_id"] for row in samples}
    if path is None:
        return {sample: "unknown" for sample in sample_ids}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample_id", "sex_chromosome_configuration"]:
            fail("ploidy_config header must be sample_id<TAB>sex_chromosome_configuration")
        configured: dict[str, str] = {}
        for line_number, row in enumerate(reader, 2):
            sample = (row.get("sample_id") or "").strip()
            value = (row.get("sex_chromosome_configuration") or "").strip()
            if not sample or sample in configured:
                fail(f"ploidy_config line {line_number} has an empty or duplicate sample_id")
            if value not in {"XX", "XY", "unknown"}:
                fail(f"ploidy_config line {line_number} has unsupported configuration {value!r}; expected XX, XY, or unknown")
            configured[sample] = value
    missing = sorted(sample_ids - configured.keys())
    extra = sorted(configured.keys() - sample_ids)
    if missing or extra:
        fail(f"ploidy_config sample IDs do not match the samplesheet (missing={missing or 'none'}, extra={extra or 'none'})")
    return configured


def validate_interval_list(path: Path, minimum_intervals: int = 1) -> dict[str, int]:
    sequence_lengths: dict[str, int] = {}
    intervals = 0
    span = 0
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith("@SQ"):
                fields = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
                if "SN" in fields and fields.get("LN", "").isdigit():
                    sequence_lengths[fields["SN"]] = int(fields["LN"])
                continue
            if not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or not fields[1].isdigit() or not fields[2].isdigit():
                fail(f"malformed CNV interval list at line {number}")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 1 or end < start:
                fail(f"invalid CNV interval coordinates at line {number}")
            if contig not in sequence_lengths or end > sequence_lengths[contig]:
                fail(f"CNV interval at line {number} is incompatible with its sequence dictionary")
            intervals += 1
            span += end - start + 1
    if intervals < minimum_intervals:
        fail(f"CNV interval list has {intervals} intervals; at least {minimum_intervals} are required")
    if minimum_intervals >= 10000 and span < 10_000_000:
        fail("germline gCNV intervals span less than the documented 10 Mb minimum")
    return sequence_lengths


def manifest_metadata(path: Path, resource: Path, expected_assembly: str) -> tuple[str, str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in rows if Path(row.get("path_or_filename", "")).name == resource.name]
    if not matches:
        fail(f"production resource is absent from resource manifest: {resource.name}")
    row = matches[0]
    declared_assembly = (row.get("assembly") or "").strip()
    if declared_assembly and declared_assembly != expected_assembly:
        fail(f"resource {resource.name} declares assembly {declared_assembly!r}, expected {expected_assembly!r}")
    return row.get("version") or "NA", row.get("source") or "NA", row.get("checksum") or "NA"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--resource-mode", required=True)
    parser.add_argument("--resource-manifest", required=True)
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--reference-fai", default="")
    parser.add_argument("--run-sv", required=True)
    parser.add_argument("--run-cnv", required=True)
    parser.add_argument("--run-coverage", required=True)
    parser.add_argument("--coverage-thresholds", required=True)
    parser.add_argument("--ploidy-config", default="")
    parser.add_argument("--cnv-intervals", default="")
    parser.add_argument("--cnv-ploidy-model", default="")
    parser.add_argument("--cnv-germline-model", default="")
    parser.add_argument("--cnv-read-count-pon", default="")
    parser.add_argument("--cnv-common-sites", default="")
    parser.add_argument("--provenance-out", required=True)
    parser.add_argument("--ploidy-out", required=True)
    args = parser.parse_args()

    run_sv = boolean(args.run_sv, "run_sv")
    run_cnv = boolean(args.run_cnv, "run_cnv")
    run_coverage = boolean(args.run_coverage, "run_coverage_qc")
    parsed_thresholds = thresholds(args.coverage_thresholds)
    if args.analysis not in {"germline", "somatic"}:
        fail(f"unsupported analysis mode for SV/CNV/coverage: {args.analysis!r}")

    samplesheet = require_file(args.samplesheet, "samplesheet")
    manifest = require_file(args.resource_manifest, "resource manifest")
    samples = read_samples(samplesheet)
    ploidy_path = require_file(args.ploidy_config, "ploidy config") if args.ploidy_config else None
    configurations = read_ploidy(ploidy_path, samples)

    resources: list[tuple[str, Path]] = []
    if run_cnv:
        if args.resource_mode != "production":
            fail("CNV requires production mode with external intervals and trained/model resources; chr20 fixtures cannot train a legitimate CNV model or PoN")
        intervals = require_file(args.cnv_intervals, "CNV interval list")
        resources.append(("cnv_intervals", intervals))
        if args.analysis == "germline":
            interval_sequences = validate_interval_list(intervals, 10000)
            resources.extend([
                ("cnv_ploidy_model", require_directory(args.cnv_ploidy_model, "germline contig-ploidy model")),
                ("cnv_germline_model", require_directory(args.cnv_germline_model, "germline gCNV model")),
            ])
        else:
            interval_sequences = validate_interval_list(intervals)
            pon = require_file(args.cnv_read_count_pon, "somatic CNV read-count PoN")
            with pon.open("rb") as handle:
                if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                    fail("somatic CNV read-count PoN is not an HDF5 file")
            sites = require_file(args.cnv_common_sites, "somatic CNV common-sites VCF")
            require_file(f"{sites}.tbi", "somatic CNV common-sites VCF index")
            resources.extend([("cnv_read_count_pon", pon), ("cnv_common_sites", sites)])

        reference_fai = require_file(args.reference_fai, "production reference FAI")
        reference_sequences = {}
        for number, line in enumerate(reference_fai.read_text(encoding="utf-8").splitlines(), 1):
            fields = line.split("\t")
            if len(fields) < 2 or not fields[1].isdigit():
                fail(f"malformed production reference FAI at line {number}")
            reference_sequences[fields[0]] = int(fields[1])
        incompatible = [name for name, length in interval_sequences.items() if reference_sequences.get(name) != length]
        if incompatible:
            fail("CNV interval dictionary is incompatible with the production reference FAI: " + ",".join(incompatible[:5]))

        project = Path(args.project_dir).resolve()
        for label, resource in resources:
            try:
                resource.resolve().relative_to(project / "test_data")
            except ValueError:
                pass
            else:
                fail(f"production {label} must be external and cannot use fixture-only test_data")
            manifest_metadata(manifest, resource, args.assembly)

    with Path(args.provenance_out).open("w", newline="", encoding="utf-8") as handle:
        fields = ["resource", "enabled", "engine", "engine_version", "container", "assembly", "path", "release", "source", "checksum", "semantics"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow({"resource": "structural_variants", "enabled": str(run_sv).lower(), "engine": "Manta", "engine_version": "1.6.0", "container": "quay.io/biocontainers/manta:1.6.0--h9ee0642_3", "assembly": args.assembly, "path": "NA", "release": "1.6.0", "source": "https://github.com/Illumina/manta", "checksum": "NA", "semantics": "native diploidSV or somaticSV final output"})
        writer.writerow({"resource": "coverage", "enabled": str(run_coverage).lower(), "engine": "mosdepth", "engine_version": "0.3.14", "container": "quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_1", "assembly": args.assembly, "path": "NA", "release": "0.3.14", "source": "https://github.com/brentp/mosdepth", "checksum": "NA", "semantics": f"reporting thresholds={','.join(map(str, parsed_thresholds))}; no per-base output"})
        ploidy_checksum = f"sha256:{hashlib.sha256(ploidy_path.read_bytes()).hexdigest()}" if ploidy_path else "NA"
        writer.writerow({"resource": "ploidy_configuration", "enabled": "true", "engine": "pipeline configuration", "engine_version": "1", "container": "NA", "assembly": args.assembly, "path": ploidy_path.name if ploidy_path else "NA", "release": "NA", "source": "user-supplied" if ploidy_path else "default unknown sex-chromosome configuration", "checksum": ploidy_checksum, "semantics": "provenance/model declaration; does not alter small-variant commands"})
        for label, resource in resources:
            version, source, checksum = manifest_metadata(manifest, resource, args.assembly)
            writer.writerow({"resource": label, "enabled": "true", "engine": "GATK CNV", "engine_version": "4.6.2.0", "container": "broadinstitute/gatk:4.6.2.0", "assembly": args.assembly, "path": resource.name, "release": version, "source": source, "checksum": checksum, "semantics": "external production resource"})

    with Path(args.ploidy_out).open("w", newline="", encoding="utf-8") as handle:
        fields = ["sample_id", "patient_id", "analysis_type", "sex_chromosome_configuration", "autosomal_ploidy", "nonpar_x_ploidy", "nonpar_y_ploidy", "par_ploidy", "application_scope", "reference_assembly"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in samples:
            config = configurations[row["sample_id"]]
            if args.analysis == "somatic":
                values = ("NA", "NA", "NA", "NA", "recorded_not_imposed_on_somatic_callers")
            elif config == "XX":
                values = ("2", "2", "0", "2", "configuration_only; small-variant per-region ploidy routing deferred")
            elif config == "XY":
                values = ("2", "1", "1", "2", "configuration_only; small-variant per-region ploidy routing deferred")
            else:
                values = ("2", "NA", "NA", "NA", "autosomal_model_only; sex-chromosome configuration unknown")
            writer.writerow(dict(zip(fields, [row["sample_id"], row["patient_id"], args.analysis, config, *values, args.assembly])))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        print(error, file=__import__("sys").stderr)
        raise SystemExit(2)
