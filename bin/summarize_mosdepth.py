#!/usr/bin/env python3
"""Create stable coverage tables from native mosdepth summaries and thresholds."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def fmt(value: float) -> str:
    return f"{value:.8f}".rstrip("0").rstrip(".")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--thresholds", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--patient", required=True)
    parser.add_argument("--role", required=True)
    parser.add_argument("--analysis", choices=["germline", "somatic"], required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()
    with args.summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"chrom", "length", "bases", "mean"}.issubset(rows[0]):
        raise ValueError("mosdepth summary is empty or missing chrom/length/bases/mean columns")
    totals = [row for row in rows if row["chrom"] == "total"]
    if len(totals) != 1:
        raise ValueError("mosdepth summary must contain exactly one total row")
    chromosome_rows = [row for row in rows if row["chrom"] != "total" and not row["chrom"].endswith("_region")]
    weighted_depth = 0.0
    weighted_bases = 0
    with gzip.open(args.regions, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            row = raw.rstrip("\n").split("\t")
            if len(row) < 5 or not row[1].isdigit() or not row[2].isdigit():
                raise ValueError(f"malformed mosdepth regions row at line {line_number}")
            width = int(row[2]) - int(row[1])
            if width <= 0:
                raise ValueError(f"non-positive mosdepth regions interval at line {line_number}")
            try:
                mean = float(row[4])
            except ValueError as error:
                raise ValueError(f"invalid mosdepth region mean at line {line_number}") from error
            if mean < 0:
                raise ValueError(f"negative mosdepth region mean at line {line_number}")
            weighted_bases += width
            weighted_depth += width * mean
    with gzip.open(args.thresholds, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 5 or header[:4] != ["#chrom", "start", "end", "region"]:
            raise ValueError("mosdepth thresholds output has an unexpected header")
        try:
            threshold_values = [int(value.removesuffix("X").removesuffix("x")) for value in header[4:]]
        except ValueError as error:
            raise ValueError("mosdepth thresholds header contains a non-integer threshold") from error
        covered = {value: 0 for value in threshold_values}
        region_bases = 0
        for line_number, row in enumerate(reader, 2):
            if len(row) != len(header) or not row[1].isdigit() or not row[2].isdigit():
                raise ValueError(f"malformed mosdepth threshold row at line {line_number}")
            width = int(row[2]) - int(row[1])
            if width <= 0:
                raise ValueError(f"non-positive mosdepth threshold interval at line {line_number}")
            region_bases += width
            for value, count in zip(threshold_values, row[4:]):
                if not count.isdigit() or int(count) > width:
                    raise ValueError(f"impossible threshold count at line {line_number}")
                covered[value] += int(count)
    if region_bases == 0:
        raise ValueError("mosdepth threshold regions contain zero bases")
    if weighted_bases != region_bases:
        raise ValueError("mosdepth regions and threshold outputs cover different numbers of bases")
    common = [args.sample, args.patient, args.role, args.analysis]
    with args.prefix.with_suffix(".coverage_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "patient_id", "role", "analysis_type", "mean_depth", "median_depth", "region_bases", "engine", "engine_version", "reference_assembly"])
        writer.writerow([*common, fmt(weighted_depth / weighted_bases), "NA", region_bases, "mosdepth", "0.3.14", args.assembly])
    with args.prefix.with_suffix(".coverage_threshold_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "patient_id", "role", "analysis_type", "threshold", "covered_bases", "region_bases", "breadth_fraction", "engine", "engine_version", "reference_assembly"])
        for value in threshold_values:
            writer.writerow([*common, value, covered[value], region_bases, fmt(covered[value] / region_bases), "mosdepth", "0.3.14", args.assembly])
    with args.prefix.with_suffix(".coverage_chromosome_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "patient_id", "role", "analysis_type", "chrom", "reference_length", "observed_bases", "mean_depth", "engine", "engine_version", "reference_assembly"])
        for row in chromosome_rows:
            writer.writerow([*common, row["chrom"], row["length"], row["bases"], row["mean"], "mosdepth", "0.3.14", args.assembly])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, StopIteration) as error:
        print(f"Coverage summary error: {error}", file=__import__("sys").stderr)
        raise SystemExit(2)
