#!/usr/bin/env python3
"""Summarize native GATK germline segment VCF or somatic called SEG output."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--format", choices=["germline_vcf", "somatic_seg"], required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--analysis", choices=["germline", "somatic"], required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    segments = gains = losses = neutral = altered = largest = 0
    neutral_defined = args.format == "somatic_seg"
    seg_columns: dict[str, int] | None = None
    with args.input.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#") or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if args.format == "germline_vcf":
                if fields[0] == "#CHROM":
                    continue
                if len(fields) < 8 or not fields[1].isdigit():
                    raise ValueError(f"malformed germline CNV VCF at line {number}")
                info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
                end = int(info.get("END", fields[1]))
                width = end - int(fields[1]) + 1
                kind = info.get("SVTYPE", "")
                is_gain, is_loss, is_neutral = kind == "DUP", kind == "DEL", False
            else:
                if "CONTIG" in fields and "START" in fields and "END" in fields and "CALL" in fields:
                    seg_columns = {name: index for index, name in enumerate(fields)}
                    continue
                if seg_columns is None:
                    raise ValueError("somatic called SEG is missing a CONTIG/START/END/CALL header")
                start_index, end_index, call_index = seg_columns["START"], seg_columns["END"], seg_columns["CALL"]
                if max(start_index, end_index, call_index) >= len(fields) or not fields[start_index].isdigit() or not fields[end_index].isdigit():
                    raise ValueError(f"malformed somatic called SEG at line {number}")
                width = int(fields[end_index]) - int(fields[start_index]) + 1
                call = fields[call_index].strip()
                is_gain, is_loss, is_neutral = call == "+", call == "-", call == "0"
            if width <= 0:
                raise ValueError(f"non-positive CNV segment width at line {number}")
            segments += 1
            gains += int(is_gain)
            losses += int(is_loss)
            neutral += int(is_neutral)
            if is_gain or is_loss:
                altered += width
                largest = max(largest, width)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_or_patient", "analysis_type", "segments", "gain_segments", "loss_segments", "neutral_segments", "total_altered_bases", "largest_event_bp", "reference_assembly", "engine", "engine_version"])
        writer.writerow([args.id, args.analysis, segments, gains, losses, neutral if neutral_defined else "NA", altered, largest if gains + losses else "NA", args.assembly, "GATK CNV", "4.6.2.0"])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as error:
        print(f"CNV summary error: {error}", file=__import__("sys").stderr)
        raise SystemExit(2)
