#!/usr/bin/env python3
"""Summarize a native final Manta VCF without small-variant normalization."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path


CLASSES = ("DEL", "DUP", "INS", "INV", "BND", "OTHER")
SIZE_BINS = ("50-99", "100-499", "500-999", "1000-9999", ">=10000", "<50", "UNDEFINED")


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def info_map(text: str) -> dict[str, str]:
    result = {}
    for field in text.split(";"):
        key, _, value = field.partition("=")
        result[key] = value
    return result


def classify(alt: str, info: dict[str, str]) -> str:
    kind = info.get("SVTYPE", "").upper()
    if not kind and ("[" in alt or "]" in alt):
        kind = "BND"
    return kind if kind in CLASSES[:-1] else "OTHER"


def size_of(kind: str, pos: int, info: dict[str, str]) -> int | None:
    if kind == "BND":
        return None
    raw = info.get("SVLEN", "").split(",", 1)[0]
    if re.fullmatch(r"-?[0-9]+", raw):
        return abs(int(raw))
    end = info.get("END", "")
    if kind in {"DEL", "DUP", "INV"} and end.isdigit() and int(end) >= pos:
        return int(end) - pos + 1
    return None


def size_bin(size: int | None) -> str:
    if size is None:
        return "UNDEFINED"
    if size < 50:
        return "<50"
    if size < 100:
        return "50-99"
    if size < 500:
        return "100-499"
    if size < 1000:
        return "500-999"
    if size < 10000:
        return "1000-9999"
    return ">=10000"


def write(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--fai", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--patient", required=True)
    parser.add_argument("--analysis", choices=["germline", "somatic"], required=True)
    parser.add_argument("--semantics", choices=["diploidSV", "somaticSV"], required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()
    if not args.vcf.is_file() or args.vcf.stat().st_size == 0:
        raise ValueError(f"Manta final VCF is missing or empty: {args.vcf}")
    contigs = {fields[0]: int(fields[1]) for line in args.fai.read_text().splitlines() if line for fields in [line.split("\t")] if len(fields) >= 2 and fields[1].isdigit()}
    if not contigs:
        raise ValueError("reference index is empty or malformed")
    total = passed = 0
    classes, pass_classes, sizes, pass_sizes, chromosomes, pass_chromosomes = (Counter() for _ in range(6))
    saw_header = False
    with open_text(args.vcf) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith("#CHROM\t"):
                saw_header = True
                continue
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[1].isdigit() or int(fields[1]) < 1:
                raise ValueError(f"Malformed Manta VCF record at line {number}")
            chrom, pos, _id, _ref, alt, _qual, filt, info_text = fields[:8]
            if chrom not in contigs:
                raise ValueError(f"Manta VCF contig {chrom!r} is absent from the reference index")
            if int(pos) > contigs[chrom]:
                raise ValueError(f"Manta VCF position at line {number} exceeds reference contig length")
            info = info_map(info_text)
            if info.get("END", "").isdigit() and int(info["END"]) > contigs[chrom]:
                raise ValueError(f"Manta VCF END at line {number} exceeds reference contig length")
            kind = classify(alt, info)
            bin_name = size_bin(size_of(kind, int(pos), info))
            is_pass = filt == "PASS"
            total += 1
            passed += int(is_pass)
            classes[kind] += 1
            sizes[bin_name] += 1
            chromosomes[chrom] += 1
            if is_pass:
                pass_classes[kind] += 1
                pass_sizes[bin_name] += 1
                pass_chromosomes[chrom] += 1
    if not saw_header:
        raise ValueError("Malformed Manta VCF: missing #CHROM header")
    common = {"sample_or_patient": args.sample if args.analysis == "germline" else args.patient, "analysis_type": args.analysis, "output_semantics": args.semantics, "engine": "Manta", "engine_version": "1.6.0", "reference_assembly": args.assembly}
    write(Path(f"{args.prefix}.sv_summary.tsv"), [*common, "total_sv", "pass_sv"], [{**common, "total_sv": total, "pass_sv": passed}])
    write(Path(f"{args.prefix}.sv_class_summary.tsv"), [*common, "sv_class", "total", "pass"], [{**common, "sv_class": kind, "total": classes[kind], "pass": pass_classes[kind]} for kind in CLASSES])
    write(Path(f"{args.prefix}.sv_size_summary.tsv"), [*common, "size_bin_bp", "total", "pass"], [{**common, "size_bin_bp": value, "total": sizes[value], "pass": pass_sizes[value]} for value in SIZE_BINS])
    write(Path(f"{args.prefix}.sv_chromosome_summary.tsv"), [*common, "chrom", "total", "pass"], [{**common, "chrom": chrom, "total": chromosomes[chrom], "pass": pass_chromosomes[chrom]} for chrom in sorted(chromosomes)])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        print(f"SV summary error: {error}", file=__import__("sys").stderr)
        raise SystemExit(2)
