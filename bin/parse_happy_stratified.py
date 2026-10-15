#!/usr/bin/env python3
"""Extract native hap.py subset rows and record-level GA4GH decisions."""

import argparse
import csv
import gzip
import sys

from stratification_utils import (NA, intersect_intervals, interval_bases, metric_values,
                                  overlaps, read_bed, read_manifest, sha256, variant_class)


def number(row, field):
    try:
        return int(float(row[field]))
    except (KeyError, TypeError, ValueError):
        raise ValueError("hap.py stratified row has non-numeric '%s'" % field)


def select_row(rows, subset, variant_type):
    candidates = [row for row in rows if row.get("Subset") == subset and row.get("Type") == variant_type]
    for field, expected in (("Filter", "ALL"), ("Subtype", "*"), ("Genotype", "*"), ("QQ", "*")):
        narrowed = [row for row in candidates if row.get(field) == expected]
        if narrowed:
            candidates = narrowed
    if not candidates:
        return None
    if len(candidates) != 1:
        raise ValueError("could not identify one native hap.py %s/%s aggregate row" % (subset, variant_type))
    row = candidates[0]
    return {"truth_tp": number(row, "TRUTH.TP"), "query_tp": number(row, "QUERY.TOTAL") - number(row, "QUERY.FP") - number(row, "QUERY.UNK"),
            "fp": number(row, "QUERY.FP"), "fn": number(row, "TRUTH.FN")}


def parse_format(format_text, sample_text):
    names, values = format_text.split(":"), sample_text.split(":")
    return dict(zip(names, values))


def main():
    parser = argparse.ArgumentParser()
    for name in ("extended", "happy_vcf", "manifest", "benchmark_regions", "output", "errors", "provenance", "sample", "caller",
                 "benchmark_name", "benchmark_version", "engine_version", "assembly"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    args = parser.parse_args()
    with open(args.extended, "r") as handle:
        rows = list(csv.DictReader(handle))
    required = {"Type", "Subset", "Filter", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("hap.py extended CSV is empty or missing native stratified count columns")
    benchmark = read_bed(args.benchmark_regions, label="benchmark region BED")
    strata = []
    for name, bed, configured in read_manifest(args.manifest):
        regions = intersect_intervals(read_bed(bed, label="stratum '%s' BED" % name), benchmark)
        bases = interval_bases(regions)
        if not bases:
            raise ValueError("stratum '%s' has no benchmark-region overlap" % name)
        strata.append((name, bed, configured, regions, bases))

    header = ["sample_or_patient", "analysis_type", "caller", "benchmark_name", "benchmark_version", "stratum", "variant_class",
              "truth_tp", "query_tp", "fp", "fn", "precision", "recall", "f1", "benchmark_region_bases", "engine", "engine_version", "reference_assembly"]
    with open(args.output, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(header)
        for name, _bed, _configured, _regions, bases in strata:
            native = {klass: select_row(rows, name, happy_type) for klass, happy_type in (("SNV", "SNP"), ("INDEL", "INDEL"))}
            for klass in ("OVERALL", "SNV", "INDEL"):
                if klass == "OVERALL":
                    if native["SNV"] is None or native["INDEL"] is None:
                        continue
                    counts = {key: native["SNV"][key] + native["INDEL"][key] for key in native["SNV"]}
                else:
                    counts = native[klass]
                    if counts is None:
                        continue
                precision, recall, f1 = metric_values(counts["truth_tp"], counts["query_tp"], counts["fp"], counts["fn"])
                out.writerow([args.sample, "germline", args.caller, args.benchmark_name, args.benchmark_version, name, klass,
                              counts["truth_tp"], counts["query_tp"], counts["fp"], counts["fn"], precision, recall, f1, bases,
                              "hap.py/vcfeval", args.engine_version, args.assembly])

    # The annotated VCF is the authoritative record-level classification source. Its
    # row count need not equal summary FN counts because hap.py counts haplotypes.
    with gzip.open(args.happy_vcf, "rt") as source, open(args.errors, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "caller", "chrom", "pos", "ref", "alt", "variant_class",
                      "benchmark_classification", "strata"])
        samples = []
        for line in source:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 11 or len(samples) < 2:
                continue
            truth = parse_format(fields[8], fields[9])
            query = parse_format(fields[8], fields[10])
            decisions = []
            if truth.get("BD") == "FN": decisions.append("FN")
            if query.get("BD") == "FP": decisions.append("FP")
            if not decisions:
                continue
            chrom, pos, ref = fields[0], int(fields[1]), fields[3]
            for alt in fields[4].split(","):
                names = [name for name, _bed, _configured, regions, _bases in strata if overlaps(regions, chrom, pos, ref)]
                if names:
                    out.writerow([args.sample, "germline", args.caller, chrom, pos, ref, alt, variant_class(ref, alt),
                                  ";".join(decisions), ";".join(names)])

    with open(args.provenance, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "caller", "stratum", "configured_bed_path", "staged_bed_path", "bed_checksum",
                      "reference_assembly", "benchmark_name", "benchmark_version", "engine", "engine_version", "effective_benchmark_bases",
                      "overlap_semantics", "within_stratum_overlap_semantics", "error_classification_semantics"])
        for name, bed, configured, _regions, bases in strata:
            out.writerow([args.sample, "germline", args.caller, name, configured, bed, sha256(bed), args.assembly, args.benchmark_name,
                          args.benchmark_version, "hap.py/vcfeval", args.engine_version, bases,
                          "BED union intersect benchmark-region BED; hap.py native Subset counts", "merged union; input rows not double-counted",
                          "record-level GA4GH BD annotations; not assumed to reproduce haplotype-level summary counts"])


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.stderr.write("HAPPY_STRATIFIED ERROR: %s\n" % error)
        sys.exit(1)
