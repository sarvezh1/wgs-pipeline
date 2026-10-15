#!/usr/bin/env python
from __future__ import print_function

import csv
import sys


def fail(message):
    sys.stderr.write("HAPPY_SUMMARY ERROR: %s\n" % message)
    sys.exit(1)


if len(sys.argv) != 11:
    fail("expected SUMMARY_CSV OUTPUT_TSV SAMPLE ANALYSIS CALLER BENCHMARK_NAME BENCHMARK_VERSION BENCHMARK_MODE REGION_BASES ASSEMBLY")

summary_csv, output_tsv, sample, analysis, caller, benchmark_name, benchmark_version, benchmark_mode, region_bases, assembly = sys.argv[1:]

with open(summary_csv, "r") as handle:
    rows = list(csv.DictReader(handle))

required_columns = set(["Type", "Filter", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK"])
if not rows or not required_columns.issubset(set(rows[0].keys())):
    fail("hap.py summary is empty or missing required GA4GH count columns")


def choose(variant_type):
    candidates = [row for row in rows if row.get("Type") == variant_type]
    for field, expected in [("Filter", "ALL"), ("Subset", "*"), ("Genotype", "*")]:
        narrowed = [row for row in candidates if row.get(field) == expected]
        if narrowed:
            candidates = narrowed
    if len(candidates) != 1:
        fail("could not identify one aggregate %s row in hap.py summary" % variant_type)
    row = candidates[0]
    try:
        return {
            "truth_tp": int(float(row["TRUTH.TP"])),
            "query_tp": int(float(row["QUERY.TOTAL"])) - int(float(row["QUERY.FP"])) - int(float(row["QUERY.UNK"])),
            "fp": int(float(row["QUERY.FP"])),
            "fn": int(float(row["TRUTH.FN"])),
        }
    except (KeyError, TypeError, ValueError):
        fail("non-numeric GA4GH counts in aggregate %s row" % variant_type)


def metrics(counts):
    precision_denominator = counts["query_tp"] + counts["fp"]
    recall_denominator = counts["truth_tp"] + counts["fn"]
    precision = float(counts["query_tp"]) / precision_denominator if precision_denominator else 0.0
    recall = float(counts["truth_tp"]) / recall_denominator if recall_denominator else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


snv = choose("SNP")
indel = choose("INDEL")
overall = dict((key, snv[key] + indel[key]) for key in ["truth_tp", "query_tp", "fp", "fn"])

header = [
    "sample_or_patient", "analysis_type", "caller", "benchmark_name", "benchmark_version",
    "benchmark_mode", "variant_class", "truth_tp", "query_tp", "fp", "fn",
    "precision", "recall", "f1", "benchmark_region_bases", "engine", "engine_version",
    "reference_assembly"
]

with open(output_tsv, "w") as handle:
    handle.write("\t".join(header) + "\n")
    for variant_class, counts in [("OVERALL", overall), ("SNV", snv), ("INDEL", indel)]:
        precision, recall, f1 = metrics(counts)
        values = [
            sample, analysis, caller, benchmark_name, benchmark_version, benchmark_mode, variant_class,
            str(counts["truth_tp"]), str(counts["query_tp"]), str(counts["fp"]), str(counts["fn"]),
            "%.12g" % precision, "%.12g" % recall, "%.12g" % f1, str(region_bases),
            "hap.py/vcfeval", "hap.py 0.3.15; RTG Tools 3.12.1", assembly
        ]
        handle.write("\t".join(values) + "\n")
