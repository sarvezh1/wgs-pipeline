#!/usr/bin/env python3
"""Create standardized stratified metrics for normalized exact-allele benchmarks."""

import argparse
import csv
import sys

from stratification_utils import (NA, class_matches, classes_for, fmt_fraction,
                                  intersect_intervals, interval_bases, metric_values,
                                  overlaps, read_allele_tsv, read_bed, read_manifest,
                                  sha256, variant_class)


def main():
    parser = argparse.ArgumentParser()
    for name in ("manifest", "benchmark_regions", "truth_alleles", "query_alleles", "output", "errors", "provenance",
                 "sample", "analysis", "caller", "benchmark_name", "benchmark_version", "engine", "engine_version", "assembly"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    args = parser.parse_args()
    truth, query = read_allele_tsv(args.truth_alleles), read_allele_tsv(args.query_alleles)
    benchmark = read_bed(args.benchmark_regions, label="benchmark region BED")
    strata = []
    for name, bed, configured in read_manifest(args.manifest):
        effective = intersect_intervals(read_bed(bed, label="stratum '%s' BED" % name), benchmark)
        bases = interval_bases(effective)
        if not bases:
            raise ValueError("stratum '%s' has no benchmark-region overlap" % name)
        strata.append((name, bed, configured, effective, bases))

    header = ["sample_or_patient", "analysis_type", "caller", "benchmark_name", "benchmark_version", "stratum",
              "variant_class", "truth_tp", "query_tp", "fp", "fn", "precision", "recall", "f1",
              "benchmark_region_bases", "engine", "engine_version", "reference_assembly"]
    with open(args.output, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(header)
        all_keys = truth | query
        for stratum, _bed, _configured, regions, bases in strata:
            in_truth = {key for key in truth if overlaps(regions, *key[:3])}
            in_query = {key for key in query if overlaps(regions, *key[:3])}
            for wanted in classes_for(all_keys):
                klass = "OVERALL" if wanted == "ALL" else wanted
                truth_keys = {key for key in in_truth if class_matches(wanted, key[2], key[3])}
                query_keys = {key for key in in_query if class_matches(wanted, key[2], key[3])}
                shared = truth_keys & query_keys
                truth_tp = query_tp = len(shared)
                fp, fn = len(query_keys - truth_keys), len(truth_keys - query_keys)
                precision, recall, f1 = metric_values(truth_tp, query_tp, fp, fn)
                out.writerow([args.sample, args.analysis, args.caller, args.benchmark_name, args.benchmark_version,
                              stratum, klass, truth_tp, query_tp, fp, fn, precision, recall, f1, bases,
                              args.engine, args.engine_version, args.assembly])

    with open(args.errors, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "caller", "chrom", "pos", "ref", "alt",
                      "variant_class", "benchmark_classification", "strata"])
        for key in sorted(truth ^ query, key=lambda x: (x[0], x[1], x[2], x[3])):
            names = [name for name, _bed, _configured, regions, _bases in strata if overlaps(regions, *key[:3])]
            if names:
                out.writerow([args.sample, args.analysis, args.caller, *key, variant_class(key[2], key[3]),
                              "FN" if key in truth else "FP", ";".join(names)])

    with open(args.provenance, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "caller", "stratum", "configured_bed_path", "staged_bed_path",
                      "bed_checksum", "reference_assembly", "benchmark_name", "benchmark_version", "engine", "engine_version",
                      "effective_benchmark_bases", "overlap_semantics", "within_stratum_overlap_semantics"])
        for name, bed, configured, _regions, bases in strata:
            out.writerow([args.sample, args.analysis, args.caller, name, configured, bed, sha256(bed), args.assembly,
                          args.benchmark_name, args.benchmark_version, args.engine, args.engine_version, bases,
                          "BED union intersect benchmark-region BED; variant REF span overlaps effective interval",
                          "merged union; overlapping input rows are not double-counted"])


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.stderr.write("STRATIFIED_BENCHMARK ERROR: %s\n" % error)
        sys.exit(1)
