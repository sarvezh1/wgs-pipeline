#!/usr/bin/env python3
"""Intersect immutable caller-membership outputs with benchmark strata."""

import argparse
import csv
import sys

from stratification_utils import (NA, class_matches, classes_for, fmt_fraction,
                                  intersect_intervals, interval_bases, overlaps,
                                  read_bed, read_manifest, sha256)


def read_table(path):
    with open(path, "r") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key(row):
    return row["chrom"], int(row["pos"]), row["ref"], row["alt"]


def comparison_category(row, analysis, callers):
    present = [caller for caller in callers if row[caller] == "1"]
    if analysis == "germline":
        if len(present) == 2:
            return "shared"
        return present[0] + "_only"
    if len(present) == 3:
        return "all_three_shared"
    if len(present) == 2:
        return "_".join(present) + "_only"
    return present[0] + "_specific"


def main():
    parser = argparse.ArgumentParser()
    for name in ("manifest", "benchmark_regions", "membership", "discordant", "output_prefix", "sample", "analysis", "benchmark_name", "benchmark_version"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    parser.add_argument("--genotype-details")
    args = parser.parse_args()
    if args.analysis not in ("germline", "somatic"):
        raise ValueError("unsupported analysis mode '%s'" % args.analysis)
    if args.analysis == "somatic" and args.genotype_details:
        raise ValueError("germline-only genotype stratification path cannot be used for somatic analysis")
    rows = read_table(args.membership)
    discordant = read_table(args.discordant)
    benchmark = read_bed(args.benchmark_regions, label="benchmark region BED")
    strata = []
    resources = read_manifest(args.manifest)
    for name, bed, _configured in resources:
        regions = intersect_intervals(read_bed(bed, label="stratum '%s' BED" % name), benchmark)
        bases = interval_bases(regions)
        if not bases:
            raise ValueError("stratum '%s' has no benchmark-region overlap" % name)
        strata.append((name, regions, bases))
    callers = ["haplotypecaller", "deepvariant"] if args.analysis == "germline" else ["mutect2", "strelka2", "deepsomatic"]
    pairs = [(callers[i], callers[j]) for i in range(len(callers)) for j in range(i + 1, len(callers))]
    categories = (["union", "shared", "haplotypecaller_only", "deepvariant_only"] if args.analysis == "germline" else
                  ["union", "all_three_shared", "mutect2_strelka2_only", "mutect2_deepsomatic_only", "strelka2_deepsomatic_only",
                   "mutect2_specific", "strelka2_specific", "deepsomatic_specific"])
    all_keys = [key(row) for row in rows]
    classes = classes_for(all_keys)

    with open(args.output_prefix + ".stratified_concordance_summary.tsv", "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "benchmark_name", "benchmark_version", "stratum", "variant_class",
                      "comparison_group", "count", "union_count", "fraction_of_union", "truth_match_count", "truth_match_fraction", "benchmark_region_bases"])
        for name, regions, bases in strata:
            selected = [row for row in rows if overlaps(regions, *key(row)[:3])]
            for klass in classes:
                class_rows = [row for row in selected if class_matches(klass, row["ref"], row["alt"])]
                for group_name in categories:
                    group = class_rows if group_name == "union" else [row for row in class_rows if comparison_category(row, args.analysis, callers) == group_name]
                    truth_matches = sum(row["truth_allele_match"] == "1" for row in group)
                    out.writerow([args.sample, args.analysis, args.benchmark_name, args.benchmark_version, name, klass, group_name,
                                  len(group), len(class_rows), fmt_fraction(len(group), len(class_rows)), truth_matches,
                                  fmt_fraction(truth_matches, len(group)), bases])

    with open(args.output_prefix + ".stratified_jaccard.tsv", "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "stratum", "caller_a", "caller_b", "variant_class", "intersection_count", "union_count", "jaccard", "benchmark_region_bases"])
        for name, regions, bases in strata:
            selected = [row for row in rows if overlaps(regions, *key(row)[:3])]
            for caller_a, caller_b in pairs:
                for klass in classes:
                    class_rows = [row for row in selected if class_matches(klass, row["ref"], row["alt"])]
                    union = [row for row in class_rows if row[caller_a] == "1" or row[caller_b] == "1"]
                    intersection = [row for row in union if row[caller_a] == "1" and row[caller_b] == "1"]
                    out.writerow([args.sample, args.analysis, name, caller_a, caller_b, klass, len(intersection), len(union),
                                  fmt_fraction(len(intersection), len(union)), bases])

    with open(args.output_prefix + ".stratified_discordant_variants.tsv", "w", newline="") as handle:
        if discordant:
            fields = list(discordant[0]) + ["strata"]
        else:
            fields = ["sample_or_patient", "analysis_type", "chrom", "pos", "ref", "alt", "variant_class", "strata"]
        out = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        out.writeheader()
        for row in discordant:
            names = [name for name, regions, _bases in strata if overlaps(regions, *key(row)[:3])]
            if names:
                row["strata"] = ";".join(names)
                out.writerow(row)

    if args.analysis == "germline":
        if not args.genotype_details:
            raise ValueError("germline stratification requires genotype-details")
        genotypes = read_table(args.genotype_details)
        with open(args.output_prefix + ".stratified_genotype_concordance.tsv", "w", newline="") as handle:
            out = csv.writer(handle, delimiter="\t", lineterminator="\n")
            out.writerow(["sample_or_patient", "stratum", "variant_class", "shared_alleles", "genotype_concordant", "genotype_discordant",
                          "missing_or_nocall", "callable_shared", "genotype_concordance_rate", "benchmark_region_bases"])
            for name, regions, bases in strata:
                selected = [row for row in genotypes if overlaps(regions, *key(row)[:3])]
                for klass in classes:
                    group = [row for row in selected if class_matches(klass, row["ref"], row["alt"])]
                    concordant = sum(row["comparison"] == "concordant" for row in group)
                    genotype_discordant = sum(row["comparison"] == "discordant" for row in group)
                    missing = sum(row["comparison"] == "missing_or_nocall" for row in group)
                    callable_count = concordant + genotype_discordant
                    out.writerow([args.sample, name, klass, len(group), concordant, genotype_discordant, missing, callable_count,
                                  fmt_fraction(concordant, callable_count), bases])

    bases_by_name = {name: bases for name, _regions, bases in strata}
    with open(args.output_prefix + ".stratified_comparison_provenance.tsv", "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["sample_or_patient", "analysis_type", "benchmark_name", "benchmark_version", "stratum", "configured_bed_path",
                      "staged_bed_path", "bed_checksum", "effective_benchmark_bases", "overlap_semantics", "stratum_relationship_semantics"])
        for name, bed, configured in resources:
            out.writerow([args.sample, args.analysis, args.benchmark_name, args.benchmark_version, name, configured, bed, sha256(bed),
                          bases_by_name[name], "merged BED union intersect benchmark-region union; variant REF span overlap",
                          "strata evaluated independently and may overlap; counts are not additive partitions"])


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.stderr.write("STRATIFIED_COMPARISON ERROR: %s\n" % error)
        sys.exit(1)
