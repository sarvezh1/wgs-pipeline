#!/usr/bin/env python3
"""Dependency-free helpers shared by Block 5 stratification scripts."""

import gzip
import hashlib
import os


NA = "NA"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def read_manifest(path):
    strata = []
    seen = set()
    base = os.path.dirname(os.path.abspath(path))
    with open(path, "r") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2 or not all(fields):
                raise ValueError("stratification manifest line %d must be stratum_name<TAB>bed_path" % number)
            name, configured = fields
            if name in seen:
                raise ValueError("duplicate stratum name '%s'" % name)
            seen.add(name)
            bed = configured if os.path.isabs(configured) else os.path.join(base, configured)
            strata.append((name, os.path.abspath(bed), configured))
    if not strata:
        raise ValueError("stratification manifest contains no resource rows")
    return strata


def read_bed(path, contig_lengths=None, label="BED"):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError("missing or empty %s: '%s'" % (label, path))
    intervals = {}
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError("malformed %s line %d: expected at least three tab-delimited columns" % (label, number))
            chrom, start_text, end_text = fields[:3]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                raise ValueError("malformed %s line %d: coordinates must be integers" % (label, number))
            if start < 0 or end <= start:
                raise ValueError("malformed %s line %d: require 0 <= start < end" % (label, number))
            if contig_lengths is not None:
                if chrom not in contig_lengths:
                    raise ValueError("%s line %d uses unknown reference contig '%s'" % (label, number, chrom))
                if end > contig_lengths[chrom]:
                    raise ValueError("%s line %d ends beyond reference contig '%s' (%d > %d)" %
                                     (label, number, chrom, end, contig_lengths[chrom]))
            intervals.setdefault(chrom, []).append((start, end))
    if not intervals:
        raise ValueError("%s contains no positive-width intervals" % label)
    return merge_intervals(intervals)


def merge_intervals(intervals):
    merged = {}
    for chrom, values in intervals.items():
        result = []
        for start, end in sorted(values):
            if not result or start > result[-1][1]:
                result.append([start, end])
            elif end > result[-1][1]:
                result[-1][1] = end
        merged[chrom] = [(start, end) for start, end in result]
    return merged


def intersect_intervals(left, right):
    output = {}
    for chrom in set(left).intersection(right):
        a, b, i, j = left[chrom], right[chrom], 0, 0
        values = []
        while i < len(a) and j < len(b):
            start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if start < end:
                values.append((start, end))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
        if values:
            output[chrom] = values
    return merge_intervals(output)


def interval_bases(intervals):
    return sum(end - start for values in intervals.values() for start, end in values)


def overlaps(intervals, chrom, pos, ref):
    start, end = int(pos) - 1, int(pos) - 1 + max(1, len(ref))
    return any(start < right and end > left for left, right in intervals.get(chrom, ()))


def variant_class(ref, alt):
    if not ref or not alt or alt == "*" or alt.startswith("<") or "[" in alt or "]" in alt:
        return "OTHER"
    if not all(base.upper() in "ACGTN" for base in ref + alt):
        return "OTHER"
    return "SNV" if len(ref) == len(alt) == 1 else "INDEL"


def indel_size_class(ref, alt):
    if variant_class(ref, alt) != "INDEL":
        return None
    size = abs(len(ref) - len(alt))
    if size == 0:
        return None
    if size == 1:
        return "INDEL_1BP"
    if size <= 5:
        return "INDEL_2_5BP"
    if size <= 15:
        return "INDEL_6_15BP"
    return "INDEL_GT15BP"


def class_matches(wanted, ref, alt):
    observed = variant_class(ref, alt)
    if wanted == "ALL":
        return observed in ("SNV", "INDEL")
    if wanted in ("SNV", "INDEL"):
        return wanted == observed
    return wanted == indel_size_class(ref, alt)


def classes_for(keys):
    result = ["ALL", "SNV", "INDEL"]
    present = {indel_size_class(ref, alt) for _chrom, _pos, ref, alt in keys}
    result.extend(name for name in ("INDEL_1BP", "INDEL_2_5BP", "INDEL_6_15BP", "INDEL_GT15BP") if name in present)
    return result


def fmt_fraction(numerator, denominator):
    return NA if denominator == 0 else "%.12g" % (float(numerator) / denominator)


def metric_values(truth_tp, query_tp, fp, fn):
    precision = fmt_fraction(query_tp, query_tp + fp)
    recall = fmt_fraction(truth_tp, truth_tp + fn)
    if precision == NA or recall == NA or float(precision) + float(recall) == 0:
        f1 = NA if precision == NA or recall == NA else "0"
    else:
        f1 = "%.12g" % (2 * float(precision) * float(recall) / (float(precision) + float(recall)))
    return precision, recall, f1


def read_allele_tsv(path):
    keys = set()
    with open(path, "r") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 4:
                raise ValueError("allele table '%s' line %d does not have four columns" % (path, number))
            chrom, pos, ref, alt = fields
            keys.add((chrom, int(pos), ref, alt))
    return keys


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()
