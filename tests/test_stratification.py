#!/usr/bin/env python3
import csv
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from stratification_utils import intersect_intervals, interval_bases, read_bed, read_manifest


class StratificationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="block5-tests-")
        self.directory = self.temp.name

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, text):
        path = os.path.join(self.directory, name)
        with open(path, "w") as handle:
            handle.write(text)
        return path

    def test_missing_manifest_and_bed(self):
        with self.assertRaises(OSError):
            read_manifest(os.path.join(self.directory, "missing.tsv"))
        with self.assertRaisesRegex(ValueError, "missing or empty"):
            read_bed(os.path.join(self.directory, "missing.bed"))

    def test_duplicate_name(self):
        manifest = self.write("strata.tsv", "same\ta.bed\nsame\tb.bed\n")
        with self.assertRaisesRegex(ValueError, "duplicate stratum"):
            read_manifest(manifest)

    def test_malformed_and_invalid_coordinates(self):
        for name, contents in (("columns.bed", "chr20\t1\n"), ("negative.bed", "chr20\t-1\t2\n"),
                               ("zero.bed", "chr20\t2\t2\n"), ("text.bed", "chr20\tx\t2\n")):
            with self.assertRaisesRegex(ValueError, "malformed"):
                read_bed(self.write(name, contents))

    def test_unknown_contig_and_out_of_reference(self):
        with self.assertRaisesRegex(ValueError, "unknown reference contig"):
            read_bed(self.write("unknown.bed", "chrX\t0\t2\n"), {"chr20": 10})
        with self.assertRaisesRegex(ValueError, "beyond reference"):
            read_bed(self.write("long.bed", "chr20\t0\t11\n"), {"chr20": 10})

    def test_merge_intersection_base_accounting(self):
        stratum = read_bed(self.write("overlap.bed", "chr20\t0\t10\nchr20\t5\t20\n"))
        benchmark = read_bed(self.write("benchmark.bed", "chr20\t3\t8\nchr20\t7\t15\n"))
        self.assertEqual(interval_bases(intersect_intervals(stratum, benchmark)), 12)
        disjoint = read_bed(self.write("disjoint.bed", "chr20\t30\t40\n"))
        self.assertEqual(interval_bases(intersect_intervals(disjoint, benchmark)), 0)

    def benchmark_validator(self, manifest, resource_manifest=None, analysis="germline", mode="regression"):
        resource_manifest = resource_manifest or os.path.join(ROOT, "test_data/reference/resource_manifest.tsv")
        truth = os.path.join(ROOT, "test_data/reference/benchmarks/HG002_GRCh38_giab_v4.2.1_chr20_12446968_12455543.bgzip.vcf.gz")
        command = [os.path.join(ROOT, "bin/validate_benchmark_config.sh"), "true", mode, "fixture", analysis, "compare",
                   "fixture", "v1", "auto", "GRCh38", "GRCh38", ROOT, os.path.join(ROOT, "test_data/reference/GRCh38_chr20.fa.gz"), "",
                   truth, truth + ".tbi", os.path.join(ROOT, "test_data/reference/germline_test_interval.bed"), manifest,
                   resource_manifest, "engine"]
        return subprocess.run(command, cwd=self.directory, text=True, capture_output=True)

    def test_no_benchmark_overlap_is_actionable(self):
        self.write("outside.bed", "chr20\t20000000\t20000010\n")
        result = self.benchmark_validator(self.write("strata.tsv", "outside\toutside.bed\n"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no non-empty overlap with benchmark regions", result.stderr)

    def test_incompatible_stratum_assembly_is_actionable(self):
        bed = self.write("inside.bed", "chr20\t12446967\t12446970\n")
        manifest = self.write("strata.tsv", "inside\tinside.bed\n")
        resources = self.write("resources.tsv", "resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n"
                               "inside\tstratification_bed\tGRCh37\tv1\ttest\t%s\t\ttest\n" % bed)
        result = self.benchmark_validator(manifest, resources)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("declares assembly 'GRCh37'", result.stderr)

    def test_unsupported_mode_and_germline_path_on_somatic(self):
        script = os.path.join(ROOT, "bin/stratify_comparison.py")
        common = [script, "--manifest", "x", "--benchmark-regions", "x", "--membership", "x", "--discordant", "x",
                  "--output-prefix", "x", "--sample", "x", "--benchmark-name", "x", "--benchmark-version", "x"]
        unsupported = subprocess.run(common + ["--analysis", "rna"], text=True, capture_output=True)
        self.assertIn("unsupported analysis mode", unsupported.stderr)
        somatic = subprocess.run(common + ["--analysis", "somatic", "--genotype-details", "x"], text=True, capture_output=True)
        self.assertIn("germline-only genotype stratification", somatic.stderr)

    def test_malformed_resource_manifest(self):
        malformed = self.write("resources.tsv", "wrong\theader\n")
        command = [os.path.join(ROOT, "bin/validate_resource_config.sh"), "fixture", "GRCh38", ROOT,
                   os.path.join(ROOT, "test_data/reference/GRCh38_chr20.fa.gz"), "", "", "", "", malformed,
                   "false", "false", "cache", "", "", "false", "", "germline", "compare", "vcf", "basic", "", "", ""]
        result = subprocess.run(command, cwd=self.directory, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("malformed manifest header", result.stderr)

    def test_empty_metric_denominators_are_na(self):
        bed = self.write("inside.bed", "chr20\t12446967\t12446970\n")
        manifest = self.write("strata.tsv", "empty\tinside.bed\n")
        benchmark = self.write("benchmark.bed", "chr20\t12446967\t12446970\n")
        truth = self.write("truth.tsv", "")
        query = self.write("query.tsv", "")
        output, errors, provenance = [os.path.join(self.directory, name) for name in ("summary.tsv", "errors.tsv", "provenance.tsv")]
        command = [os.path.join(ROOT, "bin/stratify_exact_benchmark.py"), "--manifest", manifest, "--benchmark-regions", benchmark,
                   "--truth-alleles", truth, "--query-alleles", query, "--output", output, "--errors", errors, "--provenance", provenance,
                   "--sample", "S", "--analysis", "somatic", "--caller", "C", "--benchmark-name", "B", "--benchmark-version", "V",
                   "--engine", "E", "--engine-version", "EV", "--assembly", "GRCh38"]
        subprocess.run(command, check=True)
        with open(output) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["benchmark_region_bases"], "3")
        self.assertEqual((rows[0]["precision"], rows[0]["recall"], rows[0]["f1"]), ("NA", "NA", "NA"))


if __name__ == "__main__":
    unittest.main()
