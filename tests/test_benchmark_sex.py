#!/usr/bin/env python3
import csv
import gzip
import os
import subprocess
import tempfile
import unittest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VALIDATOR = os.path.join(ROOT, "bin", "validate_benchmark_config.sh")
REFERENCE = os.path.join(ROOT, "test_data/reference/GRCh38_chr20.fa.gz")
FIXTURE_TRUTH = os.path.join(ROOT, "test_data/reference/benchmarks/HG002_GRCh38_giab_v4.2.1_chr20_12446968_12455543.bgzip.vcf.gz")
FIXTURE_BED = os.path.join(ROOT, "test_data/reference/germline_test_interval.bed")
FIXTURE_MANIFEST = os.path.join(ROOT, "test_data/reference/resource_manifest.tsv")


class BenchmarkSexTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="benchmark-sex-")
        self.addCleanup(self.temp.cleanup)
        self.directory = self.temp.name

    def command(self, sex="auto", *, mode="germline_happy", analysis="germline",
                resource="fixture", name="HG002_GIAB", version="v4.2.1",
                truth=FIXTURE_TRUTH, bed=FIXTURE_BED, manifest=FIXTURE_MANIFEST):
        return [VALIDATOR, "true", mode, resource, analysis, "compare", name, version, sex,
                "GRCh38", "GRCh38", ROOT, REFERENCE, "", truth, truth + ".tbi", bed, "",
                manifest, "engine"]

    def run_validator(self, *args, **kwargs):
        return subprocess.run(self.command(*args, **kwargs), cwd=self.directory,
                              text=True, capture_output=True)

    def make_primary(self):
        truth = os.path.join(self.directory, "HG002_GRCh38_v5.0q_smvar.vcf.gz")
        bed = os.path.join(self.directory, "HG002_GRCh38_v5.0q_smvar.benchmark.bed")
        with gzip.open(truth, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr20\t12446969\t.\tA\tC\t.\tPASS\t.\n")
        with open(truth + ".tbi", "w") as handle:
            handle.write("synthetic plumbing index\n")
        with open(bed, "w") as handle:
            handle.write("chr20\t12446967\t12446970\n")
        manifest = os.path.join(self.directory, "manifest.tsv")
        with open(manifest, "w") as handle:
            handle.write("resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n")
            handle.write(f"truth\ttruth_vcf\tGRCh38\tGIAB HG002 v5.0q\ttest\t{truth}\t\tbenchmarking\n")
            handle.write(f"index\ttabix_index\tGRCh38\tGIAB HG002 v5.0q\ttest\t{truth}.tbi\t\tbenchmarking\n")
            handle.write(f"regions\tbenchmark_regions_bed\tGRCh38\tGIAB HG002 v5.0q\ttest\t{bed}\t\tbenchmarking\n")
        return truth, bed, manifest

    def run_primary(self, sex):
        truth, bed, manifest = self.make_primary()
        return self.run_validator(sex, resource="production", version="v5.0q",
                                  truth=truth, bed=bed, manifest=manifest)

    def read_config(self):
        with open(os.path.join(self.directory, "benchmark_config_manifest.tsv"), newline="") as handle:
            return next(csv.DictReader(handle, delimiter="\t"))

    def test_primary_auto_resolves_male_with_provenance(self):
        result = self.run_primary("auto")
        self.assertEqual(result.returncode, 0, result.stderr)
        row = self.read_config()
        self.assertEqual(row["benchmark_sex_requested"], "auto")
        self.assertEqual(row["benchmark_sex_resolved"], "male")
        self.assertEqual(row["hap_py_gender_argument"], "male")
        self.assertEqual(row["sex_aware_benchmarking"], "true")

    def test_primary_explicit_male(self):
        result = self.run_primary("male")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.read_config()["sex_resolution_source"], "explicit_declaration")

    def test_primary_none_and_female_fail(self):
        for sex, message in (("none", "requires male-aware"), ("female", "conflicts with trusted male")):
            with self.subTest(sex=sex):
                result = self.run_primary(sex)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)

    def test_invalid_value_fails_before_resource_checks(self):
        result = self.run_validator("unknown")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unsupported benchmark_sex", result.stderr)

    def test_somatic_auto_is_unaffected_and_male_is_rejected(self):
        auto = self.run_validator("auto", mode="somatic", analysis="somatic", name="fixture", version="v1")
        self.assertEqual(auto.returncode, 0, auto.stderr)
        self.assertEqual(self.read_config()["sex_aware_benchmarking"], "false")
        male = self.run_validator("male", mode="somatic", analysis="somatic", name="fixture", version="v1")
        self.assertNotEqual(male.returncode, 0)
        self.assertIn("only to germline_happy", male.stderr)

    def test_legacy_auto_preserves_non_sex_aware_semantics(self):
        result = self.run_validator("auto")
        self.assertEqual(result.returncode, 0, result.stderr)
        row = self.read_config()
        self.assertEqual(row["benchmark_sex_resolved"], "none")
        self.assertEqual(row["hap_py_gender_argument"], "not_passed")
        self.assertEqual(row["sex_aware_benchmarking"], "false")

    def test_v5_v421_truth_region_mixing_is_rejected(self):
        truth, _, manifest = self.make_primary()
        legacy_bed = os.path.join(self.directory, "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
        with open(legacy_bed, "w") as handle:
            handle.write("chr20\t12446967\t12446970\n")
        mixed = self.run_validator("male", resource="production", version="v5.0q",
                                   truth=truth, bed=legacy_bed, manifest=manifest)
        self.assertNotEqual(mixed.returncode, 0)
        self.assertIn("v5.0q truth cannot be combined with v4.2.1", mixed.stderr)

    def test_v421_truth_cannot_be_declared_v5(self):
        result = self.run_validator("male", version="v5.0q")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("v4.2.1 resources cannot be represented", result.stderr)

    def test_happy_command_routes_validated_gender_without_changing_parser(self):
        with open(os.path.join(ROOT, "modules/BENCHMARK.nf"), encoding="utf-8") as handle:
            module = handle.read()
        self.assertIn('gender_arg="--gender \\$hap_py_gender"', module)
        self.assertIn("parse_happy_stratified.py", module)


if __name__ == "__main__":
    unittest.main()
