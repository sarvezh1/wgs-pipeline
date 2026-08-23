import csv
import os
import subprocess
import tempfile
import unittest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(ROOT, "resources", "GRCh38.production-resource-manifest.template.tsv")
VALIDATOR = os.path.join(ROOT, "bin", "validate_resource_config.sh")


class ProductionResourceContractTests(unittest.TestCase):
    def test_manifest_template_schema_and_frozen_identities(self):
        with open(TEMPLATE, newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(rows)
        self.assertEqual(
            list(rows[0]),
            ["resource_name", "resource_type", "assembly", "version", "source",
             "path_or_filename", "checksum", "required_for"],
        )
        self.assertTrue(all(row["assembly"] == "GRCh38" for row in rows))
        self.assertTrue(all(not row["checksum"] or row["checksum"].startswith("sha256:") for row in rows))
        by_name = {row["resource_name"]: row for row in rows}
        self.assertIn("GIAB HG002 v5.0q", by_name["hg002_primary_truth"]["version"])
        self.assertIn("GIAB/GA4GH v3.6", by_name["giab_stratifications_tsv"]["version"])
        self.assertIn("Ensembl VEP 116", by_name["vep_cache"]["version"])
        self.assertEqual(by_name["mutect2_panel_of_normals"]["version"], "study/platform-specific")

    def validator_command(self, directory, include_alt):
        reference = os.path.join(directory, "GRCh38.fa")
        fai = reference + ".fai"
        dictionary = os.path.join(directory, "GRCh38.dict")
        prefix = os.path.join(directory, "GRCh38")
        manifest = os.path.join(directory, "manifest.tsv")
        with open(reference, "w") as handle:
            handle.write(">chr1\nAAAAAAAAAA\n>chr1_test_alt\nAAAAA\n")
        with open(fai, "w") as handle:
            handle.write("chr1\t10\t0\t10\t11\nchr1_test_alt\t5\t0\t5\t6\n")
        with open(dictionary, "w") as handle:
            handle.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr1_test_alt\tLN:5\n")
        for suffix in (".0123", ".amb", ".bwt.2bit.64", ".pac"):
            with open(prefix + suffix, "w") as handle:
                handle.write("index\n")
        with open(prefix + ".ann", "w") as handle:
            handle.write("15 2 0\n")
        if include_alt:
            with open(prefix + ".alt", "w") as handle:
                handle.write("chr1_test_alt\n")
        with open(manifest, "w") as handle:
            handle.write("resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n")
            handle.write("reference\treference_fasta\tGRCh38\tfixture\ttest\tGRCh38.fa\t\tall\n")
        return [
            VALIDATOR, "production", "GRCh38", ROOT, reference, fai, dictionary, prefix,
            "", manifest, "false", "false", "cache", "", "", "false", "", "germline",
            "haplotypecaller", "vcf", "basic", "", "", "",
        ]

    def test_alt_reference_requires_and_accepts_bwa_alt_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = subprocess.run(
                self.validator_command(directory, include_alt=False), cwd=directory,
                text=True, capture_output=True,
            )
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("ALT metadata is missing", missing.stderr)
            accepted = subprocess.run(
                self.validator_command(directory, include_alt=True), cwd=directory,
                text=True, capture_output=True,
            )
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            with open(os.path.join(directory, "resource_manifest.tsv")) as handle:
                self.assertIn("bwa_mem2_alt", handle.read())


if __name__ == "__main__":
    unittest.main()
