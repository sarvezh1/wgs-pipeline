import csv
import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class WgsModalityTests(unittest.TestCase):
    def run_script(self, script, *args):
        return subprocess.run(["python3", str(ROOT / "bin" / script), *map(str, args)], text=True, capture_output=True)

    def test_valid_configuration_and_unknown_ploidy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sheet = root / "samples.tsv"
            sheet.write_text("sample_id\tpatient_id\trole\tfastq_r1\tfastq_r2\nS\tP\tsample\ta\tb\n")
            manifest = root / "manifest.tsv"
            manifest.write_text("resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n")
            provenance, ploidy = root / "provenance.tsv", root / "ploidy.tsv"
            result = self.run_script(
                "validate_wgs_modalities.py", "--samplesheet", sheet, "--analysis", "germline", "--assembly", "GRCh38",
                "--resource-mode", "fixture", "--resource-manifest", manifest, "--project-dir", root,
                "--run-sv", "true", "--run-cnv", "false", "--run-coverage", "true",
                "--coverage-thresholds", "10,20,30", "--provenance-out", provenance, "--ploidy-out", ploidy,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("unknown\t2\tNA\tNA\tNA", ploidy.read_text())

    def test_invalid_thresholds_and_fixture_cnv_fail_early(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sheet = root / "samples.tsv"
            sheet.write_text("sample_id\tpatient_id\trole\tfastq_r1\tfastq_r2\nS\tP\tsample\ta\tb\n")
            manifest = root / "manifest.tsv"
            manifest.write_text("resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n")
            common = ["--samplesheet", sheet, "--analysis", "germline", "--assembly", "GRCh38", "--resource-mode", "fixture", "--resource-manifest", manifest, "--project-dir", root, "--run-sv", "false", "--run-coverage", "false", "--provenance-out", root / "p.tsv", "--ploidy-out", root / "x.tsv"]
            bad = self.run_script("validate_wgs_modalities.py", *common, "--run-cnv", "false", "--coverage-thresholds", "10,-2")
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("positive integers", bad.stderr)
            guarded = self.run_script("validate_wgs_modalities.py", *common, "--run-cnv", "true", "--coverage-thresholds", "10")
            self.assertNotEqual(guarded.returncode, 0)
            self.assertIn("CNV requires production mode", guarded.stderr)

    def test_invalid_boolean_and_ploidy_declaration(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sheet = root / "samples.tsv"
            sheet.write_text("sample_id\tpatient_id\trole\tfastq_r1\tfastq_r2\nS\tP\tsample\ta\tb\n")
            manifest = root / "manifest.tsv"
            manifest.write_text("resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n")
            ploidy = root / "ploidy.tsv"
            ploidy.write_text("sample_id\tsex_chromosome_configuration\nS\tXO\n")
            args = ["--samplesheet", sheet, "--analysis", "germline", "--assembly", "GRCh38", "--resource-mode", "fixture", "--resource-manifest", manifest, "--project-dir", root, "--run-cnv", "false", "--run-coverage", "false", "--coverage-thresholds", "10", "--provenance-out", root / "p.tsv", "--ploidy-out", root / "x.tsv"]
            invalid_boolean = self.run_script("validate_wgs_modalities.py", *args, "--run-sv", "yes")
            self.assertNotEqual(invalid_boolean.returncode, 0)
            self.assertIn("run_sv must be true or false", invalid_boolean.stderr)
            invalid_ploidy = self.run_script("validate_wgs_modalities.py", *args, "--run-sv", "false", "--ploidy-config", ploidy)
            self.assertNotEqual(invalid_ploidy.returncode, 0)
            self.assertIn("unsupported configuration", invalid_ploidy.stderr)

    def test_production_cnv_missing_resources_is_actionable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sheet = root / "samples.tsv"
            sheet.write_text("sample_id\tpatient_id\trole\tfastq_r1\tfastq_r2\nS\tP\tsample\ta\tb\n")
            manifest = root / "manifest.tsv"
            manifest.write_text("resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n")
            result = self.run_script(
                "validate_wgs_modalities.py", "--samplesheet", sheet, "--analysis", "germline", "--assembly", "GRCh38",
                "--resource-mode", "production", "--resource-manifest", manifest, "--project-dir", root,
                "--run-sv", "false", "--run-cnv", "true", "--run-coverage", "false", "--coverage-thresholds", "10",
                "--cnv-intervals", root / "missing.interval_list", "--provenance-out", root / "p.tsv", "--ploidy-out", root / "x.tsv",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("CNV interval list is missing", result.stderr)

    def test_sv_summary_valid_empty_and_unknown_contig(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fai = root / "ref.fa.fai"
            fai.write_text("chr20\t1000\t0\t0\t0\n")
            vcf = root / "calls.vcf.gz"
            with gzip.open(vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            result = self.run_script("summarize_sv.py", "--vcf", vcf, "--fai", fai, "--sample", "S", "--patient", "P", "--analysis", "germline", "--semantics", "diploidSV", "--assembly", "GRCh38", "--prefix", root / "S")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("\t0\t0\n", (root / "S.sv_summary.tsv").read_text())
            with gzip.open(vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrX\t2\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=10\n")
            bad = self.run_script("summarize_sv.py", "--vcf", vcf, "--fai", fai, "--sample", "S", "--patient", "P", "--analysis", "germline", "--semantics", "diploidSV", "--assembly", "GRCh38", "--prefix", root / "S")
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("absent from the reference", bad.stderr)

    def test_coverage_summary_and_impossible_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = root / "S.mosdepth.summary.txt"
            summary.write_text("chrom\tlength\tbases\tmean\tmin\tmax\nchr20\t100\t100\t5.00\t0\t10\ntotal\t100\t100\t5.00\t0\t10\n")
            regions = root / "S.regions.bed.gz"
            threshold = root / "S.thresholds.bed.gz"
            with gzip.open(regions, "wt") as handle:
                handle.write("chr20\t0\t100\tregion\t5.00\n")
            with gzip.open(threshold, "wt") as handle:
                handle.write("#chrom\tstart\tend\tregion\t1\t10\nchr20\t0\t100\tregion\t90\t20\n")
            result = self.run_script("summarize_mosdepth.py", "--summary", summary, "--regions", regions, "--thresholds", threshold, "--sample", "S", "--patient", "P", "--role", "sample", "--analysis", "germline", "--assembly", "GRCh38", "--prefix", root / "S")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("\t10\t20\t100\t0.2\t", (root / "S.coverage_threshold_summary.tsv").read_text())
            with gzip.open(threshold, "wt") as handle:
                handle.write("#chrom\tstart\tend\tregion\t1\nchr20\t0\t100\tregion\t101\n")
            bad = self.run_script("summarize_mosdepth.py", "--summary", summary, "--regions", regions, "--thresholds", threshold, "--sample", "S", "--patient", "P", "--role", "sample", "--analysis", "germline", "--assembly", "GRCh38", "--prefix", root / "S")
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("impossible threshold count", bad.stderr)

    def test_cnv_native_classifications(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            seg = root / "called.seg"
            seg.write_text("CONTIG\tSTART\tEND\tMEAN_LOG2_COPY_RATIO\tCALL\nchr20\t1\t100\t0.4\t+\nchr20\t101\t150\t-0.5\t-\nchr20\t151\t200\t0.0\t0\n")
            out = root / "summary.tsv"
            result = self.run_script("summarize_cnv.py", "--input", seg, "--format", "somatic_seg", "--id", "P", "--analysis", "somatic", "--assembly", "GRCh38", "--output", out)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("\t3\t1\t1\t1\t150\t100\t", out.read_text())


if __name__ == "__main__":
    unittest.main()
