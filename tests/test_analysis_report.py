#!/usr/bin/env python3

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "bin" / "generate_analysis_report.py"


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class AnalysisReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="block8-report-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.results = self.root / "results"
        self.results.mkdir()
        self.output = self.results / "report" / "case.somatic.wgs_report.html"
        self.sheet = self.root / "samples.tsv"
        write_tsv(self.sheet, ["sample_id", "patient_id", "role"], [
            {"sample_id": "normal", "patient_id": "case", "role": "normal"},
            {"sample_id": "tumor", "patient_id": "case", "role": "tumor"},
        ])

    def command(self, *extra, analysis="somatic", caller="compare", resource="fixture"):
        return [
            "python3", str(GENERATOR), "--input-root", str(self.results), "--output", str(self.output),
            "--samplesheet", str(self.sheet), "--analysis-type", analysis, "--caller", caller,
            "--resource-mode", resource, "--assembly", "GRCh38", "--generated-at", "2026-08-23T00:00:00+00:00",
            *extra,
        ]

    def run_report(self, *extra, **kwargs):
        result = subprocess.run(self.command(*extra, **kwargs), text=True, capture_output=True)
        document = self.output.read_text(encoding="utf-8") if self.output.exists() else ""
        return result, document

    def add_variant(self, caller="mutect2", sample="case"):
        write_tsv(self.results / "vcf_processing" / f"{sample}.{caller}.variant_summary.tsv",
                  ["caller", "analysis_type", "sample_or_patient", "total_records", "pass_records", "snvs", "indels", "multiallelic_records"],
                  [{"caller": caller, "analysis_type": "somatic", "sample_or_patient": sample, "total_records": "17", "pass_records": "14", "snvs": "15", "indels": "2", "multiallelic_records": "0"}])

    def test_missing_optional_sections_and_cnv_not_run(self):
        self.add_variant()
        result, document = self.run_report()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Not run", document)
        self.assertIn("CNV segments", document)
        self.assertNotIn("file://", document)

    def test_empty_sv_callset_is_valid_zero(self):
        write_tsv(self.results / "structural_variants/summaries/case.manta.sv_summary.tsv",
                  ["sample_or_patient", "analysis_type", "output_semantics", "engine", "engine_version", "reference_assembly", "total_sv", "pass_sv"],
                  [{"sample_or_patient": "case", "analysis_type": "somatic", "output_semantics": "somaticSV", "engine": "Manta", "engine_version": "1.6.0", "reference_assembly": "GRCh38", "total_sv": "0", "pass_sv": "0"}])
        result, document = self.run_report("--run-sv")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(">0<", document)
        self.assertIn("No non-empty SV class plot generated", document)

    def test_single_caller_somatic_mode(self):
        self.add_variant()
        result, document = self.run_report(caller="mutect2")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("mutect2", document)
        self.assertIn("Not run", document)

    def test_germline_mode_and_missing_benchmark_metrics(self):
        self.sheet = self.root / "germline.tsv"
        write_tsv(self.sheet, ["sample_id", "patient_id", "role"], [{"sample_id": "HG002", "patient_id": "HG002", "role": "germline"}])
        result, document = self.run_report("--run-benchmark", analysis="germline", caller="deepvariant")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Germline WGS Analysis Report", document)
        self.assertIn("Not available", document)

    def test_compare_mode_and_na_values(self):
        for caller in ("mutect2", "strelka2", "deepsomatic"):
            self.add_variant(caller)
        write_tsv(self.results / "comparison/case.pairwise_jaccard.tsv",
                  ["sample_or_patient", "analysis_type", "caller_a", "caller_b", "variant_class", "intersection_count", "union_count", "jaccard"],
                  [{"sample_or_patient": "case", "analysis_type": "somatic", "caller_a": "mutect2", "caller_b": "strelka2", "variant_class": "INDEL", "intersection_count": "0", "union_count": "0", "jaccard": "NA"}])
        result, document = self.run_report()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("mutect2 / strelka2", document)
        self.assertIn(">NA<", document)

    def test_malformed_summary_fails_actionably(self):
        write_tsv(self.results / "bad.variant_summary.tsv", ["caller", "total_records"], [{"caller": "x", "total_records": "1"}])
        result, _ = self.run_report()
        self.assertEqual(result.returncode, 2)
        self.assertIn("missing column(s)", result.stderr)

    def test_html_escaping_and_path_removal(self):
        write_tsv(self.results / "vcf_processing/unsafe.variant_summary.tsv",
                  ["caller", "analysis_type", "sample_or_patient", "total_records", "pass_records", "snvs", "indels", "multiallelic_records"],
                  [{"caller": "caller<script>alert(1)</script>", "analysis_type": "somatic", "sample_or_patient": "/mnt/d/private/user/sample", "total_records": "1", "pass_records": "1", "snvs": "1", "indels": "0", "multiallelic_records": "0"}])
        result, document = self.run_report()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("<script>alert", document)
        self.assertIn("&lt;script&gt;", document)
        self.assertNotIn("/mnt/d", document)

    def test_unsafe_local_link_is_rejected(self):
        result, _ = self.run_report("--multiqc-link", "file:///tmp/report.html")
        self.assertEqual(result.returncode, 2)
        self.assertIn("Unsafe local link", result.stderr)

    def test_relative_local_links_are_preserved(self):
        result, document = self.run_report("--multiqc-link", "../multiqc/somatic.multiqc.html")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('href="../multiqc/somatic.multiqc.html"', document)

    def test_deterministic_rendering(self):
        self.add_variant()
        first, document_a = self.run_report()
        self.assertEqual(first.returncode, 0, first.stderr)
        second, document_b = self.run_report()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(document_a, document_b)

    def test_male_aware_benchmark_provenance_is_exposed(self):
        self.sheet = self.root / "germline.tsv"
        write_tsv(self.sheet, ["sample_id", "patient_id", "role"],
                  [{"sample_id": "HG002", "patient_id": "HG002", "role": "germline"}])
        write_tsv(self.results / "benchmark/HG002.deepvariant.benchmark_provenance.tsv",
                  ["sample_or_patient", "analysis_type", "caller", "benchmark_sex_requested",
                   "benchmark_sex_resolved", "sex_resolution_source", "hap_py_gender_argument",
                   "sex_aware_benchmarking"],
                  [{"sample_or_patient": "HG002", "analysis_type": "germline", "caller": "deepvariant",
                    "benchmark_sex_requested": "auto", "benchmark_sex_resolved": "male",
                    "sex_resolution_source": "GIAB HG002 v5.0q benchmark identity and resource manifest",
                    "hap_py_gender_argument": "male", "sex_aware_benchmarking": "true"}])
        result, document = self.run_report("--run-benchmark", analysis="germline", caller="deepvariant", resource="production")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Benchmark ploidy handling", document)
        self.assertIn("Male-aware (hap.py)", document)


if __name__ == "__main__":
    unittest.main()
