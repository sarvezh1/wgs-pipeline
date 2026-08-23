#!/usr/bin/env python3

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARSER = ROOT / "bin" / "parse_nextflow_trace.py"
FIELDS = [
    "task_id", "hash", "native_id", "name", "tag", "status", "exit", "submit", "duration",
    "realtime", "%cpu", "peak_rss", "peak_vmem", "rchar", "wchar", "read_bytes",
    "write_bytes", "memory", "cpus", "container", "attempt", "workdir",
]


def record(**overrides):
    row = {
        "task_id": "1", "hash": "ab/cdef12", "native_id": "42",
        "name": "WGS_PIPELINE:DEEPVARIANT (HG002:deepvariant)",
        "tag": "HG002:deepvariant", "status": "COMPLETED", "exit": "0",
        "submit": "2026-01-01 00:00:00.000", "duration": "12s", "realtime": "10s",
        "%cpu": "150%", "peak_rss": "1 GB", "peak_vmem": "2 GB", "rchar": "1 KB",
        "wchar": "2 KB", "read_bytes": "512 MB", "write_bytes": "1 GB",
        "memory": "2.5 GB", "cpus": "1", "container": "google/deepvariant:1.8.0",
        "attempt": "1", "workdir": "/private/work/ab/cdef12",
    }
    row.update(overrides)
    return row


class PerformanceParserTests(unittest.TestCase):
    def run_parser(self, rows, *, fields=FIELDS, analysis="germline", caller="compare"):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        trace = root / "trace.txt"
        with trace.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        output = root / "performance"
        command = [
            "python3", str(PARSER), "--trace", str(trace), "--output-dir", str(output),
            "--analysis", analysis, "--caller", caller, "--resource-mode", "fixture",
            "--bqsr", "true", "--benchmark", "true", "--comparison", "true",
            "--stratification", "true", "--workflow-start", "start",
            "--workflow-completion", "complete", "--workflow-elapsed-seconds", "20",
            "--workflow-success", "true",
        ]
        return subprocess.run(command, text=True, capture_output=True), output

    def test_outputs_and_cached_tasks_are_separated(self):
        completed = record()
        cached = record(task_id="2", hash="ef/123456", status="CACHED", realtime="-", **{"%cpu": "-", "peak_rss": "-", "peak_vmem": "-", "read_bytes": "-", "write_bytes": "-"})
        result, output = self.run_parser([completed, cached])
        self.assertEqual(result.returncode, 0, result.stderr)
        with (output / "performance_summary.tsv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["workdir_id"], "ab/cdef12")
        self.assertNotIn("/private/work", "\n".join(rows[0].values()))
        self.assertEqual(rows[0]["metric_scope"], "TECHNICAL_VALIDATION_METRICS")
        with (output / "caller_performance.tsv").open(newline="", encoding="utf-8") as handle:
            caller_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(caller_rows[0]["executed_tasks"], "1")
        self.assertEqual(caller_rows[0]["cached_tasks"], "1")
        self.assertEqual(caller_rows[0]["summed_task_wall_seconds"], "10")

    def test_caller_table_excludes_downstream_benchmark_tasks(self):
        benchmark = record(
            task_id="3",
            hash="aa/bb2233",
            name="HAPPY_BENCHMARK (HG002:deepvariant:hap.py)",
            tag="HG002:deepvariant:hap.py",
            realtime="25s",
        )
        result, output = self.run_parser([record(), benchmark])
        self.assertEqual(result.returncode, 0, result.stderr)
        with (output / "caller_performance.tsv").open(newline="", encoding="utf-8") as handle:
            caller_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(caller_rows[0]["executed_tasks"], "1")
        self.assertEqual(caller_rows[0]["summed_task_wall_seconds"], "10")

    def test_block7_tasks_have_distinct_stages(self):
        manta = record(name="MANTA_GERMLINE (HG002:manta)", tag="HG002:manta")
        mosdepth = record(task_id="2", hash="cd/ef1234", name="MOSDEPTH (HG002:coverage)", tag="HG002:coverage")
        cnv = record(task_id="3", hash="ef/ab1234", name="GERMLINE_CNV_CALLER_CASE (HG002:gCNV CASE)", tag="HG002:gCNV CASE")
        result, output = self.run_parser([manta, mosdepth, cnv])
        self.assertEqual(result.returncode, 0, result.stderr)
        with (output / "performance_summary.tsv").open(newline="", encoding="utf-8") as handle:
            stages = {row["process"]: row["stage"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(stages["MANTA_GERMLINE"], "structural variation")
        self.assertEqual(stages["MOSDEPTH"], "coverage QC")
        self.assertEqual(stages["GERMLINE_CNV_CALLER_CASE"], "copy-number analysis")

    def test_missing_required_column_fails_actionably(self):
        fields = [field for field in FIELDS if field != "peak_rss"]
        row = {key: value for key, value in record().items() if key in fields}
        result, _ = self.run_parser([row], fields=fields)
        self.assertEqual(result.returncode, 2)
        self.assertIn("missing required column(s): peak_rss", result.stderr)

    def test_negative_metric_fails(self):
        result, _ = self.run_parser([record(realtime="-2s")])
        self.assertEqual(result.returncode, 2)
        self.assertIn("negative realtime", result.stderr)

    def test_malformed_metric_fails(self):
        result, _ = self.run_parser([record(peak_rss="large")])
        self.assertEqual(result.returncode, 2)
        self.assertIn("Invalid peak RSS", result.stderr)

    def test_incompatible_analysis_caller_metadata_fails(self):
        result, _ = self.run_parser([record()], analysis="somatic", caller="haplotypecaller")
        self.assertEqual(result.returncode, 2)
        self.assertIn("Unsupported caller metadata", result.stderr)

    def test_negative_workflow_elapsed_time_fails(self):
        result, _ = self.run_parser([record()])
        self.assertEqual(result.returncode, 0, result.stderr)
        # Exercise argparse-compatible invocation separately to avoid duplicating the fixture writer.
        command = [
            "python3", str(PARSER), "--trace", str(output := Path(tempfile.gettempdir()) / "missing-trace"),
            "--output-dir", str(output.parent), "--analysis", "germline", "--caller", "compare",
            "--resource-mode", "fixture", "--bqsr", "true", "--benchmark", "true",
            "--comparison", "true", "--stratification", "true", "--workflow-start", "start",
            "--workflow-completion", "complete", "--workflow-elapsed-seconds", "-1",
            "--workflow-success", "true",
        ]
        invalid = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(invalid.returncode, 2)
        self.assertIn("negative workflow elapsed time", invalid.stderr)


if __name__ == "__main__":
    unittest.main()
