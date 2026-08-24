import csv
import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "render_genome_circos.py"


class GenomeCircosVisualizationTests(unittest.TestCase):
    """Controlled visualization-only records; these make no biological claim."""

    def run_renderer(self, directory, *optional):
        root = Path(directory)
        command = [
            "python3", str(SCRIPT), "--vcf", str(root / "small.vcf"),
            "--fai", str(root / "reference.fa.fai"), "--sample", "TEST_ONLY",
            "--caller", "testcaller", "--analysis-type", "somatic",
            "--bin-size", "100", "--output", str(root / "plot.html"),
            "--track-summary", str(root / "tracks.tsv"),
            "--provenance", str(root / "provenance.tsv"), "--no-render", *map(str, optional),
        ]
        return subprocess.run(command, text=True, capture_output=True, check=False)

    @staticmethod
    def base_files(root):
        (root / "reference.fa.fai").write_text("chr2\t250\nchr1\t300\nchrX\t180\nchrUn_KI\t90\n")
        (root / "small.vcf").write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tA\tG\t.\tPASS\t.\n"
            "chr1\t100\t.\tA\tAT\t.\tPASS\t.\n"
            "chr1\t101\t.\tAC\tA\t.\tPASS\t.\n"
            "chr2\t250\t.\tC\tT\t.\tPASS\t.\n"
            "chrX\t20\t.\tA\tG\t.\tLowQual\t.\n"
            "chrUn_KI\t4\t.\tA\tG\t.\tPASS\t.\n"
            "chr1\t301\t.\tA\tG\t.\tPASS\t.\n"
        )

    def test_order_bins_optional_tracks_and_real_sv_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.base_files(root)
            (root / "tumor.regions.bed").write_text("chr1\t0\t100\twindow\t31.5\nchr2\t0\t250\twindow\t28\nchr1\t299\t301\tbad\t9\n")
            (root / "normal.regions.bed").write_text("chr1\t0\t100\twindow\t30\n")
            (root / "called.seg").write_text("CONTIG\tSTART\tEND\tCALL\nchr1\t20\t80\t+\nchr2\t40\t90\t-\nchrX\t1\t100\t0\n")
            (root / "manta.vcf").write_text(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t40\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=90\n"
                "chr1\t120\tbnd1\tN\tN]chr2:150]\t.\tPASS\tSVTYPE=BND;MATEID=bnd2\n"
                "chr2\t150\tbnd2\tN\t[chr1:120[N\t.\tPASS\tSVTYPE=BND;MATEID=bnd1\n"
            )
            result = self.run_renderer(
                root, "--coverage", f"tumor:T={root / 'tumor.regions.bed'}",
                "--coverage", f"normal:N={root / 'normal.regions.bed'}",
                "--cnv", root / "called.seg", "--cnv-format", "somatic_seg", "--sv", root / "manta.vcf",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (root / "tracks.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            density = [row for row in rows if row["track"] == "small_variant_density"]
            self.assertEqual([row["chrom"] for row in density[:6:2]], ["chr1", "chr1", "chr1"])
            self.assertEqual(next(row["value"] for row in density if row["chrom"] == "chr1" and row["start_0based"] == "0" and row["series"] == "PASS_SNV"), "1")
            self.assertEqual(next(row["value"] for row in density if row["chrom"] == "chr1" and row["start_0based"] == "0" and row["series"] == "PASS_INDEL"), "1")
            self.assertEqual(next(row["value"] for row in density if row["chrom"] == "chr1" and row["start_0based"] == "100" and row["series"] == "PASS_INDEL"), "1")
            links = [row for row in rows if row["track"] == "sv_link"]
            self.assertEqual([(row["series"], row["value"]) for row in links], [("DEL", "chr1:90"), ("BND", "chr2:150")])
            with (root / "provenance.tsv").open() as handle:
                provenance = dict(csv.reader(handle, delimiter="\t"))
            self.assertEqual(provenance["included_contigs"], "chr1,chr2,chrX")
            self.assertEqual(provenance["excluded_contigs"], "chrUn_KI")
            self.assertEqual(provenance["small_variant_invalid_coordinate_records"], "1")
            self.assertEqual(provenance["sv_pass_links"], "2")
            self.assertEqual(provenance["cnv_segments"], "3")
            self.assertEqual(provenance["coverage_records"], "3")

    def test_missing_optional_tracks_and_zero_sv_are_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.base_files(root)
            (root / "empty_manta.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            result = self.run_renderer(root, "--sv", root / "empty_manta.vcf")
            self.assertEqual(result.returncode, 0, result.stderr)
            with (root / "provenance.tsv").open() as handle:
                provenance = dict(csv.reader(handle, delimiter="\t"))
            self.assertEqual(provenance["coverage_sources"], "NOT_PROVIDED")
            self.assertEqual(provenance["cnv_source"], "NOT_PROVIDED")
            self.assertEqual(provenance["sv_pass_links"], "0")

    def test_no_canonical_chromosomes_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "reference.fa.fai").write_text("chrUn_KI\t90\n")
            (root / "small.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            result = self.run_renderer(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("no canonical human chromosomes", result.stderr)

    @unittest.skipIf(importlib.util.find_spec("plotly") is None, "Plotly runtime is supplied by the visualization container")
    def test_rendered_html_is_offline_and_portable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.base_files(root)
            command = [
                "python3", str(SCRIPT), "--vcf", str(root / "small.vcf"),
                "--fai", str(root / "reference.fa.fai"), "--sample", "TEST_ONLY",
                "--caller", "testcaller", "--analysis-type", "somatic", "--bin-size", "100",
                "--showcase",
                "--output", str(root / "plot.html"), "--track-summary", str(root / "tracks.tsv"),
                "--provenance", str(root / "provenance.tsv"),
            ]
            result = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            html = (root / "plot.html").read_text()
            self.assertIn("Plotly.newPlot", html)
            self.assertIn("Genome-mode renderer demonstration", html)
            self.assertIn("Data status: visualization demonstration", html)
            self.assertIn('"name":"SNV density"', html)
            self.assertIn('"name":"INDEL density"', html)
            for forbidden in ('<script src=', 'file://', '/tmp/', '/work/'):
                self.assertNotIn(forbidden, html)


if __name__ == "__main__":
    unittest.main()
