"""Unit tests for snpEff ANN VCF parsing, interval slicing, transcript filtering, and CLI execution."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "parse_vcf_ensembl.py"
spec = importlib.util.spec_from_file_location("parse_vcf_ensembl", SCRIPT_PATH)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


SAMPLE_VCF = """##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1
17\t7661779\t.\tG\tA\t300\tPASS\tDP=40;ANN=A|3_prime_UTR_variant|MODIFIER|TP53|TP53|transcript|NM_000546.6|protein_coding|11/11|c.*10C>T||||||\tGT\t0/1
17\t7673802\trs28934576\tC\tT\t999\tPASS\tDP=120;ANN=T|missense_variant|MODERATE|TP53|TP53|transcript|NM_000546.6|protein_coding|8/11|c.839G>A|p.Arg280Lys|1000/2500|839/1182|280/393||,T|synonymous_variant|LOW|TP53|TP53|transcript|NM_001126112.3|protein_coding|8/11|c.839G>A|p.Arg280Arg|1000/2500|839/1182|280/393||;FS=0.0\tGT\t0/1
17\t7687550\trs28934576\tA\tG\t420\tPASS\tDP=60;ANN=G|5_prime_UTR_variant|MODIFIER|TP53|TP53|transcript|ENST00000269305.9|protein_coding|1/11|c.-5T>C||||||\tGT\t0/1
17\t7690000\trs123456\tA\tG\t500\tPASS\tDP=80;ANN=G|intron_variant|MODIFIER|TP53|TP53|transcript|NM_000546.6|protein_coding|1/10|c.10+5T>C||||||\tGT\t0/1
"""


class TestParseVcfEnsembl(unittest.TestCase):
    def test_interval_boundaries_and_ann_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            vcf_file = td_path / "sample.vcf"
            vcf_file.write_text(SAMPLE_VCF, encoding="utf-8")
            out_vcf = td_path / "filtered.vcf"

            rows, rsids = mod.parse_vcf_interval(vcf_file, 7661779, 7687550, out_vcf=out_vcf)
            # '.' excluded, duplicate rsID deduplicated, out-of-range rs123456 excluded
            self.assertEqual(rsids, ["rs28934576"])
            # Boundary variants at 7661779 (1 ANN) + 7673802 (2 ANN) + 7687550 (1 ANN) = 4 rows
            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[1]["variant_type"], "missense_variant")
            self.assertEqual(rows[1]["transcript"], "NM_000546.6")
            self.assertEqual(rows[1]["aa_change"], "p.Arg280Lys")
            self.assertEqual(rows[2]["variant_type"], "synonymous_variant")
            self.assertTrue(out_vcf.exists())
            out_text = out_vcf.read_text(encoding="utf-8")
            self.assertIn("7661779", out_text)
            self.assertIn("7687550", out_text)
            self.assertNotIn("rs123456", out_text)

    def test_extract_ann_records_missing_and_malformed(self) -> None:
        self.assertEqual(mod.extract_ann_records("DP=50;MQ=60"), [])
        self.assertEqual(mod.extract_ann_records("DP=50;ANN=;MQ=60"), [])
        self.assertEqual(mod.extract_ann_records("ANN=short|field|only"), [])
        recs = mod.extract_ann_records("DP=10;ANN=T|missense_variant|MODERATE|TP53|TP53|transcript|NM_000546.6|protein_coding|8/11|c.839G>A|p.Arg280Lys,bad|entry;FS=1.2")
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0][10], "p.Arg280Lys")

    def test_filter_by_transcript_version_agnostic(self) -> None:
        rows = [
            {"transcript": "NM_000546.6", "aa_change": "p.Arg280Lys"},
            {"transcript": "ENST00000269305.9", "aa_change": "p.Arg280Lys"},
            {"transcript": "NM_001126112.3", "aa_change": "p.Arg280Arg"},
        ]
        self.assertEqual(len(mod.filter_by_transcript(rows, "NM_000546")), 1)
        self.assertEqual(len(mod.filter_by_transcript(rows, "NM_000546.2")), 1)
        self.assertEqual(len(mod.filter_by_transcript(rows, "ENST00000269305")), 1)
        self.assertEqual(mod.filter_by_transcript(rows, "NM_999999"), [])

    def test_end_to_end_cli_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            vcf_file = td_path / "in.vcf"
            outdir = td_path / "out"
            vcf_file.write_text(SAMPLE_VCF, encoding="utf-8")

            res = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--vcf", str(vcf_file),
                    "--start", "7661779",
                    "--end", "7687550",
                    "--outdir", str(outdir),
                    "--transcript", "NM_000546.6",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("Parsed 4 transcript annotations", res.stdout)
            self.assertTrue((outdir / "parsed_interval.vcf").is_file())
            self.assertTrue((outdir / "variant_type_tables" / "all_variant_types.tsv").is_file())
            self.assertTrue((outdir / "variant_type_tables" / "missense_variant.tsv").is_file())

            with (outdir / "canonical_variants.tsv").open("r", encoding="utf-8") as fh:
                canon_rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(canon_rows), 2)
            self.assertEqual({r["variant_type"] for r in canon_rows}, {"3_prime_UTR_variant", "missense_variant"})


if __name__ == "__main__":
    unittest.main()
