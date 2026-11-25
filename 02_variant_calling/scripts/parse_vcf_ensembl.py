#!/usr/bin/env python3
"""
parse_vcf_ensembl.py
Parse snpEff-annotated VCF files within a target genomic interval, export per-consequence
variant tables, and optionally query the Ensembl REST API to resolve canonical / MANE Select
RefSeq transcripts and UniProt-SwissProt protein accessions.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ENSEMBL_REST = "https://rest.ensembl.org"
ROW_FIELDS = ["chrom", "pos", "rsID", "ref", "alt", "variant_type", "impact", "gene", "transcript", "nuc_change", "aa_change"]


def extract_ann_records(info_field: str) -> List[List[str]]:
    """Extract all snpEff ANN= pipe-delimited transcript records from a VCF INFO column."""
    ann_value = None
    for entry in info_field.split(";"):
        if entry.startswith("ANN="):
            ann_value = entry[4:]
            break
    if not ann_value:
        return []
    records: List[List[str]] = []
    for raw_rec in ann_value.split(","):
        fields = raw_rec.split("|")
        if len(fields) >= 11:
            records.append(fields)
    return records


def parse_vcf_interval(
    vcf_path: Path,
    start: int,
    end: int,
    out_vcf: Optional[Path] = None,
) -> Tuple[List[Dict[str, str]], List[str]]:
    """Filter variants in [start, end] (1-based inclusive) and unpack snpEff ANN annotations."""
    rows: List[Dict[str, str]] = []
    rsids: List[str] = []
    header_lines: List[str] = []
    kept_variant_lines: List[str] = []

    with vcf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_str, rsid, ref, alt, _, _, info = parts[:8]
            pos = int(pos_str)
            if pos < start or pos > end:
                continue
            kept_variant_lines.append(line)
            if rsid != "." and rsid not in rsids:
                rsids.append(rsid)
            for fields in extract_ann_records(info):
                rows.append({
                    "chrom": chrom,
                    "pos": pos_str,
                    "rsID": rsid,
                    "ref": ref,
                    "alt": alt,
                    "variant_type": fields[1],
                    "impact": fields[2],
                    "gene": fields[3],
                    "transcript": fields[6],
                    "nuc_change": fields[9],
                    "aa_change": fields[10],
                })

    if out_vcf is not None:
        out_vcf.parent.mkdir(parents=True, exist_ok=True)
        with out_vcf.open("w", encoding="utf-8") as out:
            if header_lines:
                out.writelines(header_lines)
            else:
                out.write("##fileformat=VCFv4.2\n")
                out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1\n")
            out.writelines(kept_variant_lines)

    return rows, rsids


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _ensembl_get_json(endpoint: str, params: Optional[Dict[str, str]] = None, timeout: int = 20):
    url = f"{ENSEMBL_REST}{endpoint}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "computational-omics-parser/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_gene_transcripts(gene_symbol: str, delay_s: float = 0.08) -> List[Dict[str, str]]:
    """Resolve Ensembl transcripts, MANE Select RefSeq IDs, and UniProt-SwissProt IDs via Ensembl REST."""
    xrefs = _ensembl_get_json(f"/xrefs/symbol/homo_sapiens/{urllib.parse.quote(gene_symbol)}", {"object_type": "gene"})
    gene_hits = [x for x in xrefs if x.get("type") == "gene"]
    if not gene_hits:
        raise ValueError(f"No human Ensembl gene found for symbol '{gene_symbol}'.")
    gene_id = gene_hits[0]["id"]

    expanded = _ensembl_get_json(f"/lookup/id/{gene_id}", {"expand": "1"})
    transcripts = expanded.get("Transcript", [])
    records: List[Dict[str, str]] = []

    for tx in transcripts:
        tx_id = tx["id"]
        is_canonical = "1" if tx.get("is_canonical") else ""
        mane_select = ""

        try:
            refseq_hits = _ensembl_get_json(f"/xrefs/id/{tx_id}", {"external_db": "RefSeq_mRNA"})
            for hit in refseq_hits:
                if "MANE_Select" in str(hit.get("info_text", "")) or "MANE Select" in str(hit.get("description", "")):
                    mane_select = str(hit.get("primary_id", ""))
                    break
            if not mane_select and len(refseq_hits) == 1 and is_canonical:
                mane_select = str(refseq_hits[0].get("primary_id", ""))
        except urllib.error.URLError:
            pass

        uniprot_id = ""
        translation = tx.get("Translation") or {}
        ensp = translation.get("id")
        if ensp:
            try:
                up_hits = _ensembl_get_json(f"/xrefs/id/{ensp}", {"external_db": "Uniprot/SWISSPROT"})
                if up_hits:
                    uniprot_id = str(up_hits[0].get("primary_id", ""))
            except urllib.error.URLError:
                pass

        records.append({
            "ensembl_gene_id": gene_id,
            "ensembl_transcript_id": tx_id,
            "transcript_is_canonical": is_canonical,
            "transcript_mane_select": mane_select,
            "uniprotswissprot": uniprot_id,
        })
        time.sleep(delay_s)

    return records


def filter_by_transcript(rows: List[Dict[str, str]], transcript_id: str) -> List[Dict[str, str]]:
    """Match transcript accession ignoring optional version suffix (e.g. NM_000546.6 == NM_000546)."""
    target_base = transcript_id.split(".")[0]
    return [r for r in rows if r["transcript"].split(".")[0] == target_base]


def main() -> None:
    ap = argparse.ArgumentParser(description="Parse snpEff-annotated VCFs and optionally query Ensembl REST transcripts.")
    ap.add_argument("--vcf", required=True, help="Input snpEff-annotated VCF file.")
    ap.add_argument("--start", type=int, required=True, help="1-based genomic interval start.")
    ap.add_argument("--end", type=int, required=True, help="1-based genomic interval end.")
    ap.add_argument("--outdir", default="outputs/variants", help="Output directory for parsed VCF and TSV tables.")
    ap.add_argument("--transcript", default=None, help="Optional transcript ID (e.g. NM_000546.6 or ENST...) for canonical subset.")
    ap.add_argument("--gene", default=None, help="HGNC gene symbol for Ensembl REST transcript resolution.")
    ap.add_argument("--query-ensembl", action="store_true", help="Query Ensembl REST API for MANE Select / UniProt IDs.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    type_dir = outdir / "variant_type_tables"
    rows, rsids = parse_vcf_interval(Path(args.vcf), args.start, args.end, out_vcf=outdir / "parsed_interval.vcf")

    write_tsv(type_dir / "all_variant_types.tsv", rows, ROW_FIELDS)
    by_type: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        by_type.setdefault(row["variant_type"] or "unknown", []).append(row)
    for vtype, subset in sorted(by_type.items()):
        safe_name = "".join(ch if (ch.isalnum() or ch in ("_", "-")) else "_" for ch in vtype)
        write_tsv(type_dir / f"{safe_name}.tsv", subset, ROW_FIELDS)

    canonical_tx = args.transcript
    if args.query_ensembl:
        if not args.gene:
            ap.error("--gene is required when --query-ensembl is set.")
        tx_rows = fetch_gene_transcripts(args.gene)
        write_tsv(
            outdir / f"{args.gene.lower()}_ensembl_transcripts.tsv",
            tx_rows,
            ["ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical", "transcript_mane_select", "uniprotswissprot"],
        )
        if not canonical_tx:
            for tr in tx_rows:
                if tr["transcript_is_canonical"] == "1":
                    canonical_tx = tr["transcript_mane_select"] or tr["ensembl_transcript_id"]
                    break

    if canonical_tx:
        canon_rows = filter_by_transcript(rows, canonical_tx)
        write_tsv(outdir / "canonical_variants.tsv", canon_rows, ROW_FIELDS)
        print(f"[OK] Wrote {len(canon_rows)} canonical records for transcript {canonical_tx}")

    print(f"[OK] Parsed {len(rows)} transcript annotations across {len(rsids)} unique rsID(s) -> {outdir}")


if __name__ == "__main__":
    main()
