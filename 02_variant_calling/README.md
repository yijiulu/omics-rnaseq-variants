# Human Illumina Variant Calling & Ensembl Transcript Annotation

Unix-composable variant calling and functional annotation workflow combining **BWA-MEM**, **GATK 4**, **snpEff**, **bcftools** (dbSNP `rsID` tagging), and a Python parser with optional **Ensembl REST API** cross-referencing for MANE Select and UniProt/SwissProt accessions.

## Pipeline Stages

1. **`scripts/run_variant_pipeline.sh`**:
   - Aligns paired-end Illumina reads to GRCh38 with `bwa mem` and read-group headers.
   - Marks duplicates (`GATK MarkDuplicatesSpark`), runs 2-pass `HaplotypeCaller` with GATK hard-filtering (`QD`, `FS`, `MQ`, `SOR`, `MQRankSum`, `ReadPosRankSum`) and Base Quality Score Recalibration (`BaseRecalibrator` + `ApplyBQSR`).
   - Annotates final filtered SNPs with `snpEff` (`GRCh38.p14`).
   - Redirects all tool logs to `stderr` and writes only the final annotated VCF path to `stdout` for clean piping.
2. **`scripts/annotate_against_known.sh`**:
   - Reads the upstream VCF path from `stdin`, compresses and indexes with `bgzip`/`tabix`, annotates known dbSNP `rsID`s via `bcftools annotate -c ID`, and optionally slices variants within `[gene_start, gene_end]`.
3. **`scripts/parse_vcf_ensembl.py`**:
   - Properly parses the `ANN=` INFO subfield from snpEff-annotated VCFs (preserving other INFO keys), slices genomic coordinates, exports per-consequence TSV tables (`missense_variant.tsv`, etc.), and optionally queries the Ensembl REST API (rate-limited `< 15 req/s`) to resolve canonical/MANE Select RefSeq transcripts and UniProt-SwissProt protein IDs.

## Usage

```bash
# Run alignment + GATK calling + snpEff piped directly into dbSNP annotation
SNPEFF_JAR=/path/to/snpEff.jar \
  bash scripts/run_variant_pipeline.sh brca1 17 \
    BRCA1_paired_dat1.fq BRCA1_paired_dat2.fq Homo_sapiens.GRCh38.dna.chromosome.17.fa | \
  bash scripts/annotate_against_known.sh /path/to/dbsnp_chr17.vcf.gz 43044295 43170245

# Parse snpEff ANN annotations and filter to canonical/MANE transcript
python3 scripts/parse_vcf_ensembl.py \
  --vcf brca1_filtered_snps_final_ann_dbsnp.vcf \
  --start 43044295 \
  --end 43170245 \
  --outdir outputs/brca1 \
  --gene BRCA1 \
  --query-ensembl
```

### HPC vs Local/Conda Execution

Set `USE_MODULES=1` (default `0`) when running on an HPC cluster with Environment Modules/Lmod, and override tool executables via `BWA_CMD`, `GATK_CMD`, `GATK_SIF`, `SAMTOOLS_CMD`, `JAVA_CMD`, `SNPEFF_JAR`, `BGZIP_CMD`, `TABIX_CMD`, and `BCFTOOLS_CMD`.

## Testing

```bash
python3 -m unittest discover -s tests -v
```
