#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./run_variant_pipeline.sh <gene> <chromosome> [read1.fq] [read2.fq] [reference.fa]

Examples:
  SNPEFF_JAR=/path/to/snpEff.jar ./run_variant_pipeline.sh vegfa 6
  SNPEFF_JAR=/path/to/snpEff.jar ./run_variant_pipeline.sh brca1 17 \
    BRCA1_paired_dat1.fq BRCA1_paired_dat2.fq Homo_sapiens.GRCh38.dna.chromosome.17.fa
  SNPEFF_JAR=/path/to/snpEff.jar bash run_variant_pipeline.sh brca1 17 | \
    bash annotate_against_known.sh /path/to/dbsnp_chr17.vcf.gz 43044295 43170245

Assumptions:
  - Default reference: Homo_sapiens.GRCh38.dna.chromosome.<chromosome>.fa
  - Default reads: <GENE_UPPER>_paired_dat1.fq and <GENE_UPPER>_paired_dat2.fq
  - Outputs <gene_lower>_filtered_snps_final_ann.vcf.
  - stdout emits only the final annotated VCF path so downstream tools can pipe it cleanly.

Environment overrides:
  USE_MODULES=0|1            Load HPC Environment Modules/Lmod if set to 1 (default: 0)
  SNPEFF_JAR=/path/snpEff.jar  Path to snpEff.jar (or set SNPEFF_CMD=snpEff)
  SNPEFF_GENOME=GRCh38.p14   snpEff genome build (default: GRCh38.p14)
  BWA_CMD=bwa
  GATK_CMD=gatk
  GATK_SIF=""                Optional Singularity image path for GATK fallback
  SAMTOOLS_CMD=samtools
  JAVA_CMD=java
USAGE
}

if [[ $# -lt 2 || $# -gt 5 ]]; then
  usage >&2
  exit 1
fi

GENE_INPUT=$1
CHR=$2

gene=$(printf '%s' "$GENE_INPUT" | tr '[:upper:]' '[:lower:]')
GENE_UPPER=$(printf '%s' "$GENE_INPUT" | tr '[:lower:]' '[:upper:]')

READ1=${3:-"${GENE_UPPER}_paired_dat1.fq"}
READ2=${4:-"${GENE_UPPER}_paired_dat2.fq"}
REF=${5:-"Homo_sapiens.GRCh38.dna.chromosome.${CHR}.fa"}

SAM="${gene}_aligned_reads.sam"
DEDUP_BAM="${gene}_sorted_dedup_reads.bam"
RAW_VCF="${gene}_raw_variants.vcf"
RAW_SNPS="${gene}_raw_snps.vcf"
FILTERED_SNPS="${gene}_filtered_snps.vcf"
BQSR_SNPS="${gene}_bqsr_snps.vcf"
RECAL_TABLE="${gene}_recal_data.table"
RECAL_BAM="${gene}_recal_reads.bam"
RAW_RECAL_VCF="${gene}_raw_variants_recal.vcf"
RAW_RECAL_SNPS="${gene}_raw_snps_recal.vcf"
FINAL_SNPS="${gene}_filtered_snps_final.vcf"
FINAL_ANN="${gene}_filtered_snps_final_ann.vcf"

USE_MODULES=${USE_MODULES:-0}
BWA_CMD=${BWA_CMD:-bwa}
GATK_CMD=${GATK_CMD:-gatk}
GATK_SIF=${GATK_SIF:-}
SAMTOOLS_CMD=${SAMTOOLS_CMD:-samtools}
JAVA_CMD=${JAVA_CMD:-java}
SNPEFF_CMD=${SNPEFF_CMD:-}
SNPEFF_JAR=${SNPEFF_JAR:-}
SNPEFF_GENOME=${SNPEFF_GENOME:-GRCh38.p14}
GATK_RUN=()

# Keep stdout clean for downstream pipes; send all logs to stderr.
exec 3>&1
exec 1>&2

require_file() {
  local file=$1
  if [[ ! -f "$file" ]]; then
    echo "Missing required file: $file" >&2
    exit 1
  fi
}

log_step() {
  printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"
}

maybe_load_module() {
  local mod=$1
  if [[ "$USE_MODULES" != "1" ]]; then
    return
  fi
  if ! type module >/dev/null 2>&1; then
    local init_script
    for init_script in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
      if [[ -r "$init_script" ]]; then
        # shellcheck disable=SC1090
        source "$init_script"
        break
      fi
    done
  fi
  if type module >/dev/null 2>&1; then
    module load "$mod"
  else
    echo "USE_MODULES=1 requested but 'module' command is unavailable." >&2
    exit 127
  fi
}

require_cmd() {
  local cmd=$1
  local cmd_type
  cmd_type=$(type -t "$cmd" 2>/dev/null || true)
  if [[ -z "$cmd_type" ]]; then
    echo "Required command not found in PATH: $cmd" >&2
    exit 127
  fi
  if [[ "$cmd_type" == "alias" ]]; then
    echo "Required command resolves to a shell alias (unsafe in scripts): $cmd" >&2
    exit 127
  fi
}

init_gatk_cmd() {
  local gatk_type
  gatk_type=$(type -t "$GATK_CMD" 2>/dev/null || true)
  if [[ "$gatk_type" == "file" || "$gatk_type" == "function" ]]; then
    GATK_RUN=("$GATK_CMD")
  elif [[ -n "$GATK_SIF" && -f "$GATK_SIF" ]] && type singularity >/dev/null 2>&1; then
    GATK_RUN=(singularity exec -e "$GATK_SIF" gatk)
  else
    echo "Required GATK executable not found: $GATK_CMD (or set GATK_SIF=/path/to/gatk.sif)" >&2
    exit 127
  fi
}

run_gatk() {
  "${GATK_RUN[@]}" "$@"
}

require_file "$REF"
require_file "$READ1"
require_file "$READ2"

log_step "Inputs"
echo "Gene: $gene | Chromosome: $CHR | Reference: $REF"
echo "Reads: $READ1, $READ2 | Output: $FINAL_ANN"

log_step "Index reference and align reads with BWA-MEM"
maybe_load_module "${BWA_MODULE:-bwa}"
require_cmd "$BWA_CMD"
"$BWA_CMD" index "$REF"
"$BWA_CMD" mem -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' \
  "$REF" "$READ1" "$READ2" > "$SAM"

log_step "Mark duplicates with GATK"
maybe_load_module "${GATK_MODULE:-gatk}"
init_gatk_cmd
run_gatk MarkDuplicatesSpark -I "$SAM" -M "${gene}_dedup_metrics.txt" -O "$DEDUP_BAM"

log_step "Create reference dictionary and FASTA index"
maybe_load_module "${SAMTOOLS_MODULE:-samtools}"
require_cmd "$SAMTOOLS_CMD"
"$SAMTOOLS_CMD" dict "$REF" > "${REF%.fa}.dict"
"$SAMTOOLS_CMD" faidx "$REF"

log_step "Initial HaplotypeCaller & SNP hard-filtering"
run_gatk HaplotypeCaller -R "$REF" -I "$DEDUP_BAM" -O "$RAW_VCF"
run_gatk SelectVariants -R "$REF" -V "$RAW_VCF" -select-type SNP -O "$RAW_SNPS"
run_gatk VariantFiltration \
  -R "$REF" -V "$RAW_SNPS" -O "$FILTERED_SNPS" \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 60.0" \
  -filter-name "MQ_filter" -filter "MQ < 40.0" \
  -filter-name "SOR_filter" -filter "SOR > 4.0" \
  -filter-name "MQRankSum_filter" -filter "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
  -filter-name "ReadPosRankSum_filter" -filter "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0"

log_step "Bootstrapped BQSR from passing SNPs"
run_gatk SelectVariants --exclude-filtered -V "$FILTERED_SNPS" -O "$BQSR_SNPS"
run_gatk BaseRecalibrator -R "$REF" -I "$DEDUP_BAM" --known-sites "$BQSR_SNPS" -O "$RECAL_TABLE"
run_gatk ApplyBQSR -R "$REF" -I "$DEDUP_BAM" -bqsr "$RECAL_TABLE" -O "$RECAL_BAM"

log_step "Recalibrated HaplotypeCaller & final SNP filtering"
run_gatk HaplotypeCaller -R "$REF" -I "$RECAL_BAM" -O "$RAW_RECAL_VCF"
run_gatk SelectVariants -R "$REF" -V "$RAW_RECAL_VCF" -select-type SNP -O "$RAW_RECAL_SNPS"
run_gatk VariantFiltration \
  -R "$REF" -V "$RAW_RECAL_SNPS" -O "$FINAL_SNPS" \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 60.0" \
  -filter-name "MQ_filter" -filter "MQ < 40.0" \
  -filter-name "SOR_filter" -filter "SOR > 4.0" \
  -filter-name "MQRankSum_filter" -filter "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
  -filter-name "ReadPosRankSum_filter" -filter "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0"

log_step "Annotate final SNP VCF with snpEff (${SNPEFF_GENOME})"
if [[ -n "$SNPEFF_CMD" ]]; then
  require_cmd "$SNPEFF_CMD"
  "$SNPEFF_CMD" -Xmx4g "$SNPEFF_GENOME" "$FINAL_SNPS" > "$FINAL_ANN"
elif [[ -n "$SNPEFF_JAR" ]]; then
  maybe_load_module "${JAVA_MODULE:-java}"
  require_cmd "$JAVA_CMD"
  require_file "$SNPEFF_JAR"
  "$JAVA_CMD" -Xmx4g -jar "$SNPEFF_JAR" "$SNPEFF_GENOME" "$FINAL_SNPS" > "$FINAL_ANN"
elif command -v snpEff >/dev/null 2>&1; then
  snpEff -Xmx4g "$SNPEFF_GENOME" "$FINAL_SNPS" > "$FINAL_ANN"
else
  echo "Set SNPEFF_JAR=/path/to/snpEff.jar or ensure 'snpEff' is on PATH." >&2
  exit 127
fi

log_step "Done"
echo "$FINAL_ANN" >&3
