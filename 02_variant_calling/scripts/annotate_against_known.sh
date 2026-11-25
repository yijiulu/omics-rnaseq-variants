#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash run_variant_pipeline.sh <gene> <chromosome> | ./annotate_against_known.sh <dbsnp_vcf.gz> [gene_start] [gene_end]

Examples:
  bash run_variant_pipeline.sh vegfa 6 | ./annotate_against_known.sh /path/to/dbsnp_chr6.vcf.gz
  bash run_variant_pipeline.sh brca1 17 | ./annotate_against_known.sh /path/to/dbsnp_chr17.vcf.gz 43044295 43170245

Standalone usage:
  printf '%s\n' brca1_filtered_snps_final_ann.vcf | ./annotate_against_known.sh /path/to/dbsnp_chr17.vcf.gz

Inputs:
  stdin: path to a snpEff-annotated VCF (e.g. <gene>_filtered_snps_final_ann.vcf)

Outputs:
  <input_without_.vcf>_dbsnp.vcf
  parsed_vcf_file.vcf (written when gene_start and gene_end are supplied)

Environment overrides:
  USE_MODULES=0|1   Load HPC Environment Modules/Lmod if set to 1 (default: 0)
  BGZIP_CMD=bgzip
  TABIX_CMD=tabix
  BCFTOOLS_CMD=bcftools
  PARSED_VCF=parsed_vcf_file.vcf
USAGE
}

if [[ $# -lt 1 || $# -gt 3 ]]; then
  usage >&2
  exit 1
fi

DBSNP_FILE=$1
START=${2:-}
END=${3:-}
PARSED_VCF=${PARSED_VCF:-parsed_vcf_file.vcf}

USE_MODULES=${USE_MODULES:-0}
BGZIP_CMD=${BGZIP_CMD:-bgzip}
TABIX_CMD=${TABIX_CMD:-tabix}
BCFTOOLS_CMD=${BCFTOOLS_CMD:-bcftools}

exec 3>&1
exec 1>&2

log_step() {
  printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"
}

require_file() {
  local file=$1
  if [[ ! -f "$file" ]]; then
    echo "Missing required file: $file" >&2
    exit 1
  fi
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

if [[ -n "$START" && -z "$END" ]] || [[ -z "$START" && -n "$END" ]]; then
  echo "gene_start and gene_end must be provided together." >&2
  usage >&2
  exit 1
fi

if [[ -n "$START" ]]; then
  if [[ ! "$START" =~ ^[0-9]+$ || ! "$END" =~ ^[0-9]+$ ]]; then
    echo "gene_start and gene_end must be integer genomic coordinates." >&2
    exit 1
  fi
  if (( START > END )); then
    echo "gene_start must be <= gene_end." >&2
    exit 1
  fi
fi

if [[ -t 0 ]]; then
  echo "No upstream VCF path received on stdin." >&2
  usage >&2
  exit 1
fi

ANN_VCF=$(awk 'NF { line = $0 } END { print line }')
if [[ -z "$ANN_VCF" ]]; then
  echo "No non-empty upstream VCF path received on stdin." >&2
  exit 1
fi

if [[ "$ANN_VCF" == *.gz ]]; then
  ANN_GZ=$ANN_VCF
  base=${ANN_VCF%.gz}
else
  ANN_GZ="${ANN_VCF}.gz"
  base=$ANN_VCF
fi

if [[ "$base" == *.vcf ]]; then
  DBSNP_ANN="${base%.vcf}_dbsnp.vcf"
else
  DBSNP_ANN="${base}_dbsnp.vcf"
fi

require_file "$ANN_VCF"
require_file "$DBSNP_FILE"

log_step "Compress and index annotated VCF"
maybe_load_module "${HTSLIB_MODULE:-htslib}"
require_cmd "$BGZIP_CMD"
require_cmd "$TABIX_CMD"
if [[ "$ANN_VCF" != *.gz ]]; then
  "$BGZIP_CMD" -c "$ANN_VCF" > "$ANN_GZ"
fi
"$TABIX_CMD" -f -p vcf "$ANN_GZ"

log_step "Add dbSNP rsIDs with bcftools"
maybe_load_module "${BCFTOOLS_MODULE:-bcftools}"
require_cmd "$BCFTOOLS_CMD"
"$BCFTOOLS_CMD" annotate -c ID -a "$DBSNP_FILE" "$ANN_GZ" > "$DBSNP_ANN"

if [[ -n "$START" ]]; then
  log_step "Extract variants in genomic interval ${START}-${END}"
  awk -v start="$START" -v end="$END" '
    /^#/ { print; next }
    {
      pos = $2 + 0
      if (pos >= start && pos <= end) print
    }
  ' "$DBSNP_ANN" > "$PARSED_VCF"
fi

log_step "Done"
echo "$DBSNP_ANN" >&3
if [[ -n "$START" ]]; then
  echo "$PARSED_VCF" >&3
fi
