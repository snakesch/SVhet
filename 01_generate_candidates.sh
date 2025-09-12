#!/bin/bash
# Generate deletion candidates from a cohort VCF, apply additional SV filters,
# then split candidates by |SVLEN|, keeping only final outputs in OUTDIR.

set -euo pipefail

SV_VCF=${1:-}
OUTDIR=${2:-}
LARGE_SVLEN_ABS=${3:-1000000}   # split threshold for |SVLEN| (default 1e6)
BED=${4:-}

if [[ -z "${SV_VCF}" || -z "${OUTDIR}" ]]; then
  echo "Usage: $0 SV_VCF OUTDIR [LARGE_SVLEN_ABS] [BED]" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

# --- helpers ---
bed_chr() {
  local bed="$1"
  if [[ "$bed" =~ .gz$ ]]; then
    zcat "$bed" | cut -f1 | grep -v '^#' | sort -u
  else
    cut -f1 "$bed" | grep -v '^#' | sort -u
  fi
}

vcf_chr() {
  local vcf="$1"
  if [[ "$vcf" =~ .gz$ ]]; then
    zgrep -v '^#' "$vcf" | awk '{print $1}' | sort -u
  else
    awk '!/^#/ {print $1}' "$vcf" | sort -u
  fi
}

check_chr() {
  local bed="$1"
  local vcf="$2"
  local overlap_chr
  overlap_chr=$(comm -12 <(bed_chr "$bed") <(vcf_chr "$vcf") | wc -l)
  if [[ ${overlap_chr} -eq 0 ]]; then
    echo "No candidate found in target regions. Perhaps different CHR?" >&2
    exit 1
  fi
}

# Core expressions (use INFO/ prefix; abs() and 0-based array indexing)
cand_expr='INFO/SVTYPE="DEL" && COUNT(GT="het")>0'  # deletion with at least one het carrier
# Move from candidates -> non-candidates if |SVLEN| < 50, or large CI uncertainty
bcf_expr_move='abs(INFO/SVLEN) < 50 || (abs(INFO/CIPOS[0]) + abs(INFO/CIPOS[1]) > 150) || (abs(INFO/CIEND[0]) + abs(INFO/CIEND[1]) > 150)'
# Split thresholds for final candidates
bcf_expr_large="abs(INFO/SVLEN) >= ${LARGE_SVLEN_ABS}"
bcf_expr_small="abs(INFO/SVLEN) >= 50 && abs(INFO/SVLEN) < ${LARGE_SVLEN_ABS}"

# --- initial candidate / non-candidate (tmpdir only) ---
if [[ -z "${BED:-}" ]]; then
  echo "No target region detected."
  bcftools view -Ov -i "${cand_expr}" "${SV_VCF}" > "${tmpdir}/candidates.vcf"
  bcftools view -Ov -e "${cand_expr}" "${SV_VCF}" > "${tmpdir}/no_candidates.initial.vcf"
else
  echo "Target BED detected - ${BED}"
  check_chr "${BED}" "${SV_VCF}"

  bedtools intersect -a "${SV_VCF}" -b "${BED}" -header -u > "${tmpdir}/in_region.vcf"
  bedtools intersect -a "${SV_VCF}" -b "${BED}" -header -v > "${tmpdir}/out_region.vcf"

  bcftools view -Ov -i "${cand_expr}" "${tmpdir}/in_region.vcf" > "${tmpdir}/candidates.vcf"
  bcftools view -Ov -e "${cand_expr}" "${tmpdir}/in_region.vcf" > "${tmpdir}/no_candidates.in_region.vcf"

  bcftools view -h "${SV_VCF}" > "${tmpdir}/no_candidates.initial.vcf"
  bcftools view -H "${tmpdir}/out_region.vcf" >> "${tmpdir}/no_candidates.initial.vcf"
  bcftools view -H "${tmpdir}/no_candidates.in_region.vcf" >> "${tmpdir}/no_candidates.initial.vcf"
fi

# --- additional filters: move disqualified candidates into non-candidates ---
bcftools view -Ov -i "${bcf_expr_move}" "${tmpdir}/candidates.vcf" > "${tmpdir}/to_move.vcf"
bcftools view -Ov -e "${bcf_expr_move}" "${tmpdir}/candidates.vcf" > "${tmpdir}/candidates.filtered.vcf"

# Compose final non-candidates (single header)
bcftools view -h "${SV_VCF}" > "${tmpdir}/no_candidates.vcf"
bcftools view -H "${tmpdir}/no_candidates.initial.vcf" >> "${tmpdir}/no_candidates.vcf"
bcftools view -H "${tmpdir}/to_move.vcf" >> "${tmpdir}/no_candidates.vcf"
bcftools sort -Ob -o "${OUTDIR}/no_candidates.bcf" "${tmpdir}/no_candidates.vcf" 2/dev/null && bcftools index "${OUTDIR}/no_candidates.bcf"

# --- final split of candidates by |SVLEN| ---
bcftools view -Ob -i "${bcf_expr_large}" "${tmpdir}/candidates.filtered.vcf" > "${OUTDIR}/large_candidates.bcf"
bcftools view -Ob -i "${bcf_expr_small}" "${tmpdir}/candidates.filtered.vcf" > "${OUTDIR}/small_candidates.bcf"
bcftools index "${OUTDIR}/large_candidates.bcf"
bcftools index "${OUTDIR}/small_candidates.bcf"

# --- reporting ---
echo "Small candidates: $(bcftools view -H "${OUTDIR}/small_candidates.bcf" | wc -l)"
echo "Large candidates: $(bcftools view -H "${OUTDIR}/large_candidates.bcf" | wc -l)"
echo "Non-candidates:   $(bcftools view -H "${OUTDIR}/no_candidates.bcf" | wc -l)"