#!/bin/bash

## Entry point for SVhet

set -euo pipefail
tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT
 
REF=""
SV_VCF=""
BED=""
OUTDIR=""
MANIFEST=""
JOBS=1
KEEP_INTERMEDIATE=false
MIN_DP=5
HIGH_HET=1

usage() {
    cat <<EOF
Usage: $0 --ref <reference fasta> --sv-vcf <cohort structural variant VCF> --outdir <output directory> --manifest <manifest file> [OPTIONS]

This script processes cohort-level structural variant callsets using heterozygous sites.

Required arguments:
  --ref                Path to the reference FASTA file (e.g., /path/to/chm13.v2.fasta)
  --sv-vcf             Path to the cohort structural variant VCF file (e.g., /path/to/lumpy.vcf.gz)
  --outdir             Output directory for results (e.g., /path/to/output)
  --manifest           Path to the manifest file listing BAM paths (e.g., /path/to/manifest.txt)

Optional arguments:
  --bed                Path to the target region BED file (default: none)
  --jobs               Number of parallel jobs to run (default: 1)
  --keep-intermediate  Keep intermediate files
  --min-dp             Minimum alternate depth threshold (default: 5)
  --high-het           Minimum HET count to reject a DEL (default: 1)

  --help, -h           Print this help message

EOF
    exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
        case "$1" in
                --ref)
                        REF="$2"
                        shift 2
                        ;;
                --sv-vcf)
                        SV_VCF="$2"
                        shift 2
                        ;;
                --bed)
                        BED="$2"
                        shift 2
                        ;;
                --outdir)
                        OUTDIR="$2"
                        shift 2
                        ;;
                --manifest)
                        MANIFEST="$2"
                        shift 2
                        ;;
                --jobs)
                        JOBS="$2"
                        shift 2
                        ;;
                --keep-intermediate)
                        KEEP_INTERMEDIATE=true
                        shift 1
                        ;;
                --min-dp)
                        MIN_DP="$2"
                        shift 2
                        ;;
                --high-het)
                        HIGH_HET="$2"
                        shift 2
                        ;;
                -h|--help)
                        usage
                        ;;
                *)
                        echo "Unknown option: $1"
                        # usage
                        ;;
        esac
done

# Check required arguments
if [[ -z "$REF" || -z "$SV_VCF" || -z "$OUTDIR" || -z "$MANIFEST" ]]; then
        echo -e "Error: Missing required arguments.\n"
        usage
fi

cat << 'EOF'

███████╗██╗   ██╗██╗  ██╗███████╗████████╗
██╔════╝██║   ██║██║  ██║██╔════╝╚══██╔══╝
███████╗██║   ██║███████║█████╗     ██║   
╚════██║╚██╗ ██╔╝██╔══██║██╔══╝     ██║   
███████║ ╚████╔╝ ██║  ██║███████╗   ██║   
╚══════╝  ╚═══╝  ╚═╝  ╚═╝╚══════╝   ╚═╝   
                                          
EOF

echo -e "Version v0.2.0 - Sept 6, 2025\n"

echo -e "Running SVhet with arguments:                                  "
echo -e "       REF                     : $REF                          "
echo -e "       Cohort VCF              : $SV_VCF                       "
echo -e "       Cohort manifest         : $MANIFEST                     "
echo -e "       Output directory        : $OUTDIR                       "
echo -e "       Target region           : $BED                          "
echo -e "       Minimum DP              : $MIN_DP                       "
echo -e "       HET count allowed       : $HIGH_HET                     "
echo -e "       Jobs                    : $JOBS                         "
echo -e "       Keep intermediates      : $KEEP_INTERMEDIATE          \n"

## Switch to script directory
OLDPWD=$PWD
cd $(dirname $0)

echo -e "## - Cohort level statistics - ##\n"
bash 01_generate_candidates.sh $SV_VCF $OUTDIR 1000000 $BED

echo -e "\n## - Sample level statistics - ##"
while read -r sample bam bai
do
    echo -e "Processing $sample"
    python 02_filter_by_samples.py \
            --sample $sample \
            --bam $bam \
            --cohort-vcf $OUTDIR/small_candidates.bcf \
            --out $OUTDIR \
            --jobs $JOBS \
            --mode small # 2>/dev/null

    python 02_filter_by_samples.py \
            --sample $sample \
            --bam $bam \
            --cohort-vcf $OUTDIR/large_candidates.bcf \
            --out $OUTDIR \
            --jobs $JOBS \
            --mode large # 2>/dev/null

    ## Generate 4 short variant callsets
    echo -e "Calling variants ... "

    printf '%s\n' \
    "bash 03_call_variants.sh \"$REF\" \"${OUTDIR}/$sample-wt-reads-small.bam\" $JOBS" \
    "bash 03_call_variants.sh \"$REF\" \"${OUTDIR}/$sample-wt-reads-large.bam\" $JOBS" \
    "bash 03_call_variants.sh \"$REF\" \"${OUTDIR}/$sample-mut-reads-small.bam\" $JOBS" \
    "bash 03_call_variants.sh \"$REF\" \"${OUTDIR}/$sample-mut-reads-large.bam\" $JOBS" \
    | xargs -I {} -P $JOBS bash -c '{}'

    ## Evaluate heterozygosity
    echo -e "Evaluating heterozygous sites ... "
    python 04_het_evaluator.py \
                --candidate-vcf $OUTDIR/small_candidates.bcf \
                --wt-vcf "${OUTDIR}/$sample"-wt-reads-small.bcf \
                --mut-vcf "${OUTDIR}/$sample"-mut-reads-small.bcf \
                --sample $sample \
                --out "${OUTDIR}/$sample"-annotated-small.bcf \
                --min-dp $MIN_DP \
                --high-het $HIGH_HET

    python 04_het_evaluator.py \
                --candidate-vcf $OUTDIR/large_candidates.bcf \
                --wt-vcf "${OUTDIR}/$sample"-wt-reads-large.bcf \
                --mut-vcf "${OUTDIR}/$sample"-mut-reads-large.bcf \
                --sample $sample \
                --out "${OUTDIR}/$sample"-annotated-large.bcf \
                --min-dp $MIN_DP \
                --high-het $HIGH_HET
    
    ## Index 
    bcftools index "${OUTDIR}/$sample"-annotated-small.bcf
    bcftools index "${OUTDIR}/$sample"-annotated-large.bcf

    ## Merge annotated VCFs for current sample
    bcftools view -s $sample -Ob -o ${tmpdir}/${sample}-no_candidates.bcf "${OUTDIR}/no_candidates.bcf" && bcftools index ${tmpdir}/${sample}-no_candidates.bcf
    bcftools concat -a -Ob -o "${OUTDIR}/$sample"-annotated.bcf "${OUTDIR}/$sample"-annotated-small.bcf "${OUTDIR}/$sample"-annotated-large.bcf ${tmpdir}/${sample}-no_candidates.bcf 2>/dev/null
    bcftools index "${OUTDIR}/$sample"-annotated.bcf 
    rm -f ${tmpdir}/${sample}-no_candidates.bcf

    ## Cleanup 
    if [[ $KEEP_INTERMEDIATE != "true" ]]
    then
        rm -f "${OUTDIR}/$sample"-*.bam* "${OUTDIR}/$sample"-*reads-*.bcf* 
        rm -f "${OUTDIR}/$sample"-annotated-small.bcf* "${OUTDIR}/$sample"-annotated-large.bcf*
    fi

    echo -e "Done for sample $sample \n"

done < $MANIFEST

## Merge across samples
bcftools merge -Oz -o "${OUTDIR}"/final-annotated.vcf.gz --force-single "${OUTDIR}"/*-annotated.bcf 2>/dev/null && bcftools index --tbi "${OUTDIR}"/final-annotated.vcf.gz
echo -e "SVhet done for all samples! "

## Final cleanup step
[[ $KEEP_INTERMEDIATE != "true" ]] && rm -f "${OUTDIR}"/no_candidates.bcf*;
[[ $KEEP_INTERMEDIATE != "true" ]] && rm -f "${OUTDIR}"/small_candidates.bcf*;
[[ $KEEP_INTERMEDIATE != "true" ]] && rm -f "${OUTDIR}"/*-annotated.bcf*;

cd $OLDPWD ## Return to original directory
