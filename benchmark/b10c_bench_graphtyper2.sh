#!/bin/bash
#SBATCH --job-name=single_truvari                # 1. Job name
#SBATCH --mail-type=FAIL                        # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=snakesch@connect.hku.hk     #    Email address to receive notification
#SBATCH --partition=intel                       # 3. Request a partition
#SBATCH --mem=50G                               #    CPU memory
#SBATCH --qos=normal                            # 4. Request a QoS
#SBATCH --ntasks-per-node=5                     # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                               # 6. Request number of node(s)
#SBATCH --time=0-00:10:00                       # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=/dev/null                      # We do not need stdout/stderr from this job
#SBATCH --error=/dev/null

sample=$1

TRUTH=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/benchmark/variants_T2T-CHM13_sv_insdel_sym_HGSVC2024v1.0.annot.vcf.gz
# ORIG=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz
FINAL=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/graphtyper2-merged.vcf.gz

OUT="/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/truvari"

[[ -d $OUT ]] || mkdir -p $OUT;
[[ -f ${FINAL}.tbi ]] || tabix -p vcf $FINAL;

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

## Focus on heterozygous deletions
bcftools view -s $sample $FINAL | bcftools filter -i '(SVTYPE="DEL" & FMT/GT = "het" & FMT/FT = "PASS")' -Wtbi -Oz -o $tmpdir/$sample-final.vcf.gz
bcftools view -s $sample $TRUTH | bcftools filter -i '(SVTYPE="DEL" & FMT/GT = "het")' -Wtbi -Oz -o $tmpdir/$sample-truth.vcf.gz

echo -e "Running truvari on sample $sample" 

## truvari v5.2 
[[ -d $OUT/$sample ]] || mkdir -p $OUT/$sample;
[[ -d $OUT/$sample/final ]] && rm -rf $OUT/$sample/final ;

mamba run -n svhet_benchmark truvari bench \
            --base $tmpdir/$sample-truth.vcf.gz \
            --comp $tmpdir/$sample-final.vcf.gz \
            --pctseq 0.0 \
            --pctsize 0.7 \
            --pctovl 0.7 \
            --refdist 500 \
            --passonly \
            --sizemin 50 \
            --sizemax 2000000 \
            --no-ref a \
            --output $OUT/$sample/final/

## Builds joblibs          
mamba run -n svhet_benchmark truvari vcf2df \
          --bench-dir $OUT/$sample/final/ \
          $OUT/$sample/final/data.jl 
