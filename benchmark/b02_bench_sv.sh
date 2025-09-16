#!/bin/bash
#SBATCH --job-name=single_truvari                # 1. Job name
#SBATCH --mail-type=FAIL                        # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=snakesch@connect.hku.hk     #    Email address to receive notification
#SBATCH --partition=intel                       # 3. Request a partition
#SBATCH --mem=50G                               #    CPU memory
#SBATCH --qos=normal                            # 4. Request a QoS
#SBATCH --ntasks-per-node=5                     # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                               # 6. Request number of node(s)
#SBATCH --time=0-00:15:00                       # 7. Job execution duration limit day-hour:min:sec

sample=$1

TRUTH=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/benchmark/variants_T2T-CHM13_sv_insdel_sym_HGSVC2024v1.0.annot.vcf.gz
ORIG=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz
FINAL=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/final-annotated.vcf.gz

[[ -f ${ORIG}.tbi ]] || tabix -p vcf $ORIG;
[[ -f ${FINAL}.tbi ]] || tabix -p vcf $FINAL;

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

## TODO: Focus on deletions ?
bcftools view -s $sample $ORIG | bcftools filter -i '(SVTYPE="DEL" & FMT/GT = "het")' -Wtbi -Oz -o $tmpdir/$sample-orig.vcf.gz
bcftools view -s $sample $FINAL | bcftools filter -i '(SVTYPE="DEL" & FMT/GT = "het")' | bcftools filter -e '(SVHET="HIGH_HET")' -Wtbi -Oz -o $tmpdir/$sample-final.vcf.gz
bcftools view -s $sample $TRUTH | bcftools filter -i '(SVTYPE="DEL" & FMT/GT = "het")' -Wtbi -Oz -o $tmpdir/$sample-truth.vcf.gz

echo -e "Running truvari on sample $sample" 

## truvari v5.2 
[[ -d /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/ ]] || mkdir -p /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/ ;
[[ -d /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig ]] && rm -rf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig ;
[[ -d /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final ]] && rm -rf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final ;

mamba run -n svhet_benchmark truvari bench \
            --base $tmpdir/$sample-truth.vcf.gz \
            --comp $tmpdir/$sample-orig.vcf.gz \
            --pctseq 0.0 \
            --pctsize 0.7 \
            --pctovl 0.7 \
            --refdist 500 \
            --passonly \
            --sizemin 50 \
            --sizemax 2000000 \
            --no-ref a \
            --output /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig/

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
            --output /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final/

## Builds joblibs 
mamba run -n svhet_benchmark truvari vcf2df \
          --bench-dir /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig/ \
          /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig/data.jl 
          
mamba run -n svhet_benchmark truvari vcf2df \
          --bench-dir /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final/ \
          /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final/data.jl 

## Print summary statistics for comparison
# ORIG_FP=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig/fp.vcf.gz
# ORIG_FN=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/orig/fn.vcf.gz
# FINAL_FP=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final/fp.vcf.gz
# FINAL_FN=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/truvari/$sample/final/fn.vcf.gz

# echo -e "Correctly filtered $(bcftools isec -n-1 $ORIG_FP $FINAL_FP 2>/dev/null | grep -cv ^#) candidates. "
# echo -e "Incorrectly filtered $(bcftools isec -n-1 $ORIG_FN $FINAL_FN 2>/dev/null | grep -cv ^#) candidates. "
