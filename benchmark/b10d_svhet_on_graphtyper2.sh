#!/bin/bash
#SBATCH --job-name=svhet_on_graphtyper2                 # 1. Job name
#SBATCH --mail-type=FAIL                        # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=snakesch@connect.hku.hk     #    Email address to receive notification
#SBATCH --partition=intel                       # 3. Request a partition
#SBATCH --mem=50G                               #    CPU memory
#SBATCH --qos=normal                            # 4. Request a QoS
#SBATCH --ntasks-per-node=10                     # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                               # 6. Request number of node(s)
#SBATCH --time=0-10:00:00                       # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/log/%x_%j.out                          # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/log/%x_%j.err                           #    Standard error log as $job_name_$job_id.err

## Batch script to run SVhet on a cohort-level VCF

tmpdir=$(mktemp -d); trap 'rm -rf "$tmpdir"' EXIT

## Need to annotate original CIPOS and CIEND
ORIG=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz
FINAL=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/graphtyper2-merged.vcf.gz

bcftools query -f '%CHROM\t%POS\t%INFO/CIPOS\t%INFO/CIEND' $ORIG | bgzip > $tmpdir/cipos-ciend.txt.gz
tabix -s1 -b2 -e2 $tmpdir/cipos-ciend.txt.gz
bcftools view -h $ORIG | egrep -w 'CIPOS|CIEND' > $tmpdir/vcf.header
bcftools annotate --threads 10 \
                  -h $tmpdir/vcf.header \
                  -a $tmpdir/cipos-ciend.txt.gz \
                  -c CHROM,POS,INFO/CIPOS,INFO/CIEND \
                  -Oz -o $tmpdir/graphtyper2-merged-fixed.vcf.gz \
                  $FINAL

time /lustre1/g/paed_yangwl/snakesch/work/SVhet/core/svhet/svhet.sh \
        --ref /group/paed_yangwl/t2t-chm13/chm13.v2.fasta \
        --sv-vcf $tmpdir/graphtyper2-merged-fixed.vcf.gz \
        --outdir /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2+svhet \
        --manifest /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams/manifest.txt \
        --jobs 8 \
        --min-dp 5 \
        --high-het 1