#!/bin/bash
#SBATCH --job-name=svhet_on_duphold                 # 1. Job name
#SBATCH --mail-type=FAIL                        # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=snakesch@connect.hku.hk     #    Email address to receive notification
#SBATCH --partition=intel                       # 3. Request a partition
#SBATCH --mem=50G                               #    CPU memory
#SBATCH --qos=normal                            # 4. Request a QoS
#SBATCH --ntasks-per-node=5                     # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                               # 6. Request number of node(s)
#SBATCH --time=0-10:00:00                       # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/log/%x_%j.out                          # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/log/%x_%j.err                           #    Standard error log as $job_name_$job_id.err

## Batch script to run SVhet on a cohort-level VCF

time /lustre1/g/paed_yangwl/snakesch/work/SVhet/core/svhet/svhet.sh \
        --ref /group/paed_yangwl/t2t-chm13/chm13.v2.fasta \
        --sv-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/duphold/final_lumpy_gt-smoove.genotyped.duphold.vcf.gz \
        --outdir /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/duphold+svhet \
        --manifest /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams/manifest.txt \
        --jobs 8 \
        --min-dp 5 \
        --high-het 1