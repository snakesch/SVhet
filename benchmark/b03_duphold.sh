#!/bin/bash
#SBATCH --job-name=duphold
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=200G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=25
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=3-0:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

## This script runs other SV filtering tools for comparison

module load singularity/3.8.0

SV_VCF="/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz"

BAM_DIR="/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams"
OUT_DIR="/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/duphold"
REF="/group/paed_yangwl/t2t-chm13/chm13.v2.fasta"
SAMPLES=( $(find ${BAM_DIR} -name "*.bam" -exec basename {} \; | sort) )

singularity run \
        -B ${BAM_DIR}:/mnt/bam \
        -B $(dirname $SV_VCF):/mnt/data \
        -B $(dirname ${REF}):/mnt/ref \
        -B /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/vcfs:/mnt/snps \
        -B ${OUT_DIR}:/mnt/duphold \
        --cleanenv \
        /lustre1/g/paed_yangwl/snakesch/tools/smoove_latest.sif \
        smoove duphold \
        --vcf /mnt/data/$(basename $SV_VCF) \
        --fasta /mnt/ref/$(basename $REF) \
        -p 30 \
        --snps /mnt/snps/1kg_joint_called.vcf.gz \
        --outvcf /mnt/duphold/final_lumpy_gt-smoove.genotyped.duphold.vcf.gz \
        /mnt/bam/HG00096.sorted.markdup.bam \
        /mnt/bam/HG00171.sorted.markdup.bam \
        /mnt/bam/HG00268.sorted.markdup.bam \
        /mnt/bam/HG00358.sorted.markdup.bam \
        /mnt/bam/HG00513.sorted.markdup.bam \
        /mnt/bam/HG00731.sorted.markdup.bam \
        /mnt/bam/HG00732.sorted.markdup.bam \
        /mnt/bam/HG01596.sorted.markdup.bam \
        /mnt/bam/HG01890.sorted.markdup.bam \
        /mnt/bam/HG02282.sorted.markdup.bam \
        /mnt/bam/HG02666.sorted.markdup.bam \
        /mnt/bam/HG02769.sorted.markdup.bam \
        /mnt/bam/HG02953.sorted.markdup.bam \
        /mnt/bam/HG03009.sorted.markdup.bam \
        /mnt/bam/HG03452.sorted.markdup.bam \
        /mnt/bam/HG03520.sorted.markdup.bam \
        /mnt/bam/NA18534.sorted.markdup.bam \
        /mnt/bam/NA18939.sorted.markdup.bam \
        /mnt/bam/NA18989.sorted.markdup.bam \
        /mnt/bam/NA19036.sorted.markdup.bam \
        /mnt/bam/NA19129.sorted.markdup.bam \
        /mnt/bam/NA19238.sorted.markdup.bam \
        /mnt/bam/NA19239.sorted.markdup.bam \
        /mnt/bam/NA19317.sorted.markdup.bam \
        /mnt/bam/NA19331.sorted.markdup.bam \
        /mnt/bam/NA19347.sorted.markdup.bam \
        /mnt/bam/NA19384.sorted.markdup.bam \
        /mnt/bam/NA19434.sorted.markdup.bam \
        /mnt/bam/NA20355.sorted.markdup.bam \
        /mnt/bam/NA20509.sorted.markdup.bam \
        /mnt/bam/NA20847.sorted.markdup.bam

