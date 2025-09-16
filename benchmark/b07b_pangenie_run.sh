#!/bin/bash
#SBATCH --job-name=pangenie-run
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=200G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=30
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

## This is PanGenie v4.2.1 ran using PanGenie resources released as of v3.1.0
# https://github.com/eblerjana/pangenie/wiki/D:-Running-PanGenie-on-HPRC-data (HPRC-CHM13; 88 haplotypes)

module load singularity

sample=$1

CALLSET_VCF="/lustre1/g/paed_yangwl/snakesch/tools/pangenie/chm13_cactus_filtered_ids_biallelic.vcf.gz"
RAW_DATA="/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/raw_data"
OUT="/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/pangenie"

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

# Concatenate two FASTQ files
zcat $RAW_DATA/$sample/${sample}_1.fastq.gz $RAW_DATA/$sample/${sample}_2.fastq.gz > $tmpdir/$sample.fastq 

singularity run \
        -B $tmpdir:/tmp \
        -B $OUT:/mnt/out \
        /lustre1/g/paed_yangwl/snakesch/tools/pangenie/pangenie.sif \
        PanGenie \
        -f /mnt/out/index \
        -i /tmp/$sample.fastq \
        -o /mnt/out/$sample-pangenie \
        -j 24 \
        -t 24 

cat $OUT/$sample-pangenie_genotyping.vcf | python /lustre1/g/paed_yangwl/snakesch/tools/pangenie/convert-to-biallelic.py ${CALLSET_VCF} > $OUT/$sample-pangenie-biallelic.vcf