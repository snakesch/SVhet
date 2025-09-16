#!/bin/bash
#SBATCH --job-name=pangenie-index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=200G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=3-0:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

## This is PanGenie v4.2.1 ran using PanGenie resources released as of v3.1.0
# https://github.com/eblerjana/pangenie/wiki/D:-Running-PanGenie-on-HPRC-data (HPRC-CHM13; 88 haplotypes)

module load singularity

REF="/group/paed_yangwl/t2t-chm13/chm13.v2.fasta"
GRAPH_VCF="/lustre1/g/paed_yangwl/snakesch/tools/pangenie/chm13_cactus_filtered_ids.vcf"
CALLSET_VCF="/lustre1/g/paed_yangwl/snakesch/tools/pangenie/chm13_cactus_filtered_ids_biallelic.vcf.gz"
OUT="/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/pangenie"

## One-off indexing

singularity run \
        -B $(dirname $REF):/mnt/ref \
        -B $(dirname ${GRAPH_VCF}):/mnt/resources \
        -B $OUT:/mnt/out \
        /lustre1/g/paed_yangwl/snakesch/tools/pangenie/pangenie.sif \
        PanGenie-index \
        -r /mnt/ref/$(basename $REF) \
        -v /mnt/resources/$(basename ${GRAPH_VCF}) \
        -t 24 \
        -o /mnt/out/index

