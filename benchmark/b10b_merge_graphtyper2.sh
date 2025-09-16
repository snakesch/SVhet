#!/bin/bash
#SBATCH --job-name=merge-graphtyper2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=100G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=0-1:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

cd /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2

for i in {1..22} X Y
do
    find $PWD/chr$i/chr$i/ -type f -name "*.vcf.gz" | sort > $tmpdir/chr$i-to-merge.txt
    bcftools concat -Wcsi --file-list $tmpdir/chr$i-to-merge.txt -Ob -o $tmpdir/graphtyper-chr$i.bcf --threads 10
done

find $tmpdir -type f -name "graphtyper-chr*.bcf" | sort > $tmpdir/all-to-merge.txt
bcftools concat -Wtbi --file-list $tmpdir/all-to-merge.txt -Oz -o graphtyper2-merged.vcf.gz --threads 10