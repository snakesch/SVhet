#!/bin/bash
#SBATCH --job-name=graphtyper2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=100G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=2-0:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

# GraphTyper v2.7.7

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

chr="$1" ## chr1 .. chrY

## Change --avg_cov_by_readlen to improve speed
for i in {1..31}; do echo "0.2" >> $tmpdir/avg-cov-by-readlen.txt; done

[[ -d /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/$chr/ ]] || mkdir -p /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/$chr/; 

graphtyper genotype_sv /group/paed_yangwl/t2t-chm13/chm13.v2.fasta \
                       /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz \
                       --output /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/$chr/ \
                       --sams /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/graphtyper2/bam_list.txt \
                       --threads 20 \
                       --region $chr \
                       --avg_cov_by_readlen $tmpdir/avg-cov-by-readlen.txt \
                       --verbose



