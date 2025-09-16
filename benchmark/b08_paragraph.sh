#!/bin/bash
#SBATCH --job-name=paragraph
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=100G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=3-0:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

SV_VCF=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz
REF=/group/paed_yangwl/t2t-chm13/chm13.v2.fasta

bcftools view -i '(SVTYPE="DEL")' $SV_VCF | bcftools +fill-from-fasta - -- --fasta $REF --column REF > $tmpdir/final_lumpy_gt-smoove.genotyped.fixed.vcf.gz

echo -e "Running Paragraph ... "
multigrmpy.py --input $tmpdir/final_lumpy_gt-smoove.genotyped.fixed.vcf.gz \
              --manifest /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/paragraph/manifest.txt \
              --reference-sequence $REF \
              --output /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/paragraph \
              --max-reads-per-event 600 \
              --threads 20 \
              --read-length 150 \
              --verbose \
              --scratch-dir $tmpdir && echo -e "Paragraph completed successfully! "