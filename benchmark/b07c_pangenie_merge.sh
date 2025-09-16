#!/bin/bash
#SBATCH --job-name=pangenie-merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=snakesch@connect.hku.hk
#SBATCH --partition=amd
#SBATCH --mem=100G
#SBATCH --qos=normal
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=0-01:00:00
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.out ## CHANGE
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/log/%x_%j.err ## CHANGE

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

cd /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/pangenie

samples=( $(find . -name "*-biallelic.vcf" -type f -exec basename {} -pangenie-biallelic.vcf \; | sort) )

## Extra headers
HEADER=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/pangenie/vcf.header
echo -e "Outputting to $tmpdir" 

for sample in ${samples[@]};
do  
    if [[ -f $sample-sv-annot.bcf ]]
    then
        continue
    fi

    echo -e "Processing sample $sample "

    bcftools reheader --header $HEADER --samples-list $sample --threads 10 $sample-pangenie-biallelic.vcf | \
    bcftools norm -m - | \
    bcftools view --threads 8 -i '(abs(STRLEN(REF) - STRLEN(ALT[0])) > 50 && AC > 0)' -Wcsi -Ob -o $tmpdir/$sample-sv.bcf 

    bcftools query -f '%ID' $tmpdir/$sample-sv.bcf | \
        awk -v FS="-" -v OFS="\t" '($3 == "DEL") { print $1, $2, $3, -$5; } ($3 != "INS") { print $1, $2, $3, $5;} ' \
        | bgzip > $tmpdir/$sample-annotations.txt.gz && tabix -s1 -b2 -e2 $tmpdir/$sample-annotations.txt.gz

    bcftools annotate \
                      -a $tmpdir/$sample-annotations.txt.gz \
                      -c CHROM,POS,INFO/SVTYPE,INFO/SVLEN \
                      -Wcsi -Ob -o $sample-sv-annot.bcf \
                      --threads 10 \
                      $tmpdir/$sample-sv.bcf

done

## Merging will lead to multiallelic records
# find $PWD -name "*-sv-annot.bcf" | sort > $tmpdir/to_merge.txt
# bcftools merge -m all --file-list $tmpdir/to_merge.txt -Wtbi -Oz -o $tmpdir/pangenie-biallelic-merged.vcf.gz --threads 10

# bcftools +fill-tags -Oz -o pangenie-biallelic-merged.vcf.gz $tmpdir/pangenie-biallelic-merged.vcf.gz -- -t INFO/END
