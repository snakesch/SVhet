#!/bin/bash

# Input VCF (bgzipped and indexed)
VCF=$1

OUT_TYPE_FN=${2:-"sv-summary.txt"}
OUTD="/lustre1/g/paed_yangwl/snakesch/work/manuscript_figure_table/svhet/tables"

# Extract sample IDs
samples=$(bcftools query -l /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz)

echo -e "Sample\tHET_DEL\tHOM_DEL\tHET_INS\tHOM_INS\tHET_DUP\tHOM_DUP" > "${OUTD}/${OUT_TYPE_FN}"
for sample in $samples; do
    echo -e "Summarizing SVTYPE for sample $sample " 
    het_del=$(bcftools query -s $sample -i 'GT="het" && INFO/SVTYPE=="DEL"' -f '.' $VCF | wc -l)
    hom_del=$(bcftools query -s $sample -i 'GT="hom" && INFO/SVTYPE=="DEL"' -f '.' $VCF | wc -l)
    het_ins=$(bcftools query -s $sample -i 'GT="het" && INFO/SVTYPE=="INS"' -f '.' $VCF | wc -l)
    hom_ins=$(bcftools query -s $sample -i 'GT="hom" && INFO/SVTYPE=="INS"' -f '.' $VCF | wc -l)
    het_dup=$(bcftools query -s $sample -i 'GT="het" && INFO/SVTYPE=="DUP"' -f '.' $VCF | wc -l)
    hom_dup=$(bcftools query -s $sample -i 'GT="hom" && INFO/SVTYPE=="DUP"' -f '.' $VCF | wc -l)
    echo -e "$sample\t$het_del\t$hom_del\t$het_ins\t$hom_ins\t$het_dup\t$hom_dup" >> "${OUTD}/${OUT_TYPE_FN}"
done


