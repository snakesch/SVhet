svhet --input test.vcf.gz --bam test.bam \
        --ref /group/paed_yangwl/hs37d5/hs37d5.fa \
        --output test.svhet.vcf.gz \
        --image /lustre1/g/paed_yangwl/snakesch/tools/deepvariant/deepvariant-1.9.0.sif \
        --outdir test_output/ \
        --threads 4 \
        --cipos-tag CIPOS \
        --ciend-tag CIEND
