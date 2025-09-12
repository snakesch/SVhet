#!/bin/bash

## Given read evidences generated from step 2, perform variant calling using bcftools.

REF=$1
BAM=$2
JOBS=$3

## Variant calling using bcftools
bcftools mpileup --threads $JOBS -f $REF $BAM 2>/dev/null | bcftools call --threads $JOBS -mv --ploidy 2 | bcftools filter -i 'GT="het"' -Ob -o "${BAM/.bam/.bcf}"

bcftools index "${BAM/.bam/.bcf}"