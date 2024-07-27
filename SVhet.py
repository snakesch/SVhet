#!/usr/bin/env python3

## Driver code

import argparse
import os
import logging
import subprocess
import pandas as pd
import multiprocessing as mp

import pysam

from src.utils import loadVCF, create_variants, write_to_bed
from src.heterozygosity import filter_by_presence_of_het
from src.homozygosity import get_heterozygous_sites_from_bam

from src.constants import VALID_CONTIGS

# for benchmarking
# from perf import write_candidates_to_dataframe

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Evaluate structural and short variant calls')

    # Add the arguments
    parser.add_argument('--svcall', required=True, type=str, help='Candidate structural variant VCF')
    parser.add_argument('--short-var', required=True, type=str, help='Short variant calls')
    parser.add_argument('--bam', required=True, type=str, help='Sequence alignments')
    parser.add_argument('--output', required=True, type=str, help='Output directory')
    parser.add_argument('--problematic-regions', required=True, type=str, help="Problematic regions (BED)")
    parser.add_argument('--repeat-regions', type=str, required=True, help="Regions of known repeat elements >100bp (BED)")
    parser.add_argument("--verbose", default = "INFO", help = "level of verbosity (default: INFO)")
    parser.add_argument("--threads", type = int, default = 4, help = "Number of CPU threads (default: 4)")

    # Parse the arguments
    args = parser.parse_args()
    
    # Initialize a logger object
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = args.verbose.upper())
    logger = logging.getLogger("root")

    # Validate the input files
    if not os.path.isfile(args.svcall):
        raise FileNotFoundError(f"The file '{args.svcall}' does not exist.")
    if not os.path.isfile(args.short_var):
        raise FileNotFoundError(f"The file '{args.short_var}' does not exist.")
    if not os.path.isfile(args.bam):
        raise FileNotFoundError(f"The file '{args.bam}' does not exist.")
    if not os.path.isfile(args.problematic_regions):
        raise FileNotFoundError(f"The file '{args.problematic_regions}' does not exist.")
    if not os.path.isfile(args.repeat_regions):
        raise FileNotFoundError(f"The file '{args.repeat_regions}' does not exist.")

    # Check if the output directory exists, and create it if it doesn't
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Now you can use the parsed arguments in your code
    logger.info(f"Structural variant calls: {args.svcall}")
    logger.info(f"Short variant calls: {args.short_var}")
    logger.info(f"BAM file: {args.bam}")
    logger.info(f"Output directory: {args.output}")
    
    basename = os.path.basename(args.svcall).replace(".gz", "").replace(".vcf", "").replace(".bcf", "")
    
    # Filters out non-HDELs before loading
    proc = subprocess.run(f"bcftools filter -i 'SVTYPE~\"DEL\"' {args.svcall} | bcftools filter -i 'GT=\"het\"' -Oz -o {os.path.join(args.output, basename + '.hdel.vcf.gz')} && tabix -f -p vcf {os.path.join(args.output, basename + '.hdel.vcf.gz')} ", shell=True)

    if proc.returncode != 0:
        raise RuntimeError("Unable to filter out HDELs! ")
        
    sv_calls = create_variants(os.path.join(args.output, basename + '.hdel.vcf.gz'))
    sv_calls = [ call for call in sv_calls if call.chrom in VALID_CONTIGS]
    
    logger.info(f"Loaded {len(sv_calls)} HDEL variants. ")
    
    problematic_regions = args.problematic_regions
    repeat_regions = args.repeat_regions

    # Filter out short variants in problematic regions
    proc = subprocess.run(f"bcftools filter --threads {args.threads} --mask-file {problematic_regions} --soft-filter DIFFICULT {args.short_var} | bcftools filter --threads {args.threads} -i 'GT=\"het\"' | bcftools filter -Wtbi --threads {args.threads} -e 'FILTER~\"DIFFICULT\"' -Oz -o {os.path.join(args.output, basename + '.problematic_regions.high_confidence_hets.vcf.gz')} ", shell=True)
    
    if proc.returncode != 0:
        raise RuntimeError("Unable to filter out hets in problematic regions! ")
        
    proc = subprocess.run(f"bcftools filter --threads {args.threads} --mask-file {repeat_regions} --soft-filter DIFFICULT {os.path.join(args.output, basename + '.problematic_regions.high_confidence_hets.vcf.gz')} | bcftools filter --threads {args.threads} -i 'GT=\"het\"' | bcftools filter -Wtbi --threads {args.threads} -e 'FILTER~\"DIFFICULT\"' -Oz -o {os.path.join(args.output, basename + '.high_confidence_hets.vcf.gz')} ", shell=True)
    
    if proc.returncode != 0:
        raise RuntimeError("Unable to filter out hets in repeat regions! ")
    
    tol_variants_n = len(sv_calls)
    
    # Filters out true het events
    with pysam.VariantFile(os.path.join(args.output, basename + '.high_confidence_hets.vcf.gz')) as short_vars:
        sv_filtered_true_het = filter_by_presence_of_het(sv_calls, short_vars, problematic_regions)
    logger.info(f"Step 1 (Filter by heterozygosity): Done (Filtered {tol_variants_n - len(sv_filtered_true_het)} variants)")
    
    sv_filtered_3prime, sv_filtered_3prime_5prime = [], []
    
    for region in sv_filtered_true_het:
        if region.het_cnt == 0:
            # 3prime breakpoint
            het = get_heterozygous_sites_from_bam(args.bam, region, breakpoint = "3prime", repeat_regions = args.repeat_regions)
            if len(het) <= 1 or len(het) > 5:
                sv_filtered_3prime.append(region)
        else:
            sv_filtered_3prime.append(region)
    logger.info(f"Step 2 (Filter by 3' homozygosity): Done (Filtered {len(sv_filtered_true_het) - len(sv_filtered_3prime)} variants)")
    
    for region in sv_filtered_3prime:
        if region.het_cnt == 0:
            # 5prime breakpoint
            het = get_heterozygous_sites_from_bam(args.bam, region, breakpoint = "5prime", repeat_regions = args.repeat_regions)
            if len(het) <= 2:
                sv_filtered_3prime_5prime.append(region)
        else:
            sv_filtered_3prime_5prime.append(region)
    logger.info(f"Step 3 (Filter by 5' homozygosity): Done (Filtered {len(sv_filtered_3prime) - len(sv_filtered_3prime_5prime)} variants)")

    write_to_bed(sv_filtered_3prime_5prime, args.output, basename + ".sv_filtered_3prime_5prime.bed", sep="\t", index=False, header=False)
    # setting f=0.95 is necessary here to include HDELs spanning (start, end) only
    proc = subprocess.run(f"bedtools intersect -a {os.path.join(args.output, basename + '.hdel.vcf.gz')} -b {os.path.join(args.output, basename + '.sv_filtered_3prime_5prime.bed')} -f 0.95 -r -header | bgzip -c > {os.path.join(args.output, basename + '.filtered.hdel.vcf.gz')}", shell=True)
    
    proc = subprocess.run(f"bcftools index --tbi {os.path.join(args.output, basename + '.filtered.hdel.vcf.gz')}", shell=True)
    
    proc = subprocess.run(f"bcftools filter --threads {args.threads} -Wcsi -e '(GT=\"het\") & (SVTYPE~\"DEL\")' -Ob -o {os.path.join(args.output, basename + '.non_hdels.bcf')} {args.svcall} && bcftools concat -Wtbi -a -Oz -o {os.path.join(args.output, basename + '.filtered.vcf.gz')} {os.path.join(args.output, basename + '.filtered.hdel.vcf.gz')} {os.path.join(args.output, basename + '.non_hdels.bcf')} 2>/dev/null ", shell=True)
    logger.info("Step 4 (Consolidating candidate SV calls): Done ")
    
    os.remove(os.path.join(args.output, basename + '.filtered.hdel.vcf.gz'))
    os.remove(os.path.join(args.output, basename + '.filtered.hdel.vcf.gz.tbi'))
    os.remove(os.path.join(args.output, basename + '.non_hdels.bcf'))
    os.remove(os.path.join(args.output, basename + '.non_hdels.bcf.csi'))
    os.remove(os.path.join(args.output, basename + '.hdel.vcf.gz'))
    os.remove(os.path.join(args.output, basename + '.hdel.vcf.gz.tbi'))
    os.remove(os.path.join(args.output, basename + '.sv_filtered_3prime_5prime.bed'))
    os.remove(os.path.join(args.output, basename + '.problematic_regions.high_confidence_hets.vcf.gz'))
    os.remove(os.path.join(args.output, basename + '.problematic_regions.high_confidence_hets.vcf.gz.tbi'))
    os.remove(os.path.join(args.output, basename + '.high_confidence_hets.vcf.gz'))
    os.remove(os.path.join(args.output, basename + '.high_confidence_hets.vcf.gz.tbi'))
    
    return sv_calls, sv_filtered_true_het, sv_filtered_3prime, sv_filtered_3prime_5prime

if __name__ == "__main__":
    main()