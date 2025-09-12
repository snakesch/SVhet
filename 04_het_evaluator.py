#!/usr/bin/env python3

## This script takes a candidate VCF and its corresponding short variant callset from step 3.
## Returns an annotated candidate VCF for merging.

import os
import argparse

import pysam

def get_reliable_hets(vcf, chrom, start, stop, min_dp):
    '''
    Arguments:
        vcf: sample-level short variant VCF (parsed using pysam)
        chrom: candidate SV chrom
        start: candidate SV start
        stop: candidate SV stop

    Returns:
        Coordinates of reliable HETs: (chrom, pos, dp)

    '''

    ## Placeholder for results
    reliable_hets = []
    
    ## We consider +/- 100bp from predicted start and stop
    widened_start, widened_stop = start - 100, stop + 100
    
    last_seen = -500 ## dummy
    for variant in vcf.fetch(chrom, widened_start, widened_stop):
        sample_id = variant.samples.keys()[0]
        if min(variant.samples[sample_id]["AD"]) >= min_dp: #  and variant.pos - last_seen > 100:
            if abs(variant.pos - last_seen) > 100: # Pop densely distributed HETs as they may have originated from repeat regions
                reliable_hets.append((chrom, variant.pos))
            else:
                if (chrom, last_seen) in reliable_hets:
                    reliable_hets.remove((chrom, last_seen))
            last_seen = variant.pos
    
    return reliable_hets

def filter_main(candidate_vcf, wt_vcf, mut_vcf, min_dp, high_het, sample_id, outpath):

    wt, mut = pysam.VariantFile(wt_vcf), pysam.VariantFile(mut_vcf)
    
    candidates = pysam.VariantFile(candidate_vcf)
    candidates.subset_samples([sample_id])
    
    ## Prepare output VCF header (to be written to the same directory as wt_vcf and mut_vcf)
    out_vcf_path = outpath # os.path.join( os.path.dirname(wt_vcf), sample_id + ".annotated.small.bcf" )
    header = candidates.header
    header.formats.add("WT_HETS", 1, "Integer", "Number of observed WT heterozygous sites")
    header.formats.add("MUT_HETS", 1, "Integer", "Number of observed MUT heterozygous sites")
    header.formats.add("SVHET", 1, "String", "SVhet filter")
    
    out_vcf = pysam.VariantFile(out_vcf_path, mode="w", header=header)
    
    for candidate in candidates.fetch():
        chrom, start, stop = candidate.chrom, candidate.pos, candidate.stop
        genotype = candidate.samples[sample_id]["GT"]
        
        wt_hets, mut_hets = [], []
        if genotype == (0, 1):
            wt_hets = get_reliable_hets(wt, chrom, start, stop, min_dp=min_dp)
            mut_hets = get_reliable_hets(mut, chrom, start, stop, min_dp=min_dp)
            # print(f"Identified HET deletion : {chrom}:{start}-{stop} ")
            # print(f"wt_hets : {len(wt_hets)}; mut_hets : {len(mut_hets)} ")
        else:
            pass ## do nothing for non hets 
        
        ## Annotation
        candidate.samples[sample_id]["WT_HETS"] = len(wt_hets)
        candidate.samples[sample_id]["MUT_HETS"] = len(mut_hets)
        
        if len(wt_hets) + len(mut_hets) > high_het:
            classification = "HIGH_HET"
        else:
            classification = "PASS"
        candidate.samples[sample_id]["SVHET"] = classification
        
        ## Write to VCF
        out_vcf.write(candidate)
    
    out_vcf.close()
    wt.close()
    mut.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter deletion candidates based on short variant callsets")
    parser.add_argument("--candidate-vcf", required=True, help="Path to sample-level candidate VCF")
    parser.add_argument("--wt-vcf", required=True, help="Path to sample-level WT short variant VCF")
    parser.add_argument("--mut-vcf", required=True, help="Path to sample-level MUT short variant VCF")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--out", required=True, help="Output path with file name")
    parser.add_argument("--min-dp", default=5, type=int, help="Minimum alternate DP to consider a HET reliable")
    parser.add_argument("--high-het", default=1, type=int, help="Minimum HET count to reject a DEL")

    args = parser.parse_args()

    ## Test command:
    ## ./04_het_evaluator.py --candidate-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/small_candidates.vcf --wt-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/HG00096-wt-reads-small.bcf --mut-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/HG00096-mut-reads-small.bcf --sample-id HG00096
    ## ./04_het_evaluator.py --candidate-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/large_candidates.vcf --wt-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/HG00096-wt-reads-large.bcf --mut-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/HG00096-mut-reads-large.bcf --sample-id HG00096
    
    filter_main(args.candidate_vcf, args.wt_vcf, args.mut_vcf, args.min_dp, args.high_het, args.sample_id, args.out)

