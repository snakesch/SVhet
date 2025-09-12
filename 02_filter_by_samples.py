#!/usr/bin/env python3

# 2. For each deletion candidate (small/large), append a new FORMAT field named SVHET (values can be PASS, MUT, WT or HET).
# 2. For each sample in small/large candidate VCF, iterate over all deletion candidates:
# a. if the candidate is heterozygous and not uncalled in the sample, then proceed to step 2b. If not, the current sample should be assigned "PASS" under FORMAT/SVHET (Not filtered).
# b. For each heterozygous deletion, fetch reads corresponding to wt reads and mut reads (see later scripts) and write the two read sets to two separate BAM files in /tmp/...
# c. Call variants using the two BAM files separately.
# d. Extract heterozygous sites within candidate deletions with predefined ad2dp and min_dp (see later scripts).
# e. Determine if the deletion candidate should be rejected on the basis of excess heterozygosity observed in WT reads (FORMAT/SVHET = WT) or MUT reads (FORMAT/SVHET = MUT) or within deletion region (FORMAT/SVHET = HET). This will have to be evaluated under a Gaussian model (needs further polishing to enhance model stability).
# 3. Combine all results and the annotated VCFs. Merge the results with the non-candidate VCF and writes to final VCF output. Tabix-index final VCF output.
# The original script did not consider small and large candidates as separate files. So I'm re-implementing stuff.

## This script should be run in two modes: small and large candidates

## small samples:
## python 02_filter_by_samples.py --sample HG00096 --bam /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams/HG00096.sorted.markdup.bam --cohort-vcf /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/small_candidates.vcf --out /lustre1/g/paed_yangwl/snakesch/tmp/test_svhet --mode small

## - Constants - ##
MAPQ = 30 ## min mapping quality for read evidence

import os
import sys
import itertools
import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import argparse

import pysam

## - read-specific helpers - ##
def get_insert_dist(bam_file, num_runs=3, num_samples=2000, max_template_length=8000, sample_prob=0.001):
    """
    Runs insert size estimation multiple times, averages statistics, and handles filtering.

    Args:
        bam_file (str): Path to the BAM file.
        num_runs (int): Number of times to run the estimation.
        num_samples (int): Number of samples per run.
        max_template_length (int): Maximum template length.
        sample_prob (float): Sampling probability.

    Returns:
        tuple: (Average Mean, Average Median, Average Std), or (None, None, None) if no data.
    """
    import random

    all_stats = []
    for _ in range(num_runs):
        insert_sizes = []

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if (read.flag & 0x2 and not (read.flag & 0x90C) and
                    read.reference_id == read.next_reference_id and
                    abs(read.template_length) <= max_template_length and
                    random.random() < sample_prob):

                    insert_sizes.append(abs(read.template_length))
                    if len(insert_sizes) >= num_samples:
                        break  # Stop when num_samples reached

        if insert_sizes:
            threshold = np.percentile(insert_sizes, 99)
            filtered_sizes = [size for size in insert_sizes if size <= threshold][:num_samples]
            if filtered_sizes:
                stats = np.array([np.mean(filtered_sizes), np.median(filtered_sizes), np.std(filtered_sizes)])
                all_stats.append(stats)

    return tuple(np.mean(np.stack(all_stats), axis=0)) if all_stats else (None, None, None)

def is_clipped(read):
    '''Determine if a read is (soft/hard) clipped on either end'''
    first_op, last_op = read.cigartuples[0], read.cigartuples[-1]
    return first_op in (4, 5) or last_op in (4, 5)

def generate_variants(cohort_vcf, sample_id):

    vcf = pysam.VariantFile(cohort_vcf)
    vcf.subset_samples([sample_id])

    variant_pool = {}
    for record in vcf.fetch():
        variant_key = (
            record.chrom,
            (record.start, record.stop),
            record.ref,
            record.alts[0],
            record.info["CIPOS"],
            record.info["CIEND"]
        )
        gt = record.samples[sample_id]["GT"]
        result = "not_analysed" if gt == (0, 1) else "PASS"
        variant_pool[variant_key] = (gt, result)
    vcf.close()

    return variant_pool

def extract_wt_reads_per_candidate(candidate, bam, isize_mean, isize_std, mode="small"):
    '''
    Candidate variant should be formatted as :
    {
        ('chr1', (216567, 223202316), 'N', '<DEL>', (-5, 4), (-10, 4)), ...
    }

    Arguments:
        candidate: as above
        bam: pysam.AlignmentFile
        isize_mean, isize_std: insert size average and std

    Returns:
        list of reads
    '''
    chrom, coords, _, _, start_ci, end_ci = candidate

    rg1, rg2, rg3 = [], [], []

    ## check for duplicate entries
    read_id = set()

    ## Reads mapped to 5' breakpoint
    for read in bam.fetch(chrom, coords[0] - abs(start_ci[0]), coords[0] + start_ci[1]):
        if read.mapping_quality < MAPQ or read.is_duplicate or is_clipped(read) or read.template_length == 0: ## last case may be a TRA ?
            continue
        
        if read.has_tag("XA") or read.has_tag("SA"): ## skip if is part of a split alignment / multiple alignment
            continue
        
        if read.get_tag("NM") > 5: ## Unlikely originated from true haplotype
            continue
        
        if abs(abs(read.template_length) - isize_mean) < isize_std:
            if (read.query_name, read.flag) not in read_id:
                rg1.append(read)
                read_id.add((read.query_name, read.flag))

    ## Reads mapped to 3' breakpoint
    for read in bam.fetch(chrom, coords[1] - abs(end_ci[0]), coords[1] + end_ci[1]):
        if read.mapping_quality < MAPQ or read.is_duplicate or is_clipped(read) or read.template_length == 0:
            continue
        
        if read.has_tag("XA") or read.has_tag("SA"): ## skip if is part of a split alignment / multiple alignment
            continue
        
        if read.get_tag("NM") > 5: ## Unlikely originated from true haplotype
            continue
        
        if abs(abs(read.template_length) - isize_mean) < isize_std:
            if (read.query_name, read.flag) not in read_id:
                rg2.append(read)
                read_id.add((read.query_name, read.flag))

    ## For larger deletions
    if mode == "small": ## fetch intra-deletion reads for small candidates
        if coords[0] + abs(start_ci[1]) < coords[1] - abs(end_ci[0]):
            for read in bam.fetch(chrom, coords[0] + abs(start_ci[1]), coords[1] - abs(end_ci[0])):
                if read.mapping_quality < MAPQ or read.is_duplicate or is_clipped(read):
                    continue
                
                if read.has_tag("XA") or read.has_tag("SA"):
                    continue
                
                if read.get_tag("NM") > 5: ## Unlikely originated from true haplotype
                    continue
                
                if read.reference_start > coords[0] + abs(start_ci[1]) and read.reference_end < coords[1] - abs(end_ci[0]):
                    if (read.query_name, read.flag) not in read_id:
                        rg3.append(read)
                        read_id.add((read.query_name, read.flag))

    # if len(rg1) + len(rg2) + len(rg3) > 50:
    #     print(f"Processing candidate : {candidate}")
    #     print(f"Read group 1: {len(rg1)}")
    #     print(f"Read group 2: {len(rg2)}")
    #     print(f"Read group 3: {len(rg3)}")

    # for read in rg1:
    #     print(read.query_name, read.template_length)

    return rg1 + rg2 + rg3 ## all read groups should be disjoint here

def extract_mut_reads_per_candidate(candidate, bam, isize_mean, isize_std, mode):
    '''
    Candidate variant should be formatted as :
    {
        ('chr1', (216567, 223202316), 'N', '<DEL>', (-5, 4), (-10, 4)), ...
    }

    Returns:
        list of reads
    '''
    chrom, coords, _, _, start_ci, end_ci = candidate

    read_id = set()
    mut_reads = []

    if mode == "small":
        regions = [
            (chrom, coords[0] - abs(start_ci[0]), coords[1] + end_ci[1])
        ]
    elif mode == "large":
        regions = [
            (chrom, coords[0] - abs(start_ci[0]), coords[0] + start_ci[1]),
            (chrom, coords[1] - abs(end_ci[0]), coords[1] + end_ci[1])
        ]
    merged_reads = itertools.chain(*(bam.fetch(*region) for region in regions))

    for read in merged_reads:
        if read.mapping_quality < MAPQ or read.is_duplicate:
            continue

        ## Fetch DP and SA reads
        if abs(read.template_length) >= isize_mean + isize_std:
            if (read.query_name, read.flag) not in read_id:
                mut_reads.append(read)
                read_id.add((read.query_name, read.flag))
        elif read.has_tag("SA"):
            sa = read.get_tag("SA")
            ## SA reads should be mapped to the same chrom and split interval size similar to svlen
            sa_chrom, sa_pos, *_ = sa.split(",")
            if sa_chrom == chrom and abs(int(sa_pos) - read.reference_end) - abs(coords[1] - coords[0]) < 500:
                if (read.query_name, read.flag) not in read_id:
                    mut_reads.append(read)
                    read_id.add((read.query_name, read.flag))

    return mut_reads

def process_candidate(args):
    variant_key, bam_fp, isize_mean, isize_std, mode = args
    bam = pysam.AlignmentFile(bam_fp, "rb")
    wt_reads = extract_wt_reads_per_candidate(variant_key, bam, isize_mean, isize_std, mode)
    mut_reads = extract_mut_reads_per_candidate(variant_key, bam, isize_mean, isize_std, mode)
    bam.close()
    return [read.to_string() for read in wt_reads], [read.to_string() for read in mut_reads]

def generate_read_evidences_per_sample(cohort_vcf, sample_id, bam_fp, outdir, n_workers=4, mode="small"):
    '''
    Main function for read evidence generation.

    Arguments:
        cohort_vcf: candidate SV VCF
        sample_id: sample ID to analyse for the current script
        bam_fp: complete sample BAM file
        outdir: scratch directory

    Writes to file:
        wt_reads, mut_reads: read evidences collected for the SV candidate
    '''
    variant_pool = generate_variants(cohort_vcf=cohort_vcf, sample_id=sample_id)
    wt_bam_path = os.path.join(outdir, sample_id + f"-wt-reads-{mode}.tmp.bam")
    mut_bam_path = os.path.join(outdir, sample_id + f"-mut-reads-{mode}.tmp.bam")
    isize_mean, _, isize_std = get_insert_dist(bam_fp)
    print(f"Estimated insert size distribution: mean {isize_mean:.2f}bp, std {isize_std:.2f}bp")

    # Prepare tasks (filtering not_analysed upfront)
    tasks = [
        (variant_key, bam_fp, isize_mean, isize_std, mode)
        for variant_key, (gt, status) in variant_pool.items()
        if status == "not_analysed"
    ]

    bam = pysam.AlignmentFile(bam_fp)
    wt_bam = pysam.AlignmentFile(wt_bam_path, "wb", header=bam.header)
    mut_bam = pysam.AlignmentFile(mut_bam_path, "wb", header=bam.header)

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for wt_reads, mut_reads in tqdm(executor.map(process_candidate, tasks), desc="Collecting evidence", total=len(tasks), file=sys.stdout):
            for read in wt_reads:
                read_obj = pysam.AlignedSegment.fromstring(read, bam.header)
                wt_bam.write(read_obj)
            for read in mut_reads:
                read_obj = pysam.AlignedSegment.fromstring(read, bam.header)
                mut_bam.write(read_obj)
    wt_bam.close()
    mut_bam.close()
    bam.close()

    ## Sort BAM files and create index
    # print("Sorting alignments ... ")
    final_wt_bam_path = wt_bam_path.removesuffix(".tmp.bam") + ".bam"
    final_mut_bam_path = mut_bam_path.removesuffix(".tmp.bam") + ".bam"
    pysam.samtools.sort("-o", final_wt_bam_path, wt_bam_path, catch_stdout=False)
    pysam.samtools.sort("-o", final_mut_bam_path, mut_bam_path, catch_stdout=False)
    pysam.index(final_wt_bam_path)
    pysam.index(final_mut_bam_path)

    ## Cleanup
    os.remove(wt_bam_path)
    os.remove(mut_bam_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read evidence generator")
    parser.add_argument("--cohort-vcf", required=True, help="Path to cohort VCF")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--bam", required=True, help="Path to sample BAM file (must be indexed)")
    parser.add_argument("--outdir", required=True, help="Scratch directory")
    parser.add_argument("--jobs", default=4, type=int, help="Number of CPU workers")
    parser.add_argument("--mode", required=True, help="Operation mode (small vs large candidates)")

    args = parser.parse_args()

    generate_read_evidences_per_sample(args.cohort_vcf, args.sample, args.bam, args.outdir, n_workers=args.jobs, mode=args.mode)

