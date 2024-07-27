#!/usr/bin/env python3

# This module identifies true homozygosities by homozygous-supporting reads with heterozygous anchors at 5' and 3' breakpoints.

import pysam
from collections import deque

from utils import trim_cigar_clipped_with_seq, read_tokenizer, left_align_to_aligned_start

def get_heterozygous_sites(read_tokens, breakpoint = "5prime") -> list[int]:

    het = []
    ## Initialize an empty dict
    read_tokens.sort(key=lambda x: x[-1][0])
    if breakpoint == "5prime":
        start, end = read_tokens[0][-1][0], read_tokens[0][0][0]
    elif breakpoint == "3prime":
        start, end = read_tokens[0][0][0], read_tokens[-1][-1][0]

    nt_dict = {}
    for pos in range(start, end + 1, 1):
        nt_dict[pos] = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0, "I": 0}

    ## Populate with nt count
    for read_token in read_tokens:
        for pos, nt in read_token:
            if nt == "N":
                continue
            nt_dict[pos][nt] += 1
    
    ## Identify het sites
    for pos, counts in nt_dict.items():
        # Check if there are multiple non-zero counts / below-threshold counts
        if breakpoint == "5prime" and pos == end: # breakpoint uncertainty
            continue
        elif breakpoint == "3prime" and pos == start:
            continue

        allele = {}
        for nt, count in counts.items():
            if count > 5:
                allele[nt] = count
        if "I" in allele.keys() and len(allele.keys()) > 2:   
            het.append(pos)
        elif "I" not in allele.keys() and len(allele.keys()) == 2:
            het.append(pos)

    return het

def get_heterozygous_sites_from_bam(bam_fp, region, breakpoint, repeat_regions):
    
    import copy
    
    with pysam.AlignmentFile(bam_fp) as bam:
        read_tokens = []
        
        chrom = region.chrom
        start, end = region.start, region.end
        
        readfs = bam.fetch(str(chrom), end - 2, end - 1) if breakpoint == "3prime" else bam.fetch(str(chrom), start + 1, start + 2)
        
        for read in readfs:
            if read.is_unmapped or read.mapping_quality <= 20:
                continue
            if read.get_tag("NM") > 5:
                continue
            # discordant pairs
            if abs(read.template_length) > 1000:
                continue
            # first trim soft-clipped/hard-clipped bases on both sides
            aligned_pos = read.reference_start + 1 # we use 1-based
            cigar_block = read.cigartuples
            seq, cigar = trim_cigar_clipped_with_seq(read.seq, cigar_block)
            read_token = read_tokenizer(seq, cigar, aligned_pos, token=deque())
            
            if breakpoint == "5prime":
                if (start - 1, "D") in read_token or (start, "D") in read_token:
                    continue
            elif breakpoint == "3prime":
                if (end, "D") in read_token or (end + 1, "D") in read_token:
                    continue
            
            # store the refined reads with cigar as a deque (position, nt)
            read_token = deque(read_token)
            if breakpoint == "5prime":
                read_token.reverse()
            read_tokens.append(read_token)
    
    if len(read_tokens) == 0:
        return []
    
    # left-align refined reads to aligned_pos
    try:
        if breakpoint == "3prime":
            left_aligned_reads = list(map(lambda x: left_align_to_aligned_start(x, end, breakpoint = breakpoint), read_tokens))
        elif breakpoint == "5prime":
            left_aligned_reads = list(map(lambda x: left_align_to_aligned_start(x, start, breakpoint = breakpoint), read_tokens))
    except Exception as e:
        print(f"{chrom}:{start}-{end}")
        raise ValueError(e)
    # clipped nts may shift the aligned start position
    left_aligned_reads = list(filter(lambda x: x is not None, left_aligned_reads))
    
    if len(left_aligned_reads) > 0:
        het_sites = get_heterozygous_sites(left_aligned_reads, breakpoint = breakpoint)
    else:
        het_sites = []
    
    # Exclude HET sites within 5bp up/downstream of TRs <= 50bp
    with pysam.TabixFile(repeat_regions, parser=pysam.asTuple()) as repeats:
        search_start, search_end = (start - 100, start) if breakpoint == "5prime" else (end, end + 100)
        true_hets = het_sites.copy()
        for repeat in repeats.fetch(chrom, search_start, search_end):
            for het in het_sites:
                if het >= int(repeat[1]) and het <= int(repeat[2]):
                    true_hets.remove(het)
    
    return true_hets