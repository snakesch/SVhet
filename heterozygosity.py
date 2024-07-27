#!/usr/bin/env python3

# This module identifies highly confident subsets of heterozygous sites that suggest absence of heterzygous deletion.

from pysam import VariantFile

def filter_by_presence_of_het(callset, het_sites: VariantFile, difficult_regions: str):
    
    import copy
    
    from utils import query_region_overlap_with_bed
    
    filtered_callset = []
    for _call in callset:
        call = copy.deepcopy(_call)
        for rec in het_sites.fetch(call.chrom, call.start, call.end):
            
            if rec.pos <= call.start or rec.pos >= call.end:
                continue
            
            ## determine true het and reject candidate het deletion
            variant_ad = rec.samples[0]["AD"]
            sample_gt = rec.samples[0]["GT"]
            call_quality = rec.samples[0]["GQ"]
            major_ad = max(variant_ad[sample_gt[0]], variant_ad[sample_gt[1]])
            minor_ad = _ad = min(variant_ad[sample_gt[0]], variant_ad[sample_gt[1]])
            
            if minor_ad < 5 or call_quality <= 10:
                continue
            
            if minor_ad / major_ad > .25 :
                call.het_cnt += 1
                
        # number of nucleotides overlapped with difficult regions
        total_span = call.end - call.start
        
        valid_span = total_span
        
        if valid_span < 100 and call.het_cnt > 1:
            continue
        elif valid_span < 1000 and call.het_cnt > 2:
            continue
        elif valid_span >= 1000:
            expected_heterozygosity = 10 # min(max(10, np.ceil(0.001 * (call.end - call.start))), 20)
            if call.het_cnt > expected_heterozygosity:
                continue
        
        filtered_callset.append(call)

    return filtered_callset