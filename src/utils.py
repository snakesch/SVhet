#!/usr/bin/env python3

from collections import deque
import pandas as pd

import logging
logger = logging.getLogger("root")

def trim_cigar_clipped_with_seq(seq, cigar):
    ## here cigar is the cigartuples from pysam.AlignedSegment
    start_cigar = cigar[0][0]
    end_cigar = cigar[-1][0]
    is_clipped = False
    if start_cigar == 4 or start_cigar == 5:
        seq = seq[cigar[0][1]: ]
        cigar = cigar[1:]
        is_clipped = True
    if end_cigar == 4 or end_cigar == 5:
        seq = seq[: -cigar[-1][1]]
        cigar = cigar[:-1]
        is_clipped = True
    if not is_clipped:
        return (seq, cigar)
    else:
        return trim_cigar_clipped_with_seq(seq, cigar)

def read_tokenizer(seq, cigar, pos, token=deque()):
    
    if len(cigar) == 0:
        return token
    
    cigar_operation = cigar[0][0]
    offset = cigar[0][1]
    
    if cigar_operation == 0: # Case: "M"
        add_deque = deque(enumerate(seq[: offset], pos))
        pos = pos + offset
        seq = seq[offset: ]
    elif cigar_operation == 1: # Case: "I"
        add_deque = deque([(pos, "I")])
        seq = seq[offset: ]
    elif cigar_operation == 2: # Case: "D"
        add_deque = deque(enumerate("D" * offset, pos)) 
        pos = pos + offset
    else:
        raise ValueError(f"Unhandled CIGAR operation {cigar_operation}")
    
    token.extend(add_deque)
    
    if len(cigar) > 1:
        return read_tokenizer(seq, cigar[1:], pos, token=token)
    else:
        return token

def left_align_to_aligned_start(read_token, pos, breakpoint = "3prime"):

    if len(read_token) == 0:
        return None
    
    if read_token[0][0] == pos:
        return read_token
    
    if breakpoint == "3prime":
        if not read_token[0][0] < pos:
            # logger.warning(f"Identified reads with aligned_pos {read_token[0][0]} larger than expected {pos}. ")
            return None
    elif breakpoint == "5prime":
        if not read_token[0][0] > pos:
            # logger.warning(f"Identified reads with aligned_pos {read_token[0][0]} smaller than expected {pos}. ")
            return None
    else:
        raise ValueError(f"Unhandled direction - {breakpoint}")
    
    read_token.popleft()
    return left_align_to_aligned_start(read_token, pos, breakpoint = breakpoint)

def loadVCF(path, omit_record=False, encoding=None, resolve_info=True):
    """
    Input a VCF file and returns:
    1. header (lines begin with #)
    2. subject ID list
    3. record dataframe
    """
    import gzip
    header = []
    subjects = []
    is_gzip = path.endswith(".gz")
    
    # Function to process lines
    def process_line(line):
        nonlocal subjects
        if line.startswith("#"):
            header.append(line.strip())
            if line.startswith("##"):
                return
            subjects = line.strip().split("\t")[9:]
    
    # Read file
    with (gzip.open(path, "rt", encoding=encoding) if is_gzip else open(path, "rt", encoding=encoding)) as f:
        for line in f:
            process_line(line)
            if subjects:
                break

    if omit_record:
        return header, subjects

    vcf = pd.read_csv(path, sep="\t", na_filter=False, engine="c", comment="#", header=None, 
                     compression="gzip" if is_gzip else None, encoding=encoding)
    
    vcf.columns = header[-1].split("\t")
    
    if not header or not subjects:
        raise RuntimeError("Incorrect VCF format.")
    
    if resolve_info:
        info_fields = [h.split(",")[0][11:] for h in header if "INFO" in h and "##bcftools" not in h]
        vcf.loc[:, info_fields] = "."
        
        # Populate the columns
        for idx, row in vcf.iterrows():
            for f in row["INFO"].split(";"):
                if "=" in f:
                    key, value = f.split("=")
                    if key in info_fields:
                        row[key] = value
                elif f in info_fields:
                    row[f] = True
            vcf.iloc[idx, :] = row
        
        vcf = vcf.drop("INFO", axis=1)
    
    return header, subjects, vcf

def create_variants(vcf):
    
    from src.variant import StructuralVariant
    
    from src.constants import VALID_CONTIGS
    
    header, sample, sv_calls = loadVCF(vcf)
    
    sv_callset = [ StructuralVariant(row) for _, row in sv_calls.iterrows() ]
    
    return sv_callset

def query_region_overlap_with_bed(chrom, start, end, bed):

    import pybedtools as pb
    
    bed = pd.read_csv(bed, sep="\t", header=None)
    reference_regions = pb.BedTool.from_dataframe(bed)
    query = pb.BedTool(f"{chrom}\t{start}\t{end}", from_string=True)
    nt_overlap = reference_regions.intersect(query).count()
    
    return nt_overlap

def write_to_bed(regions, outd, outp, **kwargs) -> None:
    """
    Write a list of genomic regions to a BED file.

    Args:
        regions (list): A list of StructuralVariant objects.
        outd (str): The output directory where the BED file will be written.
        outp (str): The filename of the output BED file.
        **kwargs: Additional keyword arguments to pass to the pandas.DataFrame.to_csv() function, such as index=False, header=False.

    Returns:
        None
    """
    import os

    # Create a pandas DataFrame from the list of genomic regions
    df = pd.DataFrame.from_records([(r.chrom, r.start, r.end) for r in regions])

    # Write the DataFrame to a BED file in the specified output directory
    df.to_csv(os.path.join(outd, outp), **kwargs)