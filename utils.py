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

def loadVCF(path, omit_record=False, encoding = None, resolve_info=True):
    """
    Input a VCF file and returns:

    1. header (lines begin with #)
    2. subject ID list
    3. record dataframe

    """
    from gzip import GzipFile

    header = []
    subjects = []
    
    if path.endswith(".gz"):
        with GzipFile(path) as f:
            line = next(f)
            while line:
                if line.decode()[1] == "#":
                    header.append(line.decode().strip())
                elif line.decode()[0] == "#":
                    col_names = line.decode().strip().split("\t")
                    subjects = col_names[9:]
                else:
                    break
                try:
                    line = next(f)
                except StopIteration:
                    break
    else:
        with open(path, "rt", encoding=encoding) as f:
            line = f.readline()
            while line:
                if line[0] == "#":
                    header.append(line.strip())
                elif line[0] != "#":
                    col_names = line.strip().split("\t")
                    subjects = col_names[9:]
                    break
                line = f.readline()

    if omit_record:
        return header, subjects

    vcf = pd.read_csv(path, sep="\t", na_filter=False, engine="c", comment="#", header=None, compression="gzip", encoding = encoding)
    vcf.columns = col_names
    if header == None or subjects == None:
        raise RuntimeError("Incorrect VCF format.")
        
    if resolve_info:
        # Create new columns
        info_fields = [h.split(",")[0][11:] for h in header if "INFO" in h and "##bcftools" not in h]
        vcf.loc[:, info_fields] = "."

        # Populate the columns
        for idx, c in vcf.iterrows():
            for f in c["INFO"].split(";"):
                # if we have a key-value pair
                if "=" in f:
                    key, value = f.split("=")
                    if key in info_fields:
                        c[key] = value
                # if not key-value pair, just tags
                else:
                    if f in info_fields:
                        c[f] = True

            # Assign the updated row back to the DataFrame
            vcf.iloc[idx, :] = c
        vcf = vcf.drop("INFO", axis=1)
    
    return header, subjects, vcf

def create_variants(vcf):
    
    from variant import StructuralVariant
    
    from constants import VALID_CONTIGS
    
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