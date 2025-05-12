"""Module for gathering read evidence for SVs."""

import os
import pysam
import cyvcf2
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import pybedtools as pb
from functools import partial

from svhet.const import contigs
from svhet.sv import SV
from svhet.variants import is_het_deletion
from svhet.utils.log import setup_logger

logger = setup_logger("svhet", level="DEBUG")

def _process_chunk(variants_chunk, bam_fp, header_dict, read_length, 
                  cipos_tag, ciend_tag, vcf_header):
    regions = []
        
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.vcf') as tmp_vcf, \
         tempfile.NamedTemporaryFile(delete=False, suffix='.wt.bam') as wt_tmp, \
         tempfile.NamedTemporaryFile(delete=False, suffix='.mut.bam') as mut_tmp:

        tmp_vcf.write(vcf_header)
        tmp_vcf.writelines(variants_chunk)
        tmp_vcf.flush()
        
        with pysam.AlignmentFile(bam_fp, "rb", threads=8) as bam, \
             pysam.AlignmentFile(wt_tmp.name, "wb", header=header_dict, threads=4) as wt_bam, \
             pysam.AlignmentFile(mut_tmp.name, "wb", header=header_dict, threads=4) as mut_bam:

            vcf = cyvcf2.VCF(tmp_vcf.name)
            for variant in vcf():
                if variant.CHROM not in contigs:
                    continue
            
                if is_het_deletion(variant):
                    sv = SV(variant, read_length=read_length, 
                           cipos_tag=cipos_tag, ciend_tag=ciend_tag)

                    if sv.ambiguous_start or sv.ambiguous_end:
                        continue
                    
                    wt_reads = sv.get_wt_haplotype_evidence(bam, mapping_quality=30)
                    mut_reads = sv.get_mut_haplotype_evidence(bam, mapping_quality=30)
                    # print(f"Processing {variant.ID} with {len(wt_reads)} wt reads and {len(mut_reads)} mut reads")
                    for r in wt_reads:
                        wt_bam.write(r.read)
                    for r in mut_reads:
                        mut_bam.write(r.read)
                    
                    if wt_reads or mut_reads:
                        regions.append(f"{sv.variant.CHROM}\t{sv.start_ci[0] - read_length * 10}\t{sv.end_ci[1] + read_length * 10}")
            vcf.close()
    
    pysam.sort(wt_tmp.name, "-o", wt_tmp.name)
    pysam.sort(mut_tmp.name, "-o", mut_tmp.name)
            
    pysam.index(wt_tmp.name)
    pysam.index(mut_tmp.name)
    # print(f"Temporary files created: {wt_tmp.name}, {mut_tmp.name}")
    return {'wt_bam': wt_tmp.name, 'mut_bam': mut_tmp.name, 'regions': regions}


def gather_read_evidence(vcf_fp, bam_fp, wt_evidence_bam_fp, mut_evidence_bam_fp, 
                        candidate_regions, read_length, cipos_tag, ciend_tag, 
                        max_workers=4):
    """
    Gather read evidence for wild-type and mutant haplotypes.
    
    Args:
        vcf_fp: Path to VCF file with structural variants
        bam_fp: Path to BAM file
        wt_evidence_bam_fp: Output path for wild-type evidence BAM
        mut_evidence_bam_fp: Output path for mutant evidence BAM
        candidate_regions: Output path for candidate regions BED
        read_length: Average read length from the sequencing data
        cipos_tag: INFO field tag for CIPOS (default: CIPOS)
        ciend_tag: INFO field tag for CIEND (default: CIEND)
        max_workers: Number of threads to use for processing (default: 8)
    """
    
    with tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix='.candidate.vcf') as candidate_vcf:
        
        vcf = cyvcf2.VCF(vcf_fp)
        vcf_header = vcf.raw_header
        candidate_vcf.write(vcf_header)
        
        variant_pos = [v.POS for v in vcf()]
        vcf.close()
        
        vcf = cyvcf2.VCF(vcf_fp)
        variants = []
        for variant in vcf():
            if is_het_deletion(variant) and variant.CHROM in contigs:
                if variant.end - variant.POS > 1e6 and any(variant.POS < ovp < variant.end for ovp in variant_pos): ## complex SVs
                    continue
                elif variant.end - variant.POS > 1e6:
                    logger.warning(f"Detected large variant: {variant.ID} (size: {variant.end - variant.POS}). Processing may be slow.")
                candidate_vcf.write(str(variant))
                variants.append(str(variant))
        vcf.close()

    with pysam.AlignmentFile(bam_fp, "rb") as bam:
        header_dict = bam.header.to_dict()
        if "RG" not in header_dict:
            header_dict["RG"] = [{'ID': "SAMPLE1", 'SM': "SAMPLE1", 'LB': "LIB1", 'PL': "PL1"}]

    chunk_size = max(100, len(variants) // (max_workers * 4))  # Dynamic batching
    logger.debug(f"Chunk size: {chunk_size}")
    
    chunks = [variants[i:i+chunk_size] for i in range(0, len(variants), chunk_size)]
 
    worker = partial(_process_chunk,
                    bam_fp=bam_fp,
                    header_dict=header_dict,
                    read_length=read_length,
                    cipos_tag=cipos_tag,
                    ciend_tag=ciend_tag,
                    vcf_header=vcf_header)

    results = []
    temp_files = []
    futures = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for idx, chunk in enumerate(chunks):
            futures[executor.submit(worker, chunk)] = chunk
        
        with tqdm(total=len(variants), desc="Gathering read evidence") as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                    temp_files.extend([result['wt_bam'], result['mut_bam'], result['wt_bam'] + '.bai', result['mut_bam'] + '.bai'])
                    pbar.update(len(futures[future]))  # Update by chunk size
                except Exception as e:
                    logger.error(f"Error processing chunk: {e}")

    with ThreadPoolExecutor(2) as merge_executor:
        merge_jobs = []
        for bam_type, out_fp in [('wt_bam', wt_evidence_bam_fp), ('mut_bam', mut_evidence_bam_fp)]:
            inputs = [r[bam_type] for r in results if os.path.exists(r[bam_type])]
            if inputs:
                merge_jobs.append(merge_executor.submit(
                    pysam.merge, "-@", str(max_workers - 1), "-f", "-c", "-h", inputs[0], out_fp, *inputs
                ))
        
        for job in as_completed(merge_jobs):
            job.result()

    all_regions = []
    for r in results:
        all_regions.extend(r.get('regions', []))
    
    if all_regions:
        pb.BedTool("\n".join(all_regions), from_string=True).sort().merge().saveas(candidate_regions)
    else:
        open(candidate_regions, 'w').close()

    for f in temp_files:
        try:
            os.remove(f)
        except Exception as e:
            pass

