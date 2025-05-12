"""Main filtering module for SVhet."""

import os
import pysam
import cyvcf2

from svhet.const import contigs
from svhet.evidence import gather_read_evidence
from svhet.variants import extract_het_positions, extract_variants_within_regions, call_short_variants, is_het_deletion
from svhet.sv import SV
from svhet.utils.bam import get_read_length
from svhet.utils.log import setup_logger

logger = setup_logger("svhet")

def filter_sv_callset(vcf_fp, bam_fp, outdir, reference_fasta_path, 
                     output_filename="filtered.vcf.gz", threads=4, 
                     cipos_tag="CIPOS", ciend_tag="CIEND",
                     ad2dp=0.25, min_dp=10, fully_within=True,
                     image_path=None, region_file=None, force_rerun=False):
    """
    Main function to filter structural variants based on heterozygous SNPs.
    
    Args:
        vcf_fp: Path to input VCF file with structural variants
        bam_fp: Path to input BAM file
        outdir: Output directory
        reference_fasta_path: Path to reference genome FASTA
        output_filename: Name of output VCF file
        threads: Number of CPU threads to use
        cipos_tag: INFO field tag to search for CIPOS (default: CIPOS)
        ciend_tag: INFO field tag to search for CIEND (default: CIEND)
        image_path: Path to DeepVariant singularity image
        region_file: BED file with regions to analyze
        fully_within: If True, only consider variants fully within regions. Otherwise, consider variants partially within regions.
        force_rerun: Force reprocessing even if output files exist
    
    Returns:
        Path to filtered VCF file
    """
    
    logger.info("Starting filter_sv_callset with parameters:")
    logger.info(f"  vcf_fp: {vcf_fp}")
    logger.info(f"  bam_fp: {bam_fp}")
    logger.info(f"  outdir: {outdir}")
    logger.info(f"  reference_fasta_path: {reference_fasta_path}")
    logger.info(f"  output_filename: {output_filename}")
    logger.info(f"  threads: {threads}")
    # logger.info(f"  ad2dp: {ad2dp}")
    # logger.info(f"  min_dp: {min_dp}")
    logger.info(f"  cipos_tag: {cipos_tag}")
    logger.info(f"  ciend_tag: {ciend_tag}")
    logger.info(f"  image_path: {image_path}")
    logger.info(f"  region_file: {region_file}")
    logger.info(f"  fully_within: {fully_within}")
    logger.info(f"  force_rerun: {force_rerun}")
    
    # Define output file paths
    wt_evidence_bam_fp = f"{outdir}/wt_evidence.bam"
    sorted_wt_evidence_fp = f"{outdir}/wt_evidence.sorted.bam"
    mut_evidence_bam_fp = f"{outdir}/mut_evidence.bam"
    sorted_mut_evidence_fp = f"{outdir}/mut_evidence.sorted.bam"
    candidate_regions = f"{outdir}/candidates.bed"
    wt_vcf = f"{outdir}/wt_evidence.vcf.gz"
    mut_vcf = f"{outdir}/mut_evidence.vcf.gz"
    out_filtered_vcf = f"{outdir}/{output_filename}"
    
    # Check read length
    read_length = get_read_length(bam_fp)
    if read_length > 500:
        raise NotImplementedError("SVhet does not support long reads.")
    
    # Process region file if specified
    if region_file is not None and region_file != "":
        sub_vcf_fp = _process_region_file(region_file, vcf_fp, outdir, fully_within=fully_within)
    else:
        sub_vcf_fp = vcf_fp
        
    # Gather read evidence
    if (not os.path.exists(sorted_wt_evidence_fp) or 
        not os.path.exists(sorted_mut_evidence_fp) or 
        not os.path.exists(candidate_regions) or 
        force_rerun):
        
        logger.info(f"Gathering read evidence for {sub_vcf_fp}")
        gather_read_evidence(
            vcf_fp=sub_vcf_fp,
            bam_fp=bam_fp,
            wt_evidence_bam_fp=wt_evidence_bam_fp,
            mut_evidence_bam_fp=mut_evidence_bam_fp,
            candidate_regions=candidate_regions,
            read_length=read_length,
            cipos_tag=cipos_tag,
            ciend_tag=ciend_tag,
        )
        
        # Sort and index BAM files
        pysam.sort("-o", sorted_wt_evidence_fp, wt_evidence_bam_fp)
        pysam.sort("-o", sorted_mut_evidence_fp, mut_evidence_bam_fp)
        pysam.index(sorted_wt_evidence_fp)
        pysam.index(sorted_mut_evidence_fp)

        # Clean up unsorted BAM files
        os.remove(wt_evidence_bam_fp)
        os.remove(mut_evidence_bam_fp)
    
    # Call short variants if needed
    if force_rerun or not (os.path.exists(wt_vcf) and os.path.exists(mut_vcf)):
        call_short_variants(
            sorted_wt_evidence_fp, 
            sorted_mut_evidence_fp, 
            candidate_regions, 
            image_path, 
            reference_fasta_path, 
            threads
        )
    
    # Extract heterozygous positions
    wt_pos, mut_pos = extract_het_positions(wt_vcf, mut_vcf, ad2dp=ad2dp, min_dp=min_dp)
        
    filter_variants(vcf_fp, out_filtered_vcf, wt_pos, mut_pos, 
                    read_length, cipos_tag=cipos_tag, ciend_tag=ciend_tag)
    
    return out_filtered_vcf

def _process_region_file(region_file, vcf_fp, outdir, fully_within=True):
    """Process region file and extract variants within regions."""
    import pybedtools as pb
    if not os.path.exists(region_file):
        raise FileNotFoundError("Target region file not found.")
    
    if not os.path.isfile(vcf_fp):
        raise FileNotFoundError(f"VCF file {vcf_fp} not found. ")
    
    vcf = cyvcf2.VCF(vcf_fp)
    region = pb.BedTool(region_file)
    variants = extract_variants_within_regions(vcf, region, fully_within=fully_within)
    variants = sorted(variants, key=lambda x: (x.CHROM, x.POS))
    logger.info(f"Finding variants {'fully within' if fully_within else 'partially within'} regions.")
    output_vcf = os.path.join(outdir, "candidates.vcf.gz")
    vcf_out = cyvcf2.Writer(output_vcf, vcf)
    for v in variants:
        vcf_out.write_record(v)
    vcf_out.close()
    vcf.close()
    
    pysam.tabix_index(output_vcf, preset="vcf", force=True)
    return output_vcf

def filter_variants(vcf_fp, out_filtered_vcf, wt_pos, mut_pos, read_length, cipos_tag="CIPOS", ciend_tag="CIEND"):
    """Filter variants based on heterozygous positions."""
    import numpy as np
    
    def is_high_heterozygosity(het: float, mean: float, sd: float, tol: int, p_cutoff=.05) -> bool:
        from scipy.stats import norm
        
        z = (het - mean) / sd
        p = 1 - norm.cdf(z)
        adjusted_p = p * tol

        return adjusted_p < p_cutoff
    
    vcf = cyvcf2.VCF(vcf_fp)
    variant_pos = [v.POS for v in vcf()]
    vcf.close()
    
    vcf = cyvcf2.VCF(vcf_fp)
    het_counts = []
    alt_haplotypes_for_long_deletions = {} ## variant.ID->alt haplotypes
    for variant in vcf():
        if is_het_deletion(variant) and variant.CHROM in contigs:
            if variant.end - variant.POS > 1e6 and any(variant.POS < ovp < variant.end for ovp in variant_pos):
                alt_haplotypes_for_long_deletions[variant.ID] = list(filter(lambda ovp : variant.POS < ovp < variant.end, variant_pos))
            else:
                sv = SV(variant, read_length=read_length, cipos_tag=cipos_tag, ciend_tag=ciend_tag)
                het_counts.append(len(sv.get_hets_within_deletion(wt_pos)))
    vcf.close()
    
    tol = len(het_counts)
    mean, sd = np.mean(het_counts), np.std(het_counts)
    logger.debug(f"mean: {mean:.2f}, sd: {sd:.2f}, tol: {tol}")
    logger.debug(f"het_counts: {het_counts}")
    
    vcf = cyvcf2.VCF(vcf_fp)
    ## Add FILTER headers
    new_filters = {
        "wthap": "More than one haplotypes detected from WT evidence",
        "muthap": "More than one haplotypes detected from MUT evidence",
        "het": "Heterozygous sites within deletion"
    }
    for filter_id, description in new_filters.items():
        vcf.add_filter_to_header({'ID': filter_id, 'Description': description})
        
    out_vcf = cyvcf2.Writer(out_filtered_vcf, vcf)
    
    rejected_count = 0
    for variant in vcf():
        if not is_het_deletion(variant) or variant.CHROM not in contigs:
            variant.FILTER = "PASS" if variant.FILTER is None else variant.FILTER
            out_vcf.write_record(variant)
            continue
        elif variant.end - variant.POS > 1e6 and alt_haplotypes_for_long_deletions.get(variant.ID, None) is not None:
            if is_high_heterozygosity(len(alt_haplotypes_for_long_deletions[variant.ID]), mean, sd, tol):
                logger.info(f"Rejecting {variant.CHROM}:{variant.POS}-{variant.end} due to long deletion with HET DEL within deletion.")
                rejected_count += 1
                variant.FILTER = "het"
                out_vcf.write_record(variant)
                continue
            else:
                variant.FILTER = "PASS" if variant.FILTER is None else variant.FILTER
                out_vcf.write_record(variant)
                continue
        sv = SV(variant, read_length=read_length, cipos_tag=cipos_tag, ciend_tag=ciend_tag)
        
        ## Count heterozygous sites at breakpoints
        wt_hets = sv.associated_hets(wt_pos, distance=5)
        mut_hets = sv.associated_hets(mut_pos, distance=5)
        
        ## Validate heterozygous sites within deletion
        within_hets = sv.get_hets_within_deletion(wt_pos)
        true_het_within_deletion = is_high_heterozygosity(len(within_hets), mean, sd, tol, p_cutoff=0.05)
        
        if 0 < len(wt_hets) + len(mut_hets) < 5:
            logger.info(f"Rejecting {sv} due to {len(wt_hets)} WT hets and {len(mut_hets)} MUT hets.")
            rejected_count += 1
            variant.FILTER = 'wthap' if len(wt_hets) > 0 else 'muthap'
            out_vcf.write_record(variant)
        elif true_het_within_deletion:
            logger.info(f"Rejecting {sv} due to heterozygous sites within deletion. ")
            rejected_count += 1
            variant.FILTER = 'het'
            out_vcf.write_record(variant)
        else:
            variant.FILTER = "PASS" if variant.FILTER is None else variant.FILTER
            out_vcf.write_record(variant)
    logger.info(f"Rejected {rejected_count} variants.")
    out_vcf.close()
    vcf.close()
    pysam.tabix_index(out_filtered_vcf, preset="vcf", force=True)

