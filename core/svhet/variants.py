import warnings
import os

import cyvcf2
import pybedtools as pb
from svhet.utils.log import setup_logger

logger = setup_logger("svhet")

def extract_het_positions(wt_vcf, mut_vcf, ad2dp=0.4, min_dp=5):
    wtvcf = cyvcf2.VCF(wt_vcf)
    mutvcf = cyvcf2.VCF(mut_vcf)

    wt_pos, mut_pos = [], []
    for v in wtvcf():
        if v.FILTER is not None or v.num_het == 0:
            continue
        if min(v.format("AD")[0][0], v.format("AD")[0][1]) / v.format("DP")[0] < ad2dp:
            continue
        if v.format("DP")[0] <= min_dp:
            continue
        wt_pos.append((v.CHROM, v.POS, min(v.format("AD")[0][0], v.format("AD")[0][1]), v.format("DP")[0][0]))
    for v in mutvcf():
        if v.FILTER is not None or v.num_het == 0:
            continue
        if min(v.format("AD")[0][0], v.format("AD")[0][1]) / v.format("DP")[0] < ad2dp:
            continue
        if v.format("DP")[0] <= min_dp:
            continue
        mut_pos.append((v.CHROM, v.POS, min(v.format("AD")[0][0], v.format("AD")[0][1]), v.format("DP")[0][0]))

    wtvcf.close()
    mutvcf.close()
    
    return wt_pos, mut_pos

def extract_variants_within_regions(vcf: cyvcf2.VCF, bed: pb.BedTool, fully_within=True):
    """Extracts DEL HET variants from VCF fully within BED regions."""
    
    print("Extracting variants with fully_within set to ", fully_within)
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
 
        filtered_variants_tmp = {}
        for region in bed:
            chrom = region.chrom
            start_1based = region.start + 1  # BED is 0-based, VCF is 1-based
            end_1based = region.stop + 1

            for variant in vcf(f"{chrom}:{start_1based}-{end_1based}"):
                # if variant.INFO.get("SVTYPE") != "DEL":
                #     continue
                # elif variant.gt_types[0] != 1:
                #     continue
                
                if fully_within:
                    if variant.POS >= start_1based and variant.end <= end_1based:
                        filtered_variants_tmp[variant.ID] = variant
                else:
                    if start_1based < variant.POS < end_1based or start_1based < variant.end < end_1based:
                        filtered_variants_tmp[variant.ID] = variant                   

        return list(filtered_variants_tmp.values())
    
def is_het_deletion(variant: cyvcf2.Variant):
    ## Checks if the variant is a het deletion with PASS filters and size >= 50bp
    size = abs(variant.INFO.get("SVLEN", None) or variant.INFO.get("END", variant.POS) - variant.POS) or 0
    return (variant.INFO.get("SVTYPE") == "DEL" and variant.gt_types[0] == 1 and 
            variant.FILTER is None and size >= 50)

def call_short_variants(wt_bam_fp, mut_bam_fp, candidate_regions, 
                        image_path, reference_fasta_path, threads):
    import subprocess
    
    cmd = [
        "singularity",
        "run",
        "-B",
        f"{os.path.dirname(reference_fasta_path)}:/mnt/ref",
        "-B",
        f"{os.path.dirname(wt_bam_fp)}:/mnt/work",
        image_path,
        "/opt/deepvariant/bin/run_deepvariant",
        f"--model_type=WGS",
        f"--regions=/mnt/work/{os.path.basename(candidate_regions)}",
        f"--ref=/mnt/ref/{os.path.basename(reference_fasta_path)}",
        f"--reads=/mnt/work/{os.path.basename(wt_bam_fp)}",
        f"--output_vcf=/mnt/work/wt_evidence.vcf.gz",
        f"--num_shards={min(threads, 4)}",
        # "--disable_small_model"
    ]

    logger.info("Start variant calling on WT evidence ... ")
    result = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL,
        check=False  
    )
    if result.returncode != 0:
        raise RuntimeError("Error running deepvariant on WT evidence.")

    cmd = [
        "singularity",
        "run",
        "-B",
        f"{os.path.dirname(reference_fasta_path)}:/mnt/ref",
        "-B",
        f"{os.path.dirname(mut_bam_fp)}:/mnt/work",
        image_path,
        "/opt/deepvariant/bin/run_deepvariant",
        f"--model_type=WGS",
        f"--ref=/mnt/ref/{os.path.basename(reference_fasta_path)}",
        f"--reads=/mnt/work/{os.path.basename(mut_bam_fp)}",
        f"--regions=/mnt/work/{os.path.basename(candidate_regions)}",
        f"--output_vcf=/mnt/work/mut_evidence.vcf.gz",
        f"--num_shards={min(threads, 4)}",
        # "--disable_small_model"
    ]


    logger.info("Start variant calling on MUT evidence ... ")
    result = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL,
        check=False  
    )
    if result.returncode != 0:
        raise RuntimeError("Error running deepvariant on MUT evidence.")
    
    return 