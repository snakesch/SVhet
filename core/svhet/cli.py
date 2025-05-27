#!/usr/bin/env python3
"""Command-line interface for SVhet."""

import argparse
import sys

def get_version():
    from importlib.metadata import version
    return version("svhet")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(add_help=False, description='SVhet: An accurate NGS-based structural variant filtering tool using heterozygous sites')
    
    input_args = parser.add_argument_group('Input arguments')
    input_args.add_argument('-i', '--input',
                        required=True,
                        help='Input SV callset')
    input_args.add_argument('-r', '--ref',
                        required=True,
                        help='Reference genome used to call structural variants')
    input_args.add_argument('-b', '--bam',
                        required=True,
                        help='BAM file used for SV calling')
    input_args.add_argument('--image',
                        default=None,
                        help="Singularity image of DeepVariant")
    
    output_args = parser.add_argument_group('Output arguments')
    output_args.add_argument('-d', '--outdir',
                        required=True,
                        help='Output directory')
    output_args.add_argument('-o', '--output',
                        default="filtered.vcf.gz",
                        help='Output file name (default: filtered.vcf.gz)')
   
    caller_args = parser.add_argument_group('SV caller-specific arguments')
    caller_args.add_argument('--cipos-tag',
                        default='CIPOS',
                        help='INFO field tag of CIPOS (default: CIPOS)')
    caller_args.add_argument('--ciend-tag',
                        default='CIEND',
                        help="INFO field tag of CIEND (default: CIEND)")
    
    opts = parser.add_argument_group('Optional arguments')
    opts.add_argument('-t', '--threads',
                        type=int,
                        default=4,
                        help="Number of CPU threads to use")
    opts.add_argument('--region',
                        default=None,
                        help="Selected regions for SVhet filtering (default: None)")
    opts.add_argument('--fully-within',
                        action='store_true',
                        help="Filter variants fully within the target region")
    opts.add_argument('-h', '--help',
                        action='help',
                        help='Show this help message and exit')
    
    debug_args = parser.add_argument_group('Debugging arguments')
    debug_args.add_argument('--force-rerun', 
                        action='store_true',
                        help='Force reprocessing even if output files exist')
    debug_args.add_argument('--log-level',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO',
                        help='Set logging level')
    debug_args.add_argument('--version',
                        action='version',
                        version=f'%(prog)s {get_version()}',
                        help='Show version number and exit')
    
    return parser.parse_args()

def main():
    """Main entry point for the CLI."""
    args = parse_args()
    
    import os
    import logging
    from svhet.utils.log import setup_logger
    
    logging.basicConfig(level=getattr(logging, args.log_level))
    
    logger = setup_logger("svhet", level=args.log_level)
    
    os.makedirs(args.outdir, exist_ok=True)
    
    try:
        import subprocess
        subprocess.check_call(['singularity', '--help'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logger.error(f"Singularity is not installed or not found in PATH. Please install Singularity to use this tool.")
        return 1
    
    try:
        import subprocess
        subprocess.check_call(['bedtools', '--help'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logger.error(f"Bedtools is not installed or not found in PATH. Please install Bedtools to use this tool.")
        return 1

    try:
        from svhet.filtering import filter_sv_callset
        
        filter_sv_callset(
            vcf_fp=args.input,
            bam_fp=args.bam,
            outdir=args.outdir,
            reference_fasta_path=args.ref,
            output_filename=args.output,
            threads=args.threads,
            image_path=args.image,
            region_file=args.region,
            force_rerun=args.force_rerun,
            cipos_tag=args.cipos_tag, 
            ciend_tag=args.ciend_tag,
            fully_within=args.fully_within
        )
        logger.info(f"Done. Results saved to {os.path.join(args.outdir, args.output)}")
        return 0
    except Exception as e:
        logger.error(f"Error during filtering: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
