#!/usr/bin/env python3
"""SVhet — Structural Variant Filtering using Heterozygosity (Python edition v0.2.0).

Usage:
  python svhet.py --ref <ref.fasta> --sv-vcf <cohort.vcf.gz> --outdir <out/> --manifest <manifest.txt> [OPTIONS]
"""

import argparse, logging, sys

from svhet import run_svhet, __version__


def main():
    p = argparse.ArgumentParser(description="SVhet — SV Filtering using Heterozygosity (Python)")
    p.add_argument("--ref", required=True)
    p.add_argument("--sv-vcf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--bed", default=None)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--min-dp", type=int, default=5)
    p.add_argument("--high-het", type=int, default=1)
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    print(f"SVhet v{__version__}\n")
    run_svhet(args.ref, args.sv_vcf, args.outdir, args.manifest,
              args.bed, args.jobs, args.min_dp, args.high_het)


if __name__ == "__main__":
    main()
