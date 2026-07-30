"""Short variant calling on WT/MUT read BAMs."""

import subprocess

from .utils import index


def call_variants(ref: str, bam_path: str, jobs: int = 1) -> str:
    """Call short variants on a BAM using bcftools mpileup + call.

    Filters for heterozygous sites only. Returns the output BCF path.
    """
    bcf = bam_path.replace(".bam", ".bcf")
    mp = subprocess.Popen(
        ["bcftools", "mpileup", f"--threads={jobs}", "-f", ref, bam_path],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    ca = subprocess.Popen(
        ["bcftools", "call", f"--threads={jobs}", "-mv", "--ploidy", "2"],
        stdin=mp.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    mp.stdout.close()
    fi = subprocess.Popen(
        ["bcftools", "filter", "-i", 'GT="het"', "-Ob", "-o", bcf],
        stdin=ca.stdout, stderr=subprocess.DEVNULL,
    )
    ca.stdout.close(); fi.wait()
    index(bcf)
    return bcf
