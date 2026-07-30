"""Merge per-sample and cohort-level annotated VCFs."""

import os

from .utils import index, run_cmd


def merge_sample(tmpdir: str, sid: str, no_cand_bcf: str) -> str:
    """Concat small + large + non-candidates for one sample → annotated BCF in *tmpdir*."""
    small = os.path.join(tmpdir, f"{sid}-ann-small.bcf")
    large = os.path.join(tmpdir, f"{sid}-ann-large.bcf")
    nc = os.path.join(tmpdir, f"{sid}-nc.bcf")
    run_cmd("bcftools", "view", "-s", sid, "-Ob", "-o", nc, no_cand_bcf)
    index(nc)
    merged = os.path.join(tmpdir, f"{sid}-annotated.bcf")
    run_cmd("bcftools", "concat", "-a", "-Ob", "-o", merged, small, large, nc)
    index(merged)
    return merged


def merge_cohort(annotated_bcfs: list[str], final_path: str) -> str:
    """Merge all per-sample annotated BCFs → final cohort VCF."""
    run_cmd("bcftools", "merge", "-Oz", "-o", final_path, "--force-single", *annotated_bcfs)
    index(final_path)
    return final_path
