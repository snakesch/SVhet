"""Shared constants and helper utilities for SVhet."""

import gzip, os, random, subprocess
from typing import Optional

import numpy as np
import pysam

# ── Constants ────────────────────────────────────────────────────────────────
MAPQ = 30
LARGE_SVLEN = 1_000_000
MIN_SVLEN = 50
MAX_CI = 150
ISIZE_PROB = 0.001
ISIZE_RUNS = 3
ISIZE_N = 2000
ISIZE_MAXT = 8000
HET_SPACING = 100
FLANK = 100


# ── Helpers ──────────────────────────────────────────────────────────────────
def run_cmd(*cmd) -> None:
    """Run a shell command, raising on failure."""
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)


def is_clipped(r) -> bool:
    """Check if a read is soft/hard clipped on either end."""
    c = r.cigartuples
    return bool(c) and (c[0] in (4, 5) or c[-1] in (4, 5))


def isize(bam_fp) -> tuple[float, float, float]:
    """Estimate insert size distribution by sampling reads."""
    stats = []
    for _ in range(ISIZE_RUNS):
        sizes = []
        with pysam.AlignmentFile(bam_fp, "rb") as bam:
            for r in bam:
                if (r.flag & 0x2 and not r.flag & 0x90C
                        and r.reference_id == r.next_reference_id
                        and abs(r.template_length) <= ISIZE_MAXT
                        and random.random() < ISIZE_PROB):
                    sizes.append(abs(r.template_length))
                    if len(sizes) >= ISIZE_N:
                        break
        if sizes:
            thr = np.percentile(sizes, 99)
            filt = [s for s in sizes if s <= thr][:ISIZE_N]
            if filt:
                stats.append(np.array([np.mean(filt), np.median(filt), np.std(filt)]))
    return tuple(np.mean(np.stack(stats), axis=0)) if stats else (350.0, 0.0, 50.0)


def gt_het(gt) -> bool:
    """Check if a genotype tuple is heterozygous."""
    return gt is not None and gt[0] is not None and gt[1] is not None and gt[0] != gt[1]


def svlen(rec) -> int:
    """Get |SVLEN| from a VCF record."""
    v = rec.info.get("SVLEN")
    if v is not None:
        return abs(v[0] if isinstance(v, tuple) else v)
    return abs(rec.stop - rec.pos)


def ci_sum(ci) -> int:
    """Sum of absolute values of a CI tuple."""
    if ci is None:
        return 0
    return sum(abs(x) for x in ci) if isinstance(ci, (tuple, list)) else abs(ci)


def load_bed(path: str) -> list[tuple[str, int, int]]:
    """Load a BED file into (chrom, start, end) tuples (0-based)."""
    op = gzip.open if path.endswith(".gz") else open
    return [
        (p[0], int(p[1]), int(p[2]))
        for p in (l.strip().split("\t") for l in op(path, "rt") if l.strip() and not l.startswith("#"))
        if len(p) >= 3
    ]


def in_bed(rec, regions) -> bool:
    """Check if a VCF record overlaps any BED region."""
    return any(rec.chrom == c and rec.start < e and rec.stop > s for c, s, e in regions)


def index(path: str) -> None:
    """Tabix-index a VCF/BCF file in place."""
    pysam.tabix_index(path, preset="vcf", force=True)


def has_records(vcf_path: str) -> bool:
    """Check if a VCF/BCF file contains any records."""
    v = pysam.VariantFile(vcf_path)
    has = any(True for _ in v.fetch())
    v.close()
    return has
