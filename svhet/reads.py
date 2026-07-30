"""Extract per-sample WT and MUT read evidence from BAMs."""

import itertools, logging, os
from concurrent.futures import ProcessPoolExecutor

import pysam

from .utils import is_clipped, isize, MAPQ

log = logging.getLogger("svhet")


def _variant_pool(cand_vcf: str, sid: str) -> dict:
    """Read candidate VCF and build variant pool for a sample."""
    vcf = pysam.VariantFile(cand_vcf)
    vcf.subset_samples([sid])
    pool = {}
    for rec in vcf.fetch():
        if "CIPOS" not in rec.info:
            continue
        ci = rec.info["CIPOS"]
        ce = rec.info.get("CIEND", ci)
        key = (rec.chrom, (rec.start, rec.stop), rec.ref,
               rec.alts[0] if rec.alts else ".", tuple(ci), tuple(ce) if ce else tuple(ci))
        gt = rec.samples[sid]["GT"]
        pool[key] = ("not_analysed" if gt == (0, 1) else "PASS")
    vcf.close()
    return pool


def _wt_reads(cand, bam, mu, sd, mode):
    """Extract WT (normal insert size, no clip, no SA/XA, NM≤5) reads."""
    chrom, (s, e), _, _, sci, eci = cand
    rid = set(); reads = []
    for reg in [(chrom, s - abs(sci[0]), s + sci[1]), (chrom, e - abs(eci[0]), e + eci[1])]:
        for r in bam.fetch(*reg):
            if r.mapping_quality < MAPQ or r.is_duplicate or is_clipped(r) or r.template_length == 0:
                continue
            if r.has_tag("XA") or r.has_tag("SA"):
                continue
            if r.get_tag("NM") > 5:
                continue
            if abs(abs(r.template_length) - mu) < sd:
                k = (r.query_name, r.flag)
                if k not in rid:
                    reads.append(r); rid.add(k)
    if mode == "small":
        iS, iE = s + abs(sci[1]), e - abs(eci[0])
        if iS < iE:
            for r in bam.fetch(chrom, iS, iE):
                if r.mapping_quality < MAPQ or r.is_duplicate or is_clipped(r):
                    continue
                if r.has_tag("XA") or r.has_tag("SA"):
                    continue
                if r.get_tag("NM") > 5:
                    continue
                if r.reference_start > iS and r.reference_end < iE:
                    k = (r.query_name, r.flag)
                    if k not in rid:
                        reads.append(r); rid.add(k)
    return reads


def _mut_reads(cand, bam, mu, sd, mode):
    """Extract MUT (abnormal insert size or split-read) reads."""
    chrom, (s, e), _, _, sci, eci = cand
    rid = set(); reads = []
    regs = [(chrom, s - abs(sci[0]), e + eci[1])] if mode == "small" else \
           [(chrom, s - abs(sci[0]), s + sci[1]), (chrom, e - abs(eci[0]), e + eci[1])]
    for r in itertools.chain(*(bam.fetch(*rg) for rg in regs)):
        if r.mapping_quality < MAPQ or r.is_duplicate:
            continue
        k = (r.query_name, r.flag)
        if k in rid:
            continue
        if abs(r.template_length) >= mu + sd:
            reads.append(r); rid.add(k)
        elif r.has_tag("SA"):
            sa = r.get_tag("SA").split(",")
            if sa[0] == chrom and abs(int(sa[1]) - r.reference_end) - abs(e - s) < 500:
                reads.append(r); rid.add(k)
    return reads


def _worker(args):
    """Worker for parallel read extraction — returns serialised read strings."""
    cand, bam_fp, mu, sd, mode = args
    bam = pysam.AlignmentFile(bam_fp, "rb")
    wt = _wt_reads(cand, bam, mu, sd, mode)
    mut = _mut_reads(cand, bam, mu, sd, mode)
    bam.close()
    return [r.to_string() for r in wt], [r.to_string() for r in mut]


def extract_reads(cand_vcf: str, sid: str, bam_fp: str, tmpdir: str,
                  jobs: int = 4, mode: str = "small") -> tuple[str | None, str | None]:
    """Extract WT/MUT reads → sorted BAMs in *tmpdir*.

    Returns ``(wt_bam, mut_bam)`` or ``(None, None)`` if no het carriers.
    """
    pool = _variant_pool(cand_vcf, sid)
    mu, _, sd = isize(bam_fp)
    log.info("  [%s] insert size: %.0f±%.0f", mode, mu, sd)
    tasks = [k for k, st in pool.items() if st == "not_analysed"]
    if not tasks:
        return None, None

    bam_ref = pysam.AlignmentFile(bam_fp)
    wt_p = os.path.join(tmpdir, f"{sid}-wt-{mode}.bam")
    mut_p = os.path.join(tmpdir, f"{sid}-mut-{mode}.bam")
    wt = pysam.AlignmentFile(wt_p + ".tmp", "wb", header=bam_ref.header)
    mut = pysam.AlignmentFile(mut_p + ".tmp", "wb", header=bam_ref.header)

    if jobs > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for wt_s, mut_s in ex.map(_worker, [(t, bam_fp, mu, sd, mode) for t in tasks]):
                for s in wt_s:
                    wt.write(pysam.AlignedSegment.fromstring(s, bam_ref.header))
                for s in mut_s:
                    mut.write(pysam.AlignedSegment.fromstring(s, bam_ref.header))
    else:
        bam_ref.close(); bam_ref = pysam.AlignmentFile(bam_fp, "rb")
        for t in tasks:
            for r in _wt_reads(t, bam_ref, mu, sd, mode):
                wt.write(r)
            for r in _mut_reads(t, bam_ref, mu, sd, mode):
                mut.write(r)

    wt.close(); mut.close(); bam_ref.close()
    for p in (wt_p, mut_p):
        pysam.samtools.sort("-o", p, p + ".tmp", catch_stdout=False)
        pysam.index(p); os.remove(p + ".tmp")
    return wt_p, mut_p
