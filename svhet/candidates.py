"""Generate deletion candidates from a cohort VCF."""

import logging, os

import pysam

from .utils import ci_sum, gt_het, in_bed, index, load_bed, svlen, LARGE_SVLEN, MAX_CI, MIN_SVLEN

log = logging.getLogger("svhet")


def generate_candidates(sv_vcf: str, outdir: str, bed: str | None = None) -> dict[str, str]:
    """Split cohort VCF → small/large DEL candidates + non-candidates.

    Returns dict with keys 'small_candidates', 'large_candidates', 'no_candidates' → BCF paths.
    """
    bed_r = load_bed(bed) if bed else None
    vcf = pysam.VariantFile(sv_vcf)
    samples = list(vcf.header.samples)
    keys = ("small_candidates", "large_candidates", "no_candidates")
    paths = {k: os.path.join(outdir, f"{k}.bcf") for k in keys}
    outs = {k: pysam.VariantFile(p, "wb", header=vcf.header) for k, p in paths.items()}
    counts = dict.fromkeys(keys, 0)

    for rec in vcf.fetch():
        if bed_r and not in_bed(rec, bed_r):
            outs["no_candidates"].write(rec); counts["no_candidates"] += 1; continue
        is_del = rec.info.get("SVTYPE") == "DEL"
        has_het = any(gt_het(rec.samples[s]["GT"]) for s in samples if s in rec.samples)
        if not (is_del and has_het):
            outs["no_candidates"].write(rec); counts["no_candidates"] += 1; continue
        svl = svlen(rec)
        cipos = rec.info.get("CIPOS"); ciend = rec.info.get("CIEND", cipos)
        if svl < MIN_SVLEN or ci_sum(cipos) > MAX_CI or ci_sum(ciend) > MAX_CI:
            outs["no_candidates"].write(rec); counts["no_candidates"] += 1; continue
        k = "large_candidates" if svl >= LARGE_SVLEN else "small_candidates"
        outs[k].write(rec); counts[k] += 1

    for o in outs.values():
        o.close()
    vcf.close()
    for p in paths.values():
        index(p)
    log.info("Candidates: small=%d, large=%d, non=%d",
             counts["small_candidates"], counts["large_candidates"], counts["no_candidates"])
    return paths
