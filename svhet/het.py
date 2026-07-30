"""Heterozygosity evaluation and VCF annotation."""

import pysam

from .utils import index, FLANK, HET_SPACING


def reliable_hets(vcf: pysam.VariantFile, chrom: str, start: int, stop: int, min_dp: int) -> list[tuple[str, int]]:
    """Count reliable HET sites (min AD ≥ min_dp, spaced ≥ HET_SPACING) within ±FLANK."""
    hets = []; last = -500
    for v in vcf.fetch(chrom, max(0, start - FLANK), stop + FLANK):
        sid = list(v.samples.keys())[0]
        ad = v.samples[sid].get("AD")
        if ad is None or min(ad) < min_dp:
            continue
        if abs(v.pos - last) > HET_SPACING:
            hets.append((chrom, v.pos))
        elif (chrom, last) in hets:
            hets.remove((chrom, last))
        last = v.pos
    return hets


def evaluate_het(cand_vcf: str, wt_vcf: str, mut_vcf: str, sid: str,
                 outpath: str, min_dp: int = 5, high_het: int = 1) -> str:
    """Annotate candidate VCF with WT_HETS/MUT_HETS/SVHET. Returns annotated BCF path."""
    wt = pysam.VariantFile(wt_vcf); mut = pysam.VariantFile(mut_vcf)
    cand = pysam.VariantFile(cand_vcf); cand.subset_samples([sid])
    hdr = cand.header
    hdr.formats.add("WT_HETS", 1, "Integer", "Number of observed WT heterozygous sites")
    hdr.formats.add("MUT_HETS", 1, "Integer", "Number of observed MUT heterozygous sites")
    hdr.formats.add("SVHET", 1, "String", "SVhet filter")
    out = pysam.VariantFile(outpath, "w", header=hdr)

    for c in cand.fetch():
        gt = c.samples[sid]["GT"]
        wh = mh = []
        if gt == (0, 1):
            wh = reliable_hets(wt, c.chrom, c.start, c.stop, min_dp)
            mh = reliable_hets(mut, c.chrom, c.start, c.stop, min_dp)
        c.samples[sid]["WT_HETS"] = len(wh)
        c.samples[sid]["MUT_HETS"] = len(mh)
        c.samples[sid]["SVHET"] = "HIGH_HET" if len(wh) + len(mh) > high_het else "PASS"
        out.write(c)

    out.close(); wt.close(); mut.close(); cand.close()
    index(outpath)
    return outpath


def passthrough_annotate(cand_vcf: str, sid: str, outpath: str) -> str:
    """Create annotated BCF with all PASS (no WT/MUT callsets available)."""
    v = pysam.VariantFile(cand_vcf); v.subset_samples([sid])
    hdr = v.header
    for tag, n, t, d in [("WT_HETS", 1, "Integer", "WT hets"),
                         ("MUT_HETS", 1, "Integer", "MUT hets"),
                         ("SVHET", 1, "String", "SVhet filter")]:
        hdr.formats.add(tag, n, t, d)
    o = pysam.VariantFile(outpath, "w", header=hdr)
    for r in v.fetch():
        r.samples[sid]["WT_HETS"] = 0
        r.samples[sid]["MUT_HETS"] = 0
        r.samples[sid]["SVHET"] = "PASS"
        o.write(r)
    o.close(); v.close()
    index(outpath)
    return outpath
