"""SVhet pipeline orchestration — ties all steps together."""

import logging, os, shutil, tempfile

from .candidates import generate_candidates
from .reads import extract_reads
from .call import call_variants
from .het import evaluate_het, passthrough_annotate
from .merge import merge_cohort, merge_sample
from .utils import has_records

log = logging.getLogger("svhet")


def run_svhet(ref: str, sv_vcf: str, outdir: str, manifest: str,
              bed: str | None = None, jobs: int = 1,
              min_dp: int = 5, high_het: int = 1) -> str:
    """Run the full SVhet pipeline.

    Only ``final-annotated.vcf.gz`` is written to *outdir*;
    all intermediates go to a temp directory that is cleaned up automatically.
    Returns the path to the final annotated VCF.
    """
    os.makedirs(outdir, exist_ok=True)
    tmpdir = tempfile.mkdtemp(dir=outdir, prefix=".svhet_tmp_")

    try:
        # Step 1 — candidate generation
        log.info("Step 1: Generating candidates")
        cand = generate_candidates(sv_vcf, tmpdir, bed=bed)

        # Read manifest
        samples = []
        with open(manifest) as f:
            for line in f:
                p = line.strip().split("\t")
                if p and p[0]:
                    samples.append((p[0], p[1]))
        log.info("Manifest: %d samples", len(samples))

        annotated_bcfs = []

        for sid, bam_fp in samples:
            log.info("─" * 50)
            log.info("Sample: %s", sid)

            for mode in ("small", "large"):
                cv = cand[f"{mode}_candidates"]
                ann = os.path.join(tmpdir, f"{sid}-ann-{mode}.bcf")
                if not has_records(cv):
                    log.info("  No %s candidates — skipping", mode)
                    passthrough_annotate(cv, sid, ann)
                    continue

                # Step 2 — extract reads
                log.info("  Step 2: Extracting reads (%s)", mode)
                wt_bam, mut_bam = extract_reads(cv, sid, bam_fp, tmpdir, jobs, mode)
                if wt_bam is None:
                    passthrough_annotate(cv, sid, ann)
                    continue

                # Step 3 — call variants
                log.info("  Step 3: Calling variants (%s)", mode)
                wt_bcf = call_variants(ref, wt_bam, jobs)
                mut_bcf = call_variants(ref, mut_bam, jobs)

                # Step 4 — evaluate heterozygosity
                log.info("  Step 4: Evaluating heterozygosity (%s)", mode)
                evaluate_het(cv, wt_bcf, mut_bcf, sid, ann, min_dp, high_het)

            # Step 5a — merge per-sample
            log.info("  Step 5a: Merging per-sample")
            annotated_bcfs.append(merge_sample(tmpdir, sid, cand["no_candidates"]))

        # Step 5b — merge cohort
        log.info("Step 5b: Merging cohort → final-annotated.vcf.gz")
        final = merge_cohort(annotated_bcfs, os.path.join(outdir, "final-annotated.vcf.gz"))
        log.info("Done! Output: %s", final)
        return final

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
