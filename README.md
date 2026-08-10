<p align="center">
   <img src="logo/SVhet-logo.png" alt="SVhet Logo" width="250"/>
</p>



# SVhet: Structural Variant Filtering using Heterozygosity

SVhet is a pipeline for filtering **heterozygous deletion** calls in cohort-level VCFs using evidence from heterozygous sites in short-read sequencing data. It improves the reliability of heterozygous deletion calls by leveraging read-level and variant-level information across samples. SVhet does not filter other types of structural variants.

---

## Table of Contents

1. [Features](#features)
2. [Installation & Dependencies](#installation--dependencies)
3. [Usage](#usage)
4. [Pipeline Overview](#pipeline-overview)
5. [Implementation Details](#implementation-details)
6. [Output Format](#output-format)
7. [Example Files & Run](#example-files--run)
8. [Citation](#citation)
9. [Contact](#contact)

---

## Features

- Filters **heterozygous deletions** based on genotype and quality metrics
- Per-sample read evidence extraction for wild-type (WT) and mutant (MUT) alleles
- Short variant calling on extracted read sets
- Heterozygosity evaluation within deletion regions and their flanking regions to flag unreliable calls
- Produces a single, annotated cohort-level VCF for downstream analysis
- Intermediates are written to a temporary directory inside `--outdir` and removed on exit

---

## Installation & Dependencies

**Runtime requirements**

- Python 3.10+ (the code uses PEP 604 `X | None` annotations)
- [pysam](https://pysam.readthedocs.io/)
- numpy
- [bcftools](http://samtools.github.io/bcftools/) on `PATH` (used for `view`, `concat`, `merge`, `mpileup`, `call`, `filter`)

`samtools sort`/`index` are invoked through pysam's bundled htslib, so a separate samtools installation is not required.

**From source**

```bash
git clone https://github.com/snakesch/SVhet.git
cd SVhet
pip install pysam numpy
python svhet.py --help
```

**Docker**

The provided [Dockerfile](Dockerfile) builds htslib, bcftools, samtools and bedtools from source and installs `numpy==1.26.4` and `pysam==0.22.1`.

```bash
docker build -t svhet .
docker run --rm -it svhet python svhet.py --help
```

---

## Usage

```bash
python svhet.py --ref <reference.fasta> \
                --sv-vcf <cohort.vcf.gz> \
                --outdir <output_dir> \
                --manifest <manifest.txt> \
                [--bed <regions.bed>] [--jobs <N>] [--min-dp <N>] [--high-het <N>] [--verbose]
```

**Required arguments**

- `--ref` : Reference FASTA file
- `--sv-vcf` : Cohort-level SV VCF file (bgzipped and indexed)
- `--outdir` : Output directory
- `--manifest` : Tab-delimited file with sample ID and BAM path per line

**Optional arguments**

- `--bed` : BED file of target regions; records outside these regions bypass filtering
- `--jobs` : Number of parallel workers and bcftools threads (default: 1)
- `--min-dp` : Minimum allelic depth for a HET site to be considered reliable (default: 5)
- `--high-het` : `WT_HETS + MUT_HETS` above this value flags a call as `HIGH_HET` (default: 1)
- `--verbose` : DEBUG-level logging

The pipeline can also be driven from Python:

```python
from svhet import run_svhet

run_svhet(ref, sv_vcf, outdir, manifest, bed=None, jobs=1, min_dp=5, high_het=1)
```

---

## Pipeline Overview

`svhet.py` parses arguments and delegates to `run_svhet` in [svhet/pipeline.py](svhet/pipeline.py), which runs five steps.

1. **Candidate generation** ([svhet/candidates.py](svhet/candidates.py))
   Splits the cohort VCF into `small_candidates.bcf`, `large_candidates.bcf` and `no_candidates.bcf`. A record becomes a candidate only if `SVTYPE=DEL`, at least one sample is heterozygous, `|SVLEN| >= 50`, and both the `CIPOS` and `CIEND` confidence intervals sum to `<= 150`. Candidates with `|SVLEN| >= 1e6` go to the large set. Everything else, including records outside `--bed`, goes to the non-candidate set and passes through unfiltered.

2. **Read evidence extraction** ([svhet/reads.py](svhet/reads.py))
   For each sample the insert size distribution is estimated by random sampling, then WT and MUT read sets are collected per candidate and written to sorted, indexed BAMs. Candidates are distributed across a `ProcessPoolExecutor` when `--jobs > 1`; workers return reads as SAM strings that the parent re-materialises against the source BAM header.

3. **Short variant calling** ([svhet/call.py](svhet/call.py))
   `bcftools mpileup | bcftools call -mv --ploidy 2 | bcftools filter -i 'GT="het"'` is streamed through pipes for each WT and MUT BAM, producing indexed BCFs of heterozygous sites.

4. **Heterozygosity evaluation** ([svhet/het.py](svhet/het.py))
   For every heterozygous candidate, reliable HET sites are counted in the WT and MUT callsets over the deletion interval extended by `FLANK = 100` bp. A site counts only if `min(AD) >= --min-dp`; sites closer than `HET_SPACING = 100` bp to the previous site are collapsed to suppress clustered artefacts. Counts are written as `WT_HETS`/`MUT_HETS`, and `SVHET` is set to `HIGH_HET` when their sum exceeds `--high-het`. Modes with no candidates or no het carriers are annotated all-`PASS` via `passthrough_annotate`.

5. **Merging** ([svhet/merge.py](svhet/merge.py))
   Per sample, the small, large and non-candidate BCFs are concatenated with `bcftools concat -a`. All per-sample BCFs are then merged with `bcftools merge --force-single` into `final-annotated.vcf.gz`, which is tabix-indexed.

---

## Implementation Details

### Package layout

```
svhet.py              CLI entry point
svhet/
  __init__.py         exports run_svhet, __version__
  pipeline.py         orchestration, manifest parsing, temp dir lifecycle
  candidates.py       cohort VCF -> small/large/non-candidate BCFs
  reads.py            WT/MUT read extraction, parallel workers
  call.py             bcftools mpileup/call/filter wrapper
  het.py              reliable HET counting and VCF annotation
  merge.py            per-sample concat and cohort merge
  utils.py            constants and shared helpers
```

### WT read criteria

A read supports the reference (WT) allele when it passes `MAPQ >= 30`, is not a duplicate, is not soft/hard clipped, carries no `SA` or `XA` tag, has `NM <= 5`, and has a non-zero template length. Reads spanning either breakpoint (within `CIPOS`/`CIEND`) qualify when `| |TLEN| - mean_isize | < sd_isize`. For small candidates, reads fully contained inside the deletion interval are also collected, since a genuine heterozygous deletion still retains coverage from the intact haplotype.

### MUT read criteria

A read supports the deletion when it passes `MAPQ >= 30`, is not a duplicate, and either has `|TLEN| >= mean_isize + sd_isize` (discordant pair) or carries an `SA` tag whose supplementary alignment is on the same contig and whose offset from the read end matches the deletion length within 500 bp (split read). Small candidates are scanned across the whole interval; large candidates are scanned only around the two breakpoints to bound I/O.

### Insert size estimation

`utils.isize` performs 3 sampling passes over the BAM, accepting properly paired primary reads on the same contig with `|TLEN| <= 8000` at a sampling probability of 0.001, up to 2000 observations per pass. Each pass trims values above the 99th percentile before computing mean, median and standard deviation; the three passes are averaged. If no reads qualify, it falls back to `(350, 0, 50)`.

### Constants

Defined in [svhet/utils.py](svhet/utils.py):

| Constant | Value | Meaning |
| --- | --- | --- |
| `MAPQ` | 30 | Minimum mapping quality |
| `MIN_SVLEN` | 50 | Minimum deletion length |
| `LARGE_SVLEN` | 1,000,000 | Small/large candidate split |
| `MAX_CI` | 150 | Maximum summed breakpoint confidence interval |
| `HET_SPACING` | 100 | Minimum spacing between counted HET sites (bp) |
| `FLANK` | 100 | Interval padding for HET counting (bp) |
| `ISIZE_PROB` / `ISIZE_RUNS` / `ISIZE_N` / `ISIZE_MAXT` | 0.001 / 3 / 2000 / 8000 | Insert size sampling parameters |

### Concurrency

`--jobs` controls both the `ProcessPoolExecutor` used for per-candidate read extraction and the `--threads` value passed to `bcftools mpileup`/`call`. WT and MUT variant calling run sequentially within each mode.

---

## Output Format

A single bgzipped, tabix-indexed cohort-level VCF (`final-annotated.vcf.gz`) is written to `--outdir`. All intermediates live in a `.svhet_tmp_*` directory inside `--outdir` and are deleted when the run finishes.

### SVhet-specific FORMAT annotations

- `SVHET`
   - `PASS`: Variant passes SVhet filtering
   - `HIGH_HET`: Variant flagged due to high heterozygosity in the region
- `WT_HETS`: Number of reliable heterozygous sites in the wild-type (WT) allele region (per sample)
- `MUT_HETS`: Number of reliable heterozygous sites in the mutant (MUT) allele region (per sample)

### Minimal example output

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	123456	sv1	N	<DEL>	60	PASS	.	GT:WT_HETS:MUT_HETS:SVHET	0/1:3:0:PASS	0/0:0:0:PASS
chr2	234567	sv2	N	<DEL>	50	PASS	.	GT:WT_HETS:MUT_HETS:SVHET	0/1:5:2:HIGH_HET	0/1:4:0:PASS
```

In the example above, `sv2` should be excluded from downstream analysis due to high heterozygosity detected from WT and MUT read evidence. For true heterozygous deletions, `WT_HETS` and `MUT_HETS` are typically 0 since only one haplotype exists in the deleted region.

---

## Example Files & Run

**Manifest format**

Tab-delimited, one sample per line. Only the first two fields are read; a third BAI column is accepted and ignored.

```
HG00096	/path/to/HG00096.bam	/path/to/HG00096.bam.bai
HG00097	/path/to/HG00097.bam	/path/to/HG00097.bam.bai
```

**Example run**

Download the T2T reference from [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz), decompress it, then run:

```bash
python svhet.py --ref chm13.v2.fasta \
                --sv-vcf test/chr1_127510500_128695280_HG00096.vcf.gz \
                --outdir test/results \
                --manifest test/manifest.txt \
                --jobs 4
```

`test/results/final-annotated.vcf.gz` should match `test/test_output/final-annotated.vcf.gz`. Use absolute paths in the manifest if no output is produced.

---

## Citation

If you use SVhet in your research, please cite:

> She, C.H., Chan, S.HS. & Yang, W. SVhet: towards accurate detection of germline heterozygous deletions using short reads. BMC Bioinformatics (2025). https://doi.org/10.1186/s12859-025-06342-7

---

## Contact

For questions or issues, please contact Louis (snakesch@connect.hku.hk).
