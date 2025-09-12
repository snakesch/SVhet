
# SVhet: Structural Variant Filtering using Heterozygosity

SVhet is a pipeline for filtering **heterozygous deletion** calls in cohort-level VCFs using evidence from heterozygous sites in short-read sequencing data. It is designed to improve the reliability of heterozygous deletion calls by leveraging read-level and variant-level information across samples. SVhet does not filter other types of structural variants.

## Features
- Filters **only heterozygous deletions** based on genotype and quality metrics
- Per-sample read evidence extraction for wild-type (WT) and mutant (MUT) alleles
- Short variant calling on extracted read sets
- Heterozygosity evaluation within deletion regions and their flanking regions to flag unreliable calls
- Merging and annotation of final VCFs for downstream analysis

## Pipeline Overview
The main entry point is `svhet.sh`, which orchestrates the following steps:

1. **Generate SV Candidates** (`01_generate_candidates.sh`)
   - Filters cohort VCF for deletion candidates with at least one heterozygous carrier.
   - Splits candidates by SV length (small/large) and applies additional quality filters.
   - Optionally restricts to target regions using a BED file.

2. **Extract Per-Sample Read Evidence** (`02_filter_by_samples.py`)
   - For each sample and candidate, extracts WT and MUT supporting reads from the BAM file.
   - Writes these reads to temporary BAMs for downstream variant calling.
   - Handles both small and large SV candidates.

3. **Short Variant Calling** (`03_call_variants.sh`)
   - Calls short variants (SNPs/indels) on the WT and MUT BAMs using `bcftools mpileup` and `bcftools call`.
   - Filters for heterozygous sites.

4. **Heterozygosity Evaluation** (`04_het_evaluator.py`)
   - Compares the number of reliable heterozygous sites in WT and MUT callsets within each SV region.
   - Annotates the candidate VCF with the number of HETs and a filter status (PASS or HIGH_HET).

5. **Merging and Final Output**
   - Annotated VCFs for each sample are merged with non-candidate SVs.
   - All sample-level VCFs are merged into a final, cohort-level annotated VCF.

## Usage

```bash
bash svhet.sh --ref <reference.fasta> \
              --sv-vcf <cohort.vcf.gz> \
              --outdir <output_dir> \
              --manifest <manifest.txt> \
              [--bed <regions.bed>] [--jobs <N>] [--keep-intermediate] [--min-dp <N>] [--high-het <N>]
```

### Required Arguments
- `--ref` : Reference FASTA file
- `--sv-vcf` : Cohort-level SV VCF file (bgzipped and indexed)
- `--outdir` : Output directory
- `--manifest` : Tab-delimited file with sample ID (required), BAM path (required), and BAI path per line (optional)

### Optional Arguments
- `--bed` : BED file of target regions
- `--jobs` : Number of parallel jobs (default: 1)
- `--keep-intermediate` : Keep intermediate files
- `--min-dp` : Minimum depth for reliable HETs (default: 5)
- `--high-het` : Minimum HET count to reject a DEL (default: 1)

## Output
- Annotated per-sample VCFs with SVHET status and HET counts
- Final merged, annotated VCF (`final-annotated.vcf.gz`)

## Dependencies
- [bcftools](http://samtools.github.io/bcftools/)
- [bedtools](https://bedtools.readthedocs.io/)
- [pysam](https://pysam.readthedocs.io/)
- Python 3.6+
- numpy, tqdm

## File Descriptions
- `svhet.sh` : Main pipeline script (entry point)
- `01_generate_candidates.sh` : Candidate SV filtering and splitting
- `02_filter_by_samples.py` : Per-sample read evidence extraction
- `03_call_variants.sh` : Short variant calling on read sets
- `04_het_evaluator.py` : Heterozygosity evaluation and annotation

## Example Manifest File
```
HG00096\t/path/to/HG00096.bam\t/path/to/HG00096.bam.bai
HG00097\t/path/to/HG00097.bam\t/path/to/HG00097.bam.bai
```

## Example Run
```bash
bash svhet.sh --ref chm13.v2.fasta \
             --sv-vcf lumpy.vcf.gz \
             --outdir results/ \
             --manifest manifest.txt \
             --jobs 4
```

## Citation
If you use SVhet in your research, please cite:

> CH She, et al. SVhet: Heterozygosity-based filtering of structural variant calls in cohorts. (2025)

---

For questions or issues, please contact Louis (snakesch@connect.hku.hk).
