# SVhet: An accurate NGS-based filter for moderate-to-large heterozygous deletions

SVhet filters moderate-to-large heterozygous deletions called by NGS-based structural variant callers such as Manta and DELLY. While it is designed for short read sequencing data, the algorithm works with long-read sequencing data in theory. SVhet does not benefit from parallel multi-core computation and is sufficiently performant with minimal memory footprint. SVhet takes a sample-level SV callset (in gz-VCF), the corresponding alignment BAM file, and the short variant calls derived from the alignments as input. Users also need to know the genome build and provide the relevant problematic regions and repeat regions for curating the set of highly confident heterozygous sites. For output, SVhet writes to the specified output directory the set of SVs passing all SVhet filters and the set of rejected SVs.

## Installation
```bash
git clone git@github.com:snakesch/SVhet.git
```

## Dependencies
- bcftools
- bedtools
- pysam
- git-lfs (for downloading large test files)
   
## Testing
A minimal test set was provided with SVhet under `test`. Users may test the functionality of SVhet using the following commands:
```bash
SVhet.py --svcall test/test_sv.vcf.gz --short_var test/test_short_vars.vcf.gz --bam test/test_alignments.bam --output test --problematic-regions ref/problematic_regions_all_hg19.bed.gz --repeat-regions ref/repeatMasker_gt100bp_hg19.bed.gz --threads 4
```
Users can compare the output call set in `test/` with `test/test_out/test_sv.filtered.vcf.gz`.
   
## Usage
1. Input preparation
SVhet is designed for filtering moderate-to-large heterozygous deletions called from short read NGS data. To run SVhet with maximal performance, follow the steps below to prepare your input SV call set.
```bash
## Extract only SVs passed all caller-specific filters and discard unassembled breakends
bcftools filter -i '(FILTER~"PASS") && (SVTYPE!~"BND")' -R autosomes.bed -Wtbi -Oz -o VCF_OUT VCF_IN
```
where `autosomes.bed` is a BED of all autosomal regions.

2. Run SVhet
```bash
usage: SVhet.py [-h] --svcall SVCALL --short-var SHORT_VAR --bam BAM --output OUTPUT --problematic-regions PROBLEMATIC_REGIONS --repeat-regions REPEAT_REGIONS
                [--verbose VERBOSE] [--threads THREADS]

Evaluate structural and short variant calls

options:
  -h, --help            show this help message and exit
  --svcall SVCALL       Candidate structural variant VCF
  --short-var SHORT_VAR
                        Short variant calls
  --bam BAM             Sequence alignments
  --output OUTPUT       Output directory
  --problematic-regions PROBLEMATIC_REGIONS
                        Problematic regions (BED)
  --repeat-regions REPEAT_REGIONS
                        Regions of known repeat elements >100bp (BED)
  --verbose VERBOSE     level of verbosity (default: INFO)
  --threads THREADS     Number of CPU threads (default: 4)
```
## License and contact
SVhet source code is provided under the [license](https://github.com/snakesch/SVhet/blob/main/LICENSE). Correspondence to Louis SHE (snakesch@connect.hku.hk).

