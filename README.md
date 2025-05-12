# SVhet

An accurate NGS-based structural variation filtering tool based on heterozygous sites on reads mapped to deleted regions and their flanks. SVhet has been tested on common SV callers including Manta, DELLY and Lumpy (both PE150 and PE250). While SVhet works for all regions in theory, it is only tested on GIAB v0.6 Tier 1 SV regions, with either/both breakpoint(s) located within the regions. To use SVhet with minimal loss in recall and best performance, users are recommended to include only SVs within Tier 1 regions and not specify `--fully-within`.

GIAB v0.6 Tier 1 SV benchmark regions can be downloaded [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed).

# Installation

SVhet can be downloaded and installed using conda or mamba. This will automatically configure the environment and install most dependencies of SVhet.

```bash
git clone git@github.com:snakesch/SVhet.git
cd SVhet
conda env create -f environment.yml
conda activate svhet
```

Alternatively, SVhet can be downloaded from PyPI:
```bash
pip install svhet
```

Next, make sure singularity (v3.0+) is available on the system and pull the DeepVariant image:

```bash
singularity pull docker://google/deepvariant:1.9.0
```

# Running SVhet

SVhet can be run as a single command from CLI and takes as input candidate SVs proposed by state-of-the-art SV callers (both VCF and BAM) and a reference genome. Depending on the specific SV caller used, users may adopt different CIPOS and CIEND tags (e.g. CIPOS95 for Lumpy).

```bash
svhet \
--input VCF \
--bam BAM \
--ref REFERENCE \
--output VCF_OUT \
--image DEEPVARIANT_IMAGE \
--outdir OUTPUT_DIR \
--threads INT \
--cipos-tag CIPOS \
--ciend-tag CIEND
```

A list of all available arguments can be accessed from `svhet --help`.

# Test run

A minimal test set is available from `test/`. To run the test case, download hs37d5 reference and its index file from [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/)

Test command:
```bash
svhet \
--input test.vcf.gz \
--bam test.bam \
--ref REFERENCE \
--output test.svhet.vcf.gz \
--image DEEPVARIANT_IMAGE \
--outdir test_output/ \
--threads 4 \
--cipos-tag CIPOS \
--ciend-tag CIEND
```

Upon successful execution, output files should be similar to those in `test/test_output/`.

# License

SVhet is available under an [MIT license](LICENSE).

# Issues and correspondence

Issues and correspondence to Louis SHE (snakesch@connect.hku.hk).

# Citation

Preprint pending.
