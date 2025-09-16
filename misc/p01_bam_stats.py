#!/usr/bin/env python3

## Baseline statistics 
## For each BAM, we need
## coverage, insert size mean, std, read length

import subprocess
import pandas as pd

manifest = "/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams/manifest.txt"
sample_manifest = pd.read_csv(manifest, sep="\t", header=None, usecols=[0,1], names=["sample", "bam fp"])

sample_manifest["Q30 depth"] = 0.
sample_manifest["isize mean"] = 0.
sample_manifest["isize std"] = 0.
sample_manifest["read length"] = 150

def get_insert_dist(bam_file, num_runs=1, num_samples=2000, max_template_length=8000, sample_prob=0.001):
    """
    Runs insert size estimation multiple times, averages statistics, and handles filtering.

    Returns:
        tuple: (Average Mean, Average Median, Average Std), or (None, None, None) if no data.
    """
    import random
    import pysam
    import numpy as np

    all_stats = []
    for _ in range(num_runs):
        insert_sizes = []

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if (read.flag & 0x2 and not (read.flag & 0x90C) and
                    read.reference_id == read.next_reference_id and
                    abs(read.template_length) <= max_template_length and
                    random.random() < sample_prob):

                    insert_sizes.append(abs(read.template_length))
                    if len(insert_sizes) >= num_samples:
                        break  # Stop when num_samples reached

        if insert_sizes:
            threshold = np.percentile(insert_sizes, 99)
            filtered_sizes = [size for size in insert_sizes if size <= threshold][:num_samples]
            if filtered_sizes:
                stats = np.array([np.mean(filtered_sizes), np.median(filtered_sizes), np.std(filtered_sizes)])
                all_stats.append(stats)

    return tuple(np.mean(np.stack(all_stats), axis=0)) if all_stats else (None, None, None)

for idx, row in sample_manifest.iterrows():
    outd="/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams/mosdepth"
    sample = row["sample"]
    bam = row["bam fp"]
    
    ## Q30 coverage
    subprocess.run(f"mosdepth --chrom chr21 --threads 4 --no-per-base -x -Q 30 {outd}/{sample}.stat {bam}", shell=True)
    p = subprocess.run(f"tail -n1 {outd}/{sample}.stat.mosdepth.summary.txt | cut -f4", shell=True, capture_output=True)
    sample_manifest.loc[sample_manifest["sample"] == sample, "Q30 depth"] = float(p.stdout.strip())
    
    mean, median, std = get_insert_dist(bam)
    sample_manifest.loc[sample_manifest["sample"] == sample, "isize mean"] = round(mean, 2)
    sample_manifest.loc[sample_manifest["sample"] == sample, "isize std"] = round(std, 2)

outp = "/lustre1/g/paed_yangwl/snakesch/work/manuscript_figure_table/svhet/tables/bam_statistics.tsv"
sample_manifest[['sample', 'Q30 depth', 'isize mean', 'isize std', 'read length']].to_csv(outp, index=False, sep="\t")