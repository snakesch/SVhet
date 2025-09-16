#!/bin/bash
#SBATCH --job-name=bench-time                   # 1. Job name
#SBATCH --mail-type=FAIL                        # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=snakesch@connect.hku.hk     #    Email address to receive notification
#SBATCH --partition=intel                       # 3. Request a partition
#SBATCH --mem=50G                               #    CPU memory
#SBATCH --qos=normal                            # 4. Request a QoS
#SBATCH --ntasks-per-node=5                     # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                               # 6. Request number of node(s)
#SBATCH --time=0-10:00:00                       # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/log/%x_%j.out                          # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/log/%x_%j.err                           #    Standard error log as $job_name_$job_id.err

## HG00096 has 2688 heterozygous deletions in total

## Time profiled on: 500, 1000, 1500, 2000

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

BINS=( 500 1000 1500 2000 2500 2688 )
SV_VCF=/lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/sv/lumpy/final_lumpy_gt-smoove.genotyped.vcf.gz
HEADER_LINES=$(bcftools view -h $SV_VCF | wc -l)

OUT=/lustre1/g/paed_yangwl/snakesch/tmp/test_svhet/time.txt

[[ -f $OUT ]] && rm -f $OUT; 

touch $OUT;

for bin in ${BINS[@]}
do
    for nthr in 1 2 4 8
    do
        bcftools filter -s HG00096 -i '(SVTYPE="DEL" & FMT/GT[0] = "het")' $SV_VCF | \
            head -n $((HEADER_LINES + bin)) | \
            bcftools view -Wtbi -Oz -o $tmpdir/HG00096-$bin.vcf.gz

        ## Wall clock time, Peak RSS (kb)  -o $OUT 
        /usr/bin/time -f "$nthr\t$bin\t%E\t%M" -a -o $OUT /lustre1/g/paed_yangwl/snakesch/work/SVhet/core/svhet/svhet.sh \
                --ref /group/paed_yangwl/t2t-chm13/chm13.v2.fasta \
                --sv-vcf $tmpdir/HG00096-$bin.vcf.gz \
                --outdir $tmpdir/$bin-$nthr \
                --manifest /lustre1/g/paed_yangwl/snakesch/data/svhet2/1kg/aligned_bams/test_manifest.txt \
                --jobs $nthr \
                --min-dp 5 \
                --high-het 1
    done
done