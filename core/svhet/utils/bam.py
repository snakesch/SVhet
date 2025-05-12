import pysam

def get_read_length(bam: str) -> int:
    import random
    """Estimates average read length from BAM file"""
    with pysam.AlignmentFile(bam) as bamfs:
        rls = []
        for aln in bamfs:
            if random.random() > 0.02: 
                continue
            if aln.query_length > 0:
                rls.append(aln.query_length)
            if len(rls) >= 1000: break
    return round(sum(rls) / len(rls))