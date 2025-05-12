from cyvcf2 import Variant
from pysam import AlignmentFile
from svhet.utils.log import setup_logger
from svhet.read_sets import ReadSet
from svhet.reads import Read

logger = setup_logger("svhet")

UNK = (None, None)
MAX_SV_SIZE = 50
INS = "INS"
DEL = "DEL"

def is_between(pos, start, end):
    '''Check if pos is in [start, end]'''
    return pos >= start and pos <= end

class SV(object):
    def __init__(self, variant: Variant, read_length: int, cipos_tag="CIPOS", ciend_tag="CIEND"):
        self.variant = variant
        self.vartype = variant.INFO.get("SVTYPE")
        self.cipos = variant.INFO.get(cipos_tag)
        self.ciend = variant.INFO.get(ciend_tag)
        
        # Direct assignment for faster access
        self.is_ins = (self.vartype == INS)
        self.is_del = (self.vartype == DEL)
        
        self._read_length = read_length
        
        if not self.vartype:
            raise ValueError(f"Variant {str(variant)} is not a structural variant")
        if self.vartype not in [DEL, INS]:
            raise ValueError(f"Unhandled variant type: {self.vartype}")

        # Calculate breakpoints once at initialization
        self.start_ci, self.end_ci = self.calculate_breakpoints()
        
        # Precompute commonly accessed values
        self._ambiguous_start = self._compute_ambiguous_start()
        self._ambiguous_end = self._compute_ambiguous_end()
        self._variant_size = self._compute_variant_size()

    def _compute_ambiguous_start(self):
        if self.cipos is None:
            return True
        return abs(self.cipos[0]) + abs(self.cipos[1]) >= self._read_length
        
    def _compute_ambiguous_end(self):
        end = self.ciend or self.cipos
        if end is None:
            return True
        return abs(end[0]) + abs(end[1]) >= self._read_length
        
    def _compute_variant_size(self):
        if self.is_ins:
            size = self.variant.INFO.get("SVLEN") or self.variant.INFO.get("INSLEN")
            return abs(size) if size is not None else 0
        elif self.is_del:
            size = self.variant.INFO.get("SVLEN")
            if size is None and self.variant.end and self.variant.POS:
                size = self.variant.end - self.variant.POS
                return abs(size) if size is not None else 0
            else:
                return abs(size) if size is not None else 0
        else:
            return 0

    @property
    def varsize(self) -> int:
        return self._variant_size
        
    @property
    def ambiguous_start(self):
        return self._ambiguous_start
        
    @property
    def ambiguous_end(self):
        return self._ambiguous_end
        
    def calculate_breakpoints(self) -> tuple:
        '''Calculate breakpoint coordinates with confidence intervals'''
        # Fast exit for ambiguous breakpoints
        if self._compute_ambiguous_start() and self._compute_ambiguous_end():
            return UNK, UNK
        
        start = UNK
        end = UNK
        
        # Cache frequently accessed values
        variant = self.variant
        cipos = self.cipos
        ciend = self.ciend
        
        if not self._compute_ambiguous_start():
            start = variant.POS - abs(cipos[0]), variant.POS + abs(cipos[1])
            
        if not self._compute_ambiguous_end():
            if ciend is not None and variant.end is not None:
                end = variant.end - abs(ciend[0]), variant.end + abs(ciend[1])
            elif variant.end is not None:
                end = variant.end - abs(cipos[0]), variant.end + abs(cipos[1])
            elif self.is_ins:
                end = start
            else:
                raise ValueError(f"No END for DEL variant:\n{self}\n")
        
        return start, end
    
    def get_wt_haplotype_evidence(self, bam: AlignmentFile, mapping_quality=30):
        # Fast exit if not a deletion
        if not self.is_del:
            return set()
        
        rg1, rg2, rg3 = [], [], []
        variant = self.variant  # Cache for repeated access
        
        # Fetch reads for first region
        if not self.ambiguous_start:
            start_ci = self.start_ci  # Cache for repeated access
            for read in bam.fetch(variant.CHROM, start_ci[0], start_ci[1]):
                if read.mapping_quality < mapping_quality:
                    continue
                elif read.reference_start <= start_ci[0] and read.reference_end >= start_ci[1]:
                    rg1.append(read)

        # Fetch reads for second region
        if not self.ambiguous_end:
            end_ci = self.end_ci  # Cache for repeated access
            for read in bam.fetch(variant.CHROM, end_ci[0], end_ci[1]):
                if read.mapping_quality < mapping_quality:
                    continue
                elif read.reference_start <= end_ci[0] and read.reference_end >= end_ci[1]:
                    rg2.append(read)
            
        # Fetch reads for middle region
        if self.cipos is not None:
            cipos = self.cipos  # Cache for repeated access
            ciend = self.ciend  # Cache for repeated access
            start = variant.POS + abs(cipos[1])
            end = variant.end - abs(ciend[0]) if ciend is not None else variant.end - abs(cipos[0])
            if start < end:
                for read in bam.fetch(variant.CHROM, start, end):
                    if read.mapping_quality < mapping_quality:
                        continue
                    if read.reference_start > start and read.reference_end < end:
                        rg3.append(read) 
          
        # Process each read group with cached coordinates
        start_ci = self.start_ci
        end_ci = self.end_ci
        
        rs1 = ReadSet(rg1)
        rs1.remove_clipped_reads(start_ci[0], end_ci[1])
        rs1.remove_reads_with_small_indels(start_ci[1], end_ci[0], low=0.7, high=None)
        rs1.remove_high_NM_reads(start_ci[0], end_ci[1])
            
        rs2 = ReadSet(rg2)
        rs2.remove_clipped_reads(start_ci[0], end_ci[1])
        rs2.remove_reads_with_small_indels(start_ci[1], end_ci[0], low=0.7, high=None)
        rs2.remove_high_NM_reads(start_ci[0], end_ci[1])
 
        rs3 = ReadSet(rg3)
        rs3.remove_high_NM_reads(start_ci[0], end_ci[1])
        
        # Sample reads based on coverage
        rs1.read_sampler(start_ci, end_ci)
        rs2.read_sampler(start_ci, end_ci)
        
        # Cluster reads by alignment score
        rs1.get_major_read_cluster()
        rs2.get_major_read_cluster()
        rs3.get_major_read_cluster()
        
        # Directly combine sets
        return set(rs1.get_reads()) | set(rs2.get_reads()) | set(rs3.get_reads())
    
    def get_mut_haplotype_evidence(self, bam: AlignmentFile, mapping_quality=30):
        # Fast exits for common cases
        if self._read_length < 500 and not self.is_del:
            return set()
        
        if self.ambiguous_start or self.ambiguous_end:
            return set()
                
        rg1 = []
        rg2 = []
        rg3 = []
        
        # Cache frequently accessed values
        variant = self.variant
        start_ci = self.start_ci
        end_ci = self.end_ci
        
        # Process reads in a single pass
        for pysam_read in bam.fetch(variant.CHROM, start_ci[0], end_ci[1]):
            if pysam_read.mate_is_unmapped:
                continue
          
            read = Read(pysam_read)

            if read.mapping_quality < mapping_quality:
                continue

            ref_start = read.reference_start
            ref_end = read.reference_end
            
            # Case 1: Read spanning 5' breakpoint
            if (is_between(ref_end, start_ci[0] - 1, start_ci[1]) and 
                not is_between(ref_end, end_ci[0] - 1, end_ci[1])):
                
                if read.is_right_clipped:
                    rg1.append(read)
                elif (read.has_del_in_region(start_ci[0], end_ci[1]) and 
                      read.cigartuples[-1][0] == 2):
                    rg1.append(read)
                    
            # Case 2: Read spanning 3' breakpoint
            elif (is_between(ref_start, end_ci[0] - 1, end_ci[1]) and 
                  not is_between(ref_end, end_ci[0] - 1, end_ci[1])):
                
                if read.is_left_clipped:
                    rg2.append(read)
                elif (read.has_del_in_region(start_ci[0], end_ci[1]) and 
                      read.cigartuples[0][0] == 2):
                    rg2.append(read)
                    
            # Case 3: Read spanning both 5' and 3' breakpoints
            elif (read.has_del_in_region(start_ci[0], end_ci[1]) and 
                  ref_start <= start_ci[1] and 
                  ref_end >= end_ci[0]):
                
                # Choose based on which breakpoint the read flanks more
                start_dist = start_ci[0] - ref_start
                end_dist = ref_end - end_ci[1]
                
                if end_dist < start_dist and read.is_right_clipped:
                    continue
                elif end_dist > start_dist and read.is_left_clipped:
                    continue
                    
                rg3.append(read)
        
        # Process read groups
        rs1 = ReadSet(rg1)
        rs1.get_clipped_reads(start_ci[1], end_ci[0])
        rs1.remove_high_NM_reads(start_ci[0], start_ci[1])
        rs1.remove_reads_with_small_indels(start_ci[0], end_ci[1], high=0.3)
        cleaned_rg1 = set(rs1.get_reads())
            
        rs2 = ReadSet(rg2)
        rs2.get_clipped_reads(start_ci[1], end_ci[0])
        rs2.remove_high_NM_reads(end_ci[0], end_ci[1])
        rs2.remove_reads_with_small_indels(start_ci[0], end_ci[1], high=0.3)
        cleaned_rg2 = set(rs2.get_reads())
        
        # Use direct filtering for third group
        cleaned_rg3 = set(
            r for r in rg3 
            if not r.has_del_in_region(start_ci[0], end_ci[1], high=0.3)
        )
        
        # Return union of all cleaned read groups
        return cleaned_rg1 | cleaned_rg2 | cleaned_rg3

    def get_hets_within_deletion(self, hets):
        '''Extract HETs relevant to candidate HDEL region.'''
        # Fast exit for ambiguous breakpoints
        if self.ambiguous_start and self.ambiguous_end:
            return []
        
        # Cache breakpoint coordinates
        variant = self.variant
        read_length = self._read_length
        
        # Calculate the region to search
        start = self.start_ci[1] if not self.ambiguous_start else variant.POS + read_length
        end = self.end_ci[0] if not self.ambiguous_end else variant.end - read_length
        
        # Fast exit if invalid region
        if start >= end:
            return []
        
        # Get relevant HET sites within confidence interval
        het_chrom = variant.CHROM
        
        relevant_hets = [
            het for het in hets 
            if het[0] == het_chrom and start < het[1] < end
        ]
        
        assert isinstance(relevant_hets, list), f"Incorrect dtype for relevant_hets {type(relevant_hets)}"
        return relevant_hets

    def associated_hets(self, hets, distance=5):
        '''Find heterozygous positions associated with this SV.'''
        # Fast exit for ambiguous breakpoints or small variants
        if self.ambiguous_start or self.ambiguous_end or self.varsize < MAX_SV_SIZE:
            return []
        
        # Cache values for faster access
        variant = self.variant
        read_length = self._read_length
        start_ci = self.start_ci
        end_ci = self.end_ci
        
        beg = start_ci[0] - read_length * distance
        end = end_ci[1] + read_length * distance
        
        het_chrom = variant.CHROM
        associated_hets = [
            het[1] for het in hets
            if (het[0] == het_chrom and 
                beg <= het[1] <= end and 
                (het[1] <= start_ci[0] or het[1] >= end_ci[1]))
        ]
                
        return associated_hets
        
    def __repr__(self):
        '''Returns string representation of expanded event interval'''
        return f"{self.variant.CHROM}:{self.start_ci[0]}-{self.end_ci[1]}"
