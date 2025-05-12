from collections import defaultdict
import numpy as np
import random
from functools import lru_cache
from svhet.reads import Read
from svhet.utils.log import setup_logger

logger = setup_logger("svhet")

class ReadSet:
    __slots__ = [
        'read_set', 'leftmost_pos', 'rightmost_pos',
        '_ref_starts', '_ref_ends', '_is_left_clipped', '_is_right_clipped', '_cigartuples',
        '_empty'
    ]

    def __init__(self, reads: list[Read]):
        # Fast path for empty reads
        if not reads:
            self.read_set = []
            self.leftmost_pos = self.rightmost_pos = 0
            self._ref_starts = np.array([], dtype=np.int32)  # Specify dtype
            self._ref_ends = np.array([], dtype=np.int32)
            self._is_left_clipped = np.array([], dtype=np.bool_)
            self._is_right_clipped = np.array([], dtype=np.bool_)
            self._cigartuples = []
            self._empty = True
            return

        self._empty = False
        self.read_set = [r if isinstance(r, Read) else Read(r) for r in reads]
        
        self._ref_starts = np.array([r.reference_start for r in self.read_set], dtype=np.int32)
        self._ref_ends = np.array([r.reference_end for r in self.read_set], dtype=np.int32)
        
        self.leftmost_pos = self._ref_starts.min() if len(self._ref_starts) > 0 else 0
        self.rightmost_pos = self._ref_ends.max() if len(self._ref_ends) > 0 else 0
        
        self._cigartuples = [r.cigartuples for r in self.read_set]
        
        self._is_left_clipped = np.array([
            bool(c and c[0][0] in {4, 5}) for c in self._cigartuples
        ], dtype=np.bool_)
        
        self._is_right_clipped = np.array([
            bool(c and c[-1][0] in {4, 5}) for c in self._cigartuples
        ], dtype=np.bool_)

    def _update_precomputed(self):
        # Skip if empty
        if not self.read_set:
            self._empty = True
            self._ref_starts = np.array([], dtype=np.int32)
            self._ref_ends = np.array([], dtype=np.int32)
            self._is_left_clipped = np.array([], dtype=np.bool_)
            self._is_right_clipped = np.array([], dtype=np.bool_)
            self._cigartuples = []
            self.leftmost_pos = self.rightmost_pos = 0
            return
            
        self._empty = False
        
        self._ref_starts = np.array([r.reference_start for r in self.read_set], dtype=np.int32)
        self._ref_ends = np.array([r.reference_end for r in self.read_set], dtype=np.int32)
        self._cigartuples = [r.cigartuples for r in self.read_set]
        
        # Update min/max positions
        if len(self._ref_starts) > 0:
            self.leftmost_pos = self._ref_starts.min()
            self.rightmost_pos = self._ref_ends.max()
        else:
            self.leftmost_pos = self.rightmost_pos = 0
            
        # Recalculate clipped flags
        self._is_left_clipped = np.array([
            bool(c and c[0][0] in {4, 5}) for c in self._cigartuples
        ], dtype=np.bool_)
        
        self._is_right_clipped = np.array([
            bool(c and c[-1][0] in {4, 5}) for c in self._cigartuples
        ], dtype=np.bool_)

    def get_reads(self):
        return self.read_set

    def remove_clipped_reads(self, var_start, var_end):
        # Early exit for empty sets or invalid range
        if var_start is None or var_end is None or not self.read_set:
            return
            
        mask = ~(
            ((self._ref_ends < var_end) & self._is_right_clipped) |
            ((self._ref_starts > var_start) & self._is_left_clipped)
        )
        
        self._apply_mask(mask)

    def get_clipped_reads(self, var_start, var_end):
        # Early exit for empty sets or invalid range
        if var_start is None or var_end is None or not self.read_set:
            return []
            
        right_clip_len = np.array([
            c[-1][1] if (c and c[-1][0] in {4, 5}) else 0 for c in self._cigartuples
        ], dtype=np.int32)
        
        right_clip = (
            (self._ref_ends < var_end) & 
            self._is_right_clipped & 
            ((self._ref_ends + right_clip_len) > var_start)
        )
        
        left_clip_len = np.array([
            c[0][1] if (c and c[0][0] in {4, 5}) else 0 for c in self._cigartuples
        ], dtype=np.int32)
        
        left_clip = (
            (self._ref_starts > var_start) & 
            self._is_left_clipped & 
            ((self._ref_starts - left_clip_len) < var_end)
        )
        
        mask = right_clip | left_clip
        
        self._apply_mask(mask)
        return self.read_set

    def _apply_mask(self, mask):
        # Fast path if no reads pass the filter
        if not np.any(mask):
            self.read_set = []
            self._empty = True
            self._ref_starts = np.array([], dtype=np.int32)
            self._ref_ends = np.array([], dtype=np.int32)
            self._is_left_clipped = np.array([], dtype=np.bool_)
            self._is_right_clipped = np.array([], dtype=np.bool_)
            self._cigartuples = []
            self.leftmost_pos = self.rightmost_pos = 0
            return
            
        self.read_set = [r for r, m in zip(self.read_set, mask) if m]
        self._ref_starts = self._ref_starts[mask]
        self._ref_ends = self._ref_ends[mask]
        self._is_left_clipped = self._is_left_clipped[mask]
        self._is_right_clipped = self._is_right_clipped[mask]
        self._cigartuples = [c for c, m in zip(self._cigartuples, mask) if m]
        
        # Update min/max positions if not empty
        if len(self._ref_starts) > 0:
            self.leftmost_pos = self._ref_starts.min()
            self.rightmost_pos = self._ref_ends.max()
        else:
            self.leftmost_pos = self.rightmost_pos = 0
            self._empty = True

    @lru_cache(maxsize=8)
    def has_del_in_region(self, var_start, var_end, low=None, high=None):
        return any(
            r.has_del_in_region(var_start, var_end, low=low, high=high) 
            for r in self.read_set
        )

    def remove_reads_with_small_indels(self, var_start, var_end, low=None, high=0.3):
        # Early exit for empty sets
        if not self.read_set:
            return
            
        mask = np.array([
            not r.has_del_in_region(var_start, var_end, low=low, high=high) 
            for r in self.read_set
        ], dtype=np.bool_)
        
        self._apply_mask(mask)

    def remove_high_NM_reads(self, var_start, var_end):
        # Early exit for empty sets
        if not self.read_set:
            return
            
        mask = np.array([
            (r.has_del_in_region(var_start, var_end, low=0.5, high=0.8) or 
             (r.nm is not None and r.nm <= 5))
            for r in self.read_set
        ], dtype=np.bool_)
        
        self._apply_mask(mask)

    def read_sampler(self, start_ci, end_ci):
        # Early exits for common cases
        if (not self.read_set or 
            start_ci == (None, None) or 
            end_ci == (None, None) or 
            start_ci[1] > end_ci[0]):
            return []
            
        reads = self.read_set
        sampled_reads = []
        covs = []
        
        # First pass - collect reads to process
        for read in reads:
            if read.is_from_3prime(start_ci, end_ci):
                cov = max(0, end_ci[0] - read.reference_start)
            elif read.is_from_5prime(start_ci, end_ci):
                cov = max(0, read.reference_end - start_ci[1])
            elif read.is_within_deletion(start_ci, end_ci):
                continue
            else:
                sampled_reads.append(read)
                continue
            covs.append(cov)
            
        # Fast path if no coverage calculation needed
        if not covs and sampled_reads:
            self.read_set = sampled_reads
            self._update_precomputed()
            return sampled_reads
            
        # Fast path if no covered reads
        if not covs or max(covs) == 0:
            self.read_set = []
            self._update_precomputed()
            return []
            
        # Calculate weights for sampling
        max_cov = max(covs)
        read2cov = [(read, cov / max_cov) for read, cov in zip(reads, covs)]
        total_cov = sum(cov for _, cov in read2cov)
        
        # Perform weighted sampling
        for read, cov in read2cov:
            if random.random() <= cov / total_cov:
                sampled_reads.append(read)
                
        self.read_set = sampled_reads
        self._update_precomputed()
        return sampled_reads

    def get_major_read_cluster(self):
        # Early exit for empty sets
        if not self.read_set:
            return []

        chunk_size = 150
        chunk_dict = defaultdict(list)
        
        for read in self.read_set:
            chunk_idx = (read.reference_start - self.leftmost_pos) // chunk_size
            chunk_dict[chunk_idx].append(read)
            
        valid_reads = []
        
        for reads in chunk_dict.values():
            # Fast path for small chunks or no NM variations
            if len(reads) < 5:
                valid_reads.extend(reads)
                continue
            elif all(r.nm == 0 for r in reads):
                continue
                
            scores = np.array([r.as_ for r in reads])
            mean, sd = np.mean(scores), np.std(scores)
            
            # Fast path for no variation
            if sd == 0:
                valid_reads.extend(reads)
            else:
                valid_reads.extend([
                    read for read in reads 
                    if mean - 2*sd <= read.as_ <= mean + 2*sd
                ])
                
        self.read_set = valid_reads
        self._update_precomputed()
        return self.read_set
