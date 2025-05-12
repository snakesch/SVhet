import pysam
from functools import lru_cache

class Read:
    __slots__ = (
        'read', 'query_name', 'reference_id', 'reference_start', 'reference_end',
        'mapping_quality', 'cigartuples', 'query_sequence',
        'nm', 'as_', '_has_ins', '_has_del', '_is_left_clipped',
        '_is_right_clipped', '_deletions', '_hash'
    )

    def __init__(self, read: pysam.AlignedSegment):
        self.read = read
        self.query_name = read.query_name
        self.reference_id = read.reference_id
        self.reference_start = read.reference_start
        self.reference_end = read.reference_end
        self.mapping_quality = read.mapping_quality
        self.cigartuples = read.cigartuples or ()
        self.query_sequence = read.query_sequence
        self.nm = read.get_tag("NM") if read.has_tag("NM") else None
        self.as_ = read.get_tag("AS") if read.has_tag("AS") else None

        self._precompute_flags()
        self._precompute_deletions()
        
        self._hash = hash((self.query_name, self.reference_start))

    def _precompute_flags(self):
        self._has_ins = False
        self._has_del = False
        
        if not self.cigartuples:
            self._is_left_clipped = False
            self._is_right_clipped = False
            return
            
        self._is_left_clipped = self.cigartuples[0][0] in {4, 5}
        self._is_right_clipped = self.cigartuples[-1][0] in {4, 5}
        
        for op, _ in self.cigartuples:
            if op == 1:
                self._has_ins = True
            elif op == 2:
                self._has_del = True
            if self._has_ins and self._has_del:
                break  # Exit early if both found

    def _precompute_deletions(self):
        self._deletions = []
        if not self._has_del:  # Skip if no deletions
            return
            
        pos = self.reference_start
        for op, length in self.cigartuples:
            if op == 2:  # Deletion
                self._deletions.append((pos, length))
            if op in (0, 2, 3, 7, 8):  # Ops that consume reference
                pos += length

    def __hash__(self):
        return self._hash  # Use precomputed hash

    def __eq__(self, other):
        if not isinstance(other, Read):
            return False
        return (
            self.query_name == other.query_name and
            self.reference_start == other.reference_start and
            self.query_sequence == other.query_sequence
        )

    @property
    def has_ins(self) -> bool:
        return self._has_ins

    @property
    def has_del(self) -> bool:
        return self._has_del

    @property
    def is_left_clipped(self) -> bool:
        return self._is_left_clipped

    @property
    def is_right_clipped(self) -> bool:
        return self._is_right_clipped

    @property
    def is_clipped(self) -> bool:
        return self._is_left_clipped or self._is_right_clipped

    @lru_cache(maxsize=32)
    def is_from_5prime(self, start_ci, end_ci) -> bool:
        if start_ci[0] is None or end_ci[0] is None:
            return False
        return self.reference_end < end_ci[0] and self.reference_end >= start_ci[0]

    @lru_cache(maxsize=32)
    def is_from_3prime(self, start_ci, end_ci) -> bool:
        if start_ci[0] is None or end_ci[0] is None:
            return False
        return self.reference_start > start_ci[1] and self.reference_start <= end_ci[0]

    @lru_cache(maxsize=32)
    def is_within_deletion(self, start_ci, end_ci) -> bool:
        if start_ci[0] is None or end_ci[0] is None:
            return False
        # Avoid repeated function calls
        is_from_3 = self.is_from_3prime(start_ci, end_ci)
        is_from_5 = self.is_from_5prime(start_ci, end_ci)
        if is_from_3 or is_from_5:
            return False
        return self.reference_start > start_ci[1] and self.reference_end < end_ci[0]

    @lru_cache(maxsize=64)
    def has_del_in_region(self, var_start, var_end, low=None, high=None):
        if var_start is None or var_end is None:
            return False
            
        # No need to check if there are no deletions
        if not self._has_del or not self._deletions:
            return False

        var_size = max(1, var_end - var_start)

        for del_start, del_length in self._deletions:
            del_end = del_start + del_length
            
            # Skip early if outside region
            if del_end <= var_start or del_start >= var_end:
                continue

            overlap_start = max(del_start, var_start)
            overlap_end = min(del_end, var_end)

            if overlap_start >= overlap_end:
                continue  # No overlap

            if low is None and high is None:
                return True  # Any overlap meets criteria

            del_size_ratio = del_length / var_size
            
            if (low is None or del_size_ratio > low) and (high is None or del_size_ratio < high):
                return True

        return False
