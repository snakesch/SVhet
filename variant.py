#!/usr/bin/env python3

import pandas as pd

class StructuralVariant(object):
    
    def __init__(self, row):
        
        self.chrom = str(row["#CHROM"])
        self.start = int(row["POS"])
        self.end = int(row["END"])
        self.het_cnt = 0
        
        self.het_cnt_3prime = 0
        self.het_cnt_5prime = 0
        self._record = row
        
    def __len__(self):
        return self.end - self.start
    
    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end} ({len(self)}bp; het sites: {self.het_cnt} + {self.het_cnt_3prime} (3prime) + {self.het_cnt_5prime} (5prime))"
    
    def __eq__(self, other):
        
        tolerance = 100 if len(self) >= 200 else 10 # bp tolerance for coordinate difference

        if self.chrom != other.chrom:
            return False
        if abs(self.start - other.start) < tolerance and abs(self.end - other.end) < tolerance:
            return True
        return False
    
    def to_dataframe(self):
        """
        Convert the StructuralVariant object into a pandas DataFrame.
        
        Returns:
        pandas.DataFrame: A DataFrame containing the properties of the StructuralVariant object.
        """
        data = {
            "chrom": [self.chrom],
            "start": [self.start],
            "end": [self.end],
            "length": [len(self)],
            "het_cnt": [self.het_cnt]
        }
        return pd.DataFrame(data)
