## imports
import os
import joblib

import numpy as np
import pandas as pd
from matplotlib import ticker as ticker

## Suppress annoying pandas warnings
import warnings
warnings.filterwarnings("ignore")

def metrics_by_bin(df: pd.DataFrame) -> pd.DataFrame:

    tp_states = {"tpbase"}

    def _agg(g: pd.DataFrame) -> pd.Series:
        s = g["state"].astype(str)
        tp = s.isin(tp_states).sum()
        fp = (s == "fp").sum()
        fn = (s == "fn").sum()

        precision = tp / (tp + fp) if (tp + fp) else np.nan
        recall    = tp / (tp + fn) if (tp + fn) else np.nan
        f1        = 2 * precision * recall / (precision + recall) if (precision + recall) else np.nan

        return pd.Series({"tp": tp, "fp": fp, "fn": fn,
                          "precision": precision, "recall": recall, "f1": f1})

    labels = ["[50, 100]", "(100, 500]", "(500, 1000]", "(1000, 5000]", "(5000, inf)"]
    df["bin"] = pd.cut(
                        df["svlen"],
                        bins=[50, 100, 500, 1000, 5000, np.inf],
                        labels=labels,
                        right=True,
                        include_lowest=False, 
                        ordered=True
                    )
    
    out = df.groupby("bin", observed=True, sort=True).apply(_agg)
    out.index.name = "bin"
    return out

## We need a function to create a multilevel-indexed df with P/R/F1 for each bin 
def get_metrics_dist(root, samples, compare_set, label=""):
    '''
    Computes min, max and median for all precision, recall and F1 of given comparison set
    
    Args:
        root: comparison VCF path
        samples: list of all samples
        compare_set: one of orig / final     
    '''
    binned_dfs = []
    for sample in samples:
        jl = os.path.join(root, sample, compare_set, 'data.jl')
        df = joblib.load(jl)
        binned = metrics_by_bin(df)
        binned = binned.loc[:, ["precision", "recall", "f1"]]
        binned_dfs.append(binned)

    combined = pd.concat(binned_dfs)
    combined = combined.loc[:, combined.columns.intersection(["precision", "recall", "f1"])]

    ## Only evaluates mean
    out = combined.groupby(combined.index).mean()
    
    ## min and max
    # out = combined.groupby(level=0).agg(["mean"])
    # out = (combined
    #     .groupby(level=0, sort=True, observed=True)
    #     .agg(['mean', 'min', 'max'])
    #     )
    
    out.index.name = 'length'

    ## Re-ordering
    out_stacked = (
        out.rename_axis(index="length", columns="metric")
           .stack()
           .rename(label)
           .to_frame()
    )
    
    return out_stacked

def get_counts(root, samples, compare_set, label=""):
    '''
    Computes min, max and median for all precision, recall and F1 of given comparison set
    
    Args:
        root: comparison VCF path
        samples: list of all samples
        compare_set: one of orig / final     
    '''
    binned_dfs = []
    for sample in samples:
        jl = os.path.join(root, sample, compare_set, 'data.jl')
        df = joblib.load(jl)
        binned = metrics_by_bin(df)
        binned = binned.loc[:, ["tp", "fp", "fn"]]
        binned_dfs.append(binned)

    combined = pd.concat(binned_dfs)
    combined = combined.loc[:, combined.columns.intersection(["tp", "fp", "fn"])]

    ## Only evaluates mean
    out = combined.groupby(combined.index).mean()
    
    ## min and max
    # out = combined.groupby(level=0).agg(["mean"])
    # out = (combined
    #     .groupby(level=0, sort=True, observed=True)
    #     .agg(['mean', 'min', 'max'])
    #     )
    
    out.index.name = 'length'

    ## Re-ordering
    out_stacked = (
        out.rename_axis(index="length", columns="metric")
           .stack()
           .rename(label)
           .to_frame()
    )
    
    return out_stacked

def get_fp_and_recall(root, samples, compare_set, label=""):
    '''
    Computes min, max and median for all precision, recall and F1 of given comparison set
    
    Args:
        root: comparison VCF path
        samples: list of all samples
        compare_set: one of orig / final     
    '''
    binned_dfs = []
    for sample in samples:
        jl = os.path.join(root, sample, compare_set, 'data.jl')
        df = joblib.load(jl)
        binned = metrics_by_bin(df)
        binned = binned.loc[:, ["fp", "recall"]]
        binned.index = pd.MultiIndex.from_product(
            [[sample], binned.index],
            names=["sample", "length"]
        )
        binned_dfs.append(binned)

    combined = pd.concat(binned_dfs)
    combined = combined.loc[:, combined.columns.intersection(["fp", "recall"])]

    return combined

def consolidate_across_tools(baseline, *tools, labels=[]):

    tool_perf = {}
    for label, tool in zip(labels, tools):

        delta = (tool - baseline) / baseline
        summary = (delta
                .groupby(level="length", sort=True, observed=True)
                .agg(['min', 'max', 'mean'])
                )
        summary.columns.names = ['metric', 'stat']

        tool_perf[label] = summary

    return tool_perf
