import pandas as pd
import math
import seaborn as sns
from .util import lfq_col


def import_maxquant(path, normalize=False, plot=True):
    # load proteinGroups file as dataframe
    data_raw = pd.read_csv(path, sep="\t", header=0, index_col="id")
    # add columns for uniprotID and gene; cast potential contaminant and
    # reverse to booleans
    data_raw["uniprotID"] = data_raw["Protein IDs"].apply(
        lambda val: val.split(";")[0])
    data_raw["gene"] = data_raw["Gene names"].apply(
        lambda val: val.split(";")[0] if pd.notna(val) else val
    )
    data_raw["Potential contaminant"] = data_raw["Potential contaminant"].apply(
        lambda val: val == "+"
    )
    data_raw["Reverse"] = data_raw["Reverse"].apply(lambda val: val == "+")
    # get names of samples from columns
    samples = [
        sample.split("LFQ intensity ")[1]
        for sample in data_raw.columns.to_list()
        if sample.startswith("LFQ intensity")
    ]
    # filter out potential contaminants and reverse
    data_filtered = data_raw[
        (data_raw["Potential contaminant"] == False) & (
            data_raw["Reverse"] == False)
    ]
    # log transform lfq intensities
    data_log = data_filtered.copy(False)
    data_log[lfq_col(samples)] = data_log[lfq_col(samples)].apply(
        lambda col: col.apply(lambda val: math.log2(val) if val > 0 else None)
    )
    # filter out proteins with NA in all samples
    data_notna = data_log[
        data_log[lfq_col(samples)].apply(
            lambda row: not row.apply(pd.isna).all(), axis=1
        )
    ]

    if normalize:
        # normalize sample medians
        data_normalized = data_notna.copy(False)
        max_median = data_notna[lfq_col(samples)].median(axis=0).max()
        data_normalized[lfq_col(samples)] = data_normalized[lfq_col(samples)].apply(
            lambda col: col * max_median / col.median(), axis=0
        )
    else:
        data_normalized = data_notna

    if plot:
        sns.boxplot(data=data_normalized[lfq_col(samples)], orient="h").set(
            xlabel="log2 LFQ intensity"
        )
    return data_normalized, samples
