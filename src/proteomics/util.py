def lfq_col(sample):
    if type(sample) == list:
        return list(map(lfq_col, sample))
    if type(sample) == str:
        return f"LFQ intensity {sample}"
    raise Exception(f"wrong type: {sample}")


common_cols = ["uniprotID", "gene", "Potential contaminant", "Reverse"]