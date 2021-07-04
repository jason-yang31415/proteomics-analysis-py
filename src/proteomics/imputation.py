import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
from .util import lfq_col


def impute_3_1(data):
    """
    impute missing values with small values drawn from a uniform distribution
    between replicate mean - 3 * replicate std and replicate mean - 2 *
    replicate std
    """

    def impute_col(col):
        mean, std = col.mean(), col.std()
        return col.apply(
            lambda val: val
            if pd.notna(val)
            else stats.uniform.rvs(loc=mean - 3 * std, scale=std)
        )

    return data.apply(impute_col, axis=0)


def impute_perseus(data):
    def impute_col(col):
        mean, std = col.mean(), col.std()
        return col.apply(
            lambda val: val
            if pd.notna(val)
            else stats.norm.rvs(loc=mean - 1.8 * std, scale=0.3 * std)
        )

    return data.apply(impute_col, axis=0)


def impute_4_6(data, n_rep=2):
    """
    impute missing values relatively by looking at deltas between the replicate
    with the missing value and a reference replicate where the protein was
    detected; skip proteins detected in less than `n_rep` replicates
    """
    # cache delta distributions since there are only ~N_r^2 possible
    # combinations
    delta_cache = {}

    def delta(col, refcol):
        if (col, refcol) in delta_cache:
            return delta_cache[(col, refcol)]
        # calculate column of deltas betweel col and refcol
        delta = (data[col] - data[refcol]) / data[[col, refcol]].mean(axis=1)
        # compute stats for deltas
        mean, std = delta.mean(), delta.std()
        # put result in cache then return result
        delta_cache[(col, refcol)] = mean, std
        return mean, std

    # calculate correlation between replicates
    corr_mat = data.corr()
    # cache mean correlations
    corr_cache = {}

    def corr(refcol, imputecols):
        if (refcol, frozenset(imputecols)) in corr_cache:
            return corr_cache[(refcol, frozenset(imputecols))]
        # calculate mean of correlations between refcol and each of imputecols
        corr_mean = corr_mat.loc[refcol, imputecols].mean()
        # put mean correlation in cache then return result
        corr_cache[(refcol, frozenset(imputecols))] = corr_mean
        return corr_mean

    def impute_row(row):
        # get list of columns to use as reference
        refcols = [c for c in data.columns if pd.notna(row[c])]
        # get list of columns to impute
        imputecols = [c for c in data.columns if pd.isna(row[c])]
        for col in imputecols:
            Inews = []  # list of imputed values based on each refcol
            for refcol in refcols:
                corr_mean = corr(refcol, imputecols)
                Dmean, Dstd = delta(col, refcol)
                Dnew = stats.norm.rvs(loc=Dmean, scale=Dstd / ((2 ** 0.5) * corr_mean))
                Inews.append(row[refcol] * abs(1 + Dnew))
            row[col] = np.array(Inews).mean()
        return row

    cp = data.copy(True)
    # only run imputation on rows where num non-zero replicates >= n_rep
    mask = pd.notna(data).sum(axis=1) >= n_rep
    # run imputation on each row then return result
    cp.loc[mask] = cp.loc[mask].apply(impute_row, axis=1)
    return cp


def impute_twostep(data, imputed_small):
    cp = data.copy()
    mask = pd.isna(data).any(axis=1)
    cp[mask] = imputed_small[mask]
    return cp


def compare(comparisons, conditions, unimputed, imputed_list):
    data_comparisons = {key: None for key in comparisons}

    for comparison in comparisons:
        dfs = []  # list of comparison statistics for each imputation run
        for df_imputed in imputed_list:
            df = pd.DataFrame.from_dict(
                {
                    # post-imputation mean intensity for each condition
                    f"mean {comparison[0]}": df_imputed[
                        lfq_col(conditions[comparison[0]])
                    ].mean(axis=1),
                    f"mean {comparison[1]}": df_imputed[
                        lfq_col(conditions[comparison[1]])
                    ].mean(axis=1),
                    # calculate p value using t test
                    "log p": np.log(
                        stats.ttest_ind(
                            df_imputed[lfq_col(conditions[comparison[0]])],
                            df_imputed[lfq_col(conditions[comparison[1]])],
                            axis=1,
                            equal_var=False,
                        ).pvalue
                    ),
                }
            )
            dfs.append(df)
        # concatenate comparison statistic tables and average across imputation
        # runs
        df_concat = pd.concat(dfs)
        df_comparison = df_concat.groupby(df_concat.index).mean()
        # calculate p value from averaged log p value
        df_comparison["p"] = np.exp(df_comparison["log p"])
        # caluclate log FC from averaged mean intensity
        df_comparison["log FC"] = (
            df[f"mean {comparison[1]}"] - df[f"mean {comparison[0]}"]
        )
        # copy over gene column
        df_comparison["gene"] = unimputed["gene"]
        # remove proteins not detected in either condition before adjusting p
        # values
        mask = pd.notna(
            unimputed[lfq_col(conditions[comparison[0]] + conditions[comparison[1]])]
        ).any(axis=1)
        df_comparison = df_comparison.loc[mask]
        # adjust p values using p values averaged across imputation runs
        _, padjusted, _, _ = multitest.multipletests(
            df_comparison["p"].values, method="fdr_bh"
        )
        df_comparison["p adjusted"] = padjusted
        # count number of replicates protein was detected in for each condition
        df_comparison[f"n {comparison[0]}"] = pd.notna(
            unimputed[lfq_col(conditions[comparison[0]])].loc[mask]
        ).sum(axis=1)
        df_comparison[f"n {comparison[1]}"] = pd.notna(
            unimputed[lfq_col(conditions[comparison[1]])].loc[mask]
        ).sum(axis=1)
        # mark proteins as significant
        df_comparison["significant"] = (
            (df_comparison["p adjusted"] < 0.05)
            & (abs(df_comparison["log FC"]) > 1)
            & (
                (df_comparison[f"n {comparison[0]}"] >= 2)
                | (df_comparison[f"n {comparison[1]}"] >= 2)
            )
        )

        data_comparisons[comparison] = df_comparison

    return data_comparisons