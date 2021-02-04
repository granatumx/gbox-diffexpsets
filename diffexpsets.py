#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
import math
from granatum_sdk import Granatum


def main():
    tic = time.perf_counter()

    gn = Granatum()

    assay = gn.pandas_from_assay(gn.get_import('assay'))
    groups = gn.get_import('groups')

    min_zscore = gn.get_arg('min_zscore')
    max_zscore = gn.get_arg('max_zscore')

    inv_map = {}
    for k, v in groups.items():
        inv_map[v] = inv_map.get(v, []) + [k]

    mean_dfs = []
    std_dfs = []
    for k, v in inv_map.items():
        mean_dfs.append(assay.loc[:, v].mean(axis=1))
        std_dfs.append(assay.loc[:, v].std(axis=1))

    mean_df = pd.concat(mean_dfs, axis=1)
    std_df = pd.concat(std_dfs, axis=1)

    mean_rest_dfs = []
    std_rest_dfs = []
    for k, v in inv_map.items():
        rest_v = list(set(assay.index.values).difference(set(v)))
        mean_rest_dfs.append(assay.loc[:, rest_v].mean(axis=1))
        std_rest_dfs.append(assay.loc[:, rest_v].std(axis=1))
    mean_rest_df = pd.concat(mean_rest_dfs, axis=1)
    std_rest_df = pd.concat(std_rest_dfs, axis=1)

    zscore_dfs = []
    colnames = []
    for coli in mean_df:
        for colj in mean_df:
            if coli != colj:
                zscore_dfs.append(((mean_df[coli]-mean_df[colj])/(std_df[colj]+1.0/max_zscore)).fillna(0).clip(-max_zscore, max_zscore))
                colnames.append("{} vs {}".format(coli, colj)) 
    for coli in mean_df:
        zscore_dfs.append(((mean_df[coli]-mean_rest_df[colj])/(std_rest_df[colj]+1.0/max_zscore)).fillna(0).clip(-max_zscore, max_zscore))
        colnames.append("{} vs rest".format(coli)) 

    zscore_df = pd.concat(zscore_dfs, axis=1)
    zscore_df.columns = colnames
    norms_df = zscore_df.apply(np.linalg.norm, axis=1)
    colsmatching = norms_df.T[(norms_df.T >= min_zscore)].index.values
    return_df = zscore_df.T[colsmatching]
    gn.export_statically(gn.assay_from_pandas(return_df), 'Differential expression sets')
    gn.export(return_df.T.to_csv(), 'differential_gene_sets.csv', kind='raw', meta=None, raw=True)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
