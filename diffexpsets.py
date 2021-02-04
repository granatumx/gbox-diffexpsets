#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
from granatum_sdk import Granatum


def main():
    tic = time.perf_counter()

    gn = Granatum()

    assay = gn.pandas_from_assay(gn.get_import('assay'))
    groups = gn.get_import('groups')

    numrows = gn.get_arg('numrows')

    inv_map = {}
    for k, v in groups.items():
        if(!(v in inv_map)): 
            inv_map[v] = []
        inv_map[v] = inv_map[v].append(k)

    mean_dfs = []
    for k, v in inv_map.items():
        mean_dfs.append(df.loc[:, v].mean(axis=1))
    mean_df = pd.concat(mean_dfs, axis=1)

    std_dfs = []
    for k, v in inv_map.items():
        std_dfs.append(df.loc[:, v].std(axis=1))
    std_df = pd.concat(std_dfs, axis=1)

    zscore_dfs = []
    for coli in mean_df:
        for colj in mean_df:
            if coli != colj:
                zscore_dfs.append((mean_df[coli]-mean_df[colj])/std_df[colj])
    zscore_df = pd.concat(zscore_dfs, axis=1)

    gn.export_statically(gn.assay_from_pandas(zscore_df.T), 'Differential expression sets')

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
