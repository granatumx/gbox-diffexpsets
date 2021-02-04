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

    numrows = gn.get_arg('numrows')

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

    zscore_dfs = []
    colnames = []
    for coli in mean_df:
        for colj in mean_df:
            if coli != colj:
                zscore_dfs.append(((mean_df[coli]-mean_df[colj])/(std_df[coli]+1.0)).fillna(0).clip(-20.0, 20.0))
                colnames.append("{} vs {}".format(coli, colj)) 
    zscore_df = pd.concat(zscore_dfs, axis=1)
    zscore_df.columns = colnames

    gn.export_statically(gn.assay_from_pandas(zscore_df.T), 'Differential expression sets')

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
