"""For Figure 5"""

import sys

import numpy as np
import pandas as pd

from polymodel.config import Config
from polymodel.consts import (DEFAULT_P, FUNG_MUTATION_SCALE,
                              HOST_MUTATION_SCALE, MUTATION_PROP)
from polymodel.run import no_joblib_multiple_run


def main(run, host_on, N_ITS=5, N_YEARS=15):

    cmh = Config(
        type='multi',
        sprays=[2],
        host_on=[host_on],
        n_years=N_YEARS,
        n_k=100,
        n_l=100,
        n_iterations=N_ITS,
        mutation_proportion=MUTATION_PROP,
        mutation_scale_fung=FUNG_MUTATION_SCALE*DEFAULT_P,
        mutation_scale_host=HOST_MUTATION_SCALE*DEFAULT_P,
    )

    cmh.beta_multi = cmh.beta_multi[run*N_ITS*N_YEARS:(1+run)*N_ITS*N_YEARS]

    multi_runs = no_joblib_multiple_run(cmh)

    if host_on:
        host_str = 'Y'
    else:
        host_str = 'N'

    output = multi_runs[f'spray_Y2_host_{host_str}']

    df_out = get_dataframe(output, run, N_ITS)

    return df_out, cmh


def get_dataframe(output, run, n_its):

    all_runs = pd.DataFrame()

    for ii in range(n_its):
        index = n_its*run + ii
        dist = output[ii]['fung_dists']
        tv = output[ii]['k_vec']

        trait_mean = np.asarray([
            np.dot(
                dist[:, yr],
                tv
            )
            for yr in range(dist.shape[1])
        ]
        )

        this_run = pd.DataFrame(
            {
                f'yld{index}': np.concatenate([[np.nan], output[ii]['yield_vec']]),
                f'sev{index}': np.concatenate([[np.nan], 100*output[ii]['dis_sev']]),
                f'fung_mean{index}': trait_mean,
            }
        )

        all_runs = pd.concat([all_runs, this_run], axis=1)

    return all_runs


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a run index")

    run_index = int(sys.argv[1])

    result, config = main(run_index, True)

    result.to_csv(
        f'../outputs/fig5_T_{run_index}_{config.n_k}_{config.n_l}.csv'
    )
