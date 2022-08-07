"""For Figure 3"""

import sys

import numpy as np
import pandas as pd

from polymodel.config import Config
from polymodel.consts import (DEFAULT_P, FUNG_MUTATION_SCALE,
                              HOST_MUTATION_SCALE, MUTATION_PROP)
from polymodel.run import no_joblib_multiple_run


def main(
    run,
    N_ITS=5,
    N_YEARS=15
):

    cmh = Config(
        type='multi',
        sprays=[0, 1, 2, 3],
        host_on=[False],
        n_years=N_YEARS,
        n_k=500,
        n_l=50,
        n_iterations=N_ITS,
        mutation_proportion=MUTATION_PROP,
        mutation_scale_fung=FUNG_MUTATION_SCALE*DEFAULT_P,
        mutation_scale_host=HOST_MUTATION_SCALE*DEFAULT_P,
    )

    cmh.beta_multi = cmh.beta_multi[run*N_ITS*N_YEARS:(1+run)*N_ITS*N_YEARS]

    multi_runs = no_joblib_multiple_run(cmh)

    output0 = multi_runs['spray_N_host_N']
    output1 = multi_runs['spray_Y1_host_N']
    output2 = multi_runs['spray_Y2_host_N']
    output3 = multi_runs['spray_Y3_host_N']

    ylds_sevs_df_0, fung_dfs_0 = get_dataframe(output0, run, N_ITS)
    ylds_sevs_df_1, fung_dfs_1 = get_dataframe(output1, run, N_ITS)
    ylds_sevs_df_2, fung_dfs_2 = get_dataframe(output2, run, N_ITS)
    ylds_sevs_df_3, fung_dfs_3 = get_dataframe(output3, run, N_ITS)

    yld_sevs_out = pd.concat([
        ylds_sevs_df_0.assign(sprays=0),
        ylds_sevs_df_1.assign(sprays=1),
        ylds_sevs_df_2.assign(sprays=2),
        ylds_sevs_df_3.assign(sprays=3),
    ], axis=0)

    fung_dfs_out = pd.concat([
        fung_dfs_0.assign(sprays=0),
        fung_dfs_1.assign(sprays=1),
        fung_dfs_2.assign(sprays=2),
        fung_dfs_3.assign(sprays=3),
    ], axis=0)

    return yld_sevs_out, fung_dfs_out, cmh


def get_dataframe(output, run, n_its):
    """Get dataframe which has 

    Parameters
    ----------
    output : list of dicts
        --
    run : int
        --
    n_its : int
        --

    Returns
    -------
    all_ylds_sevs_df : pd.DataFrame
        all yields and sevs, columns:
        - yld
        - sev
        - run
        - year

        then can groupby year and take median

    all_fung_dfs : pd.DataFrame
        all fung dists, columns:
        - fung_dist_ii for ii in range(n_k)
        - run
        - year

        then can groupby year and take median
    """

    all_ylds_sevs_df = pd.DataFrame()
    all_fung_dfs = pd.DataFrame()

    for ii in range(n_its):
        index = n_its*run + ii
        dist = output[ii]['fung_dists']

        yld_sev_df = pd.DataFrame(
            {
                'yld': np.concatenate([[np.nan], output[ii]['yield_vec']]),
                'sev': np.concatenate([[np.nan], 100*output[ii]['dis_sev']]),
                'run': index,
                'year': np.arange(1+len(output[ii]['yield_vec']))
            }
        )

        fung_df = (
            pd.DataFrame(dist)
            .T
            .rename(columns=lambda x: 'fung_dist_' + str(x))
            .assign(
                run=index,
                year=np.arange(dist.shape[1]),
            )
        )

        all_fung_dfs = pd.concat([all_fung_dfs, fung_df], axis=0)

        all_ylds_sevs_df = pd.concat([all_ylds_sevs_df, yld_sev_df], axis=0)

    return all_ylds_sevs_df, all_fung_dfs


if __name__ == '__main__':

    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a run index")

    run_index = int(sys.argv[1])

    yld_and_sevs_out, fung_dists_out, config = main(run_index)

    conf_str = f'{run_index}_{config.n_k}_{config.n_l}'
    yld_and_sevs_out.to_csv(f'../outputs/fig3_ylds_sevs_df_{conf_str}.csv')
    fung_dists_out.to_csv(f'../outputs/fig3_fung_df_{conf_str}.csv')
