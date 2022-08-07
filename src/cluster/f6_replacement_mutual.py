"""For Figure 6 Col 2"""

import sys

import numpy as np
import pandas as pd
from polymodel.config import Config
from polymodel.consts import DEFAULT_BETA, DEFAULT_I0, DEFAULT_P, FUNG_MUTATION_SCALE, HOST_MUTATION_SCALE, MUTATION_PROP
from polymodel.run import no_joblib_simulations_run

from polymodel.utils import get_dist_mean


def get_data_replacing_N_times(rep_index):
    """Get data for replacing cultivar every rep_index number of years

    Parameters
    ----------
    rep_index : int
        Replace after rep_index years, e.g. for 2 replace at end of 
        year 1, year 3, year 5 etc (starting from year 0).

    Returns
    -------
    out : pd.DataFrame
        _description_
    """

    N_K = 100
    N_L = 100
    N_YEARS = 100

    cf = Config(
        'single',
        n_k=N_K,
        n_l=N_L,
        sprays=[1, 2, 3],
        host_on=[True],
        n_years=N_YEARS,
        mutation_proportion=MUTATION_PROP,
        mutation_scale_fung=FUNG_MUTATION_SCALE*DEFAULT_P,
        mutation_scale_host=HOST_MUTATION_SCALE*DEFAULT_P,
        replace_cultivars=get_replacement_array(rep_index, N_YEARS)
    )

    res = no_joblib_simulations_run(
        cf, [DEFAULT_I0]*N_YEARS, [DEFAULT_BETA]*N_YEARS)

    out = pd.concat([
        get_df_this_host_mean(res, 1, cf),
        get_df_this_host_mean(res, 2, cf),
        get_df_this_host_mean(res, 3, cf),
    ])

    return out


def get_replacement_array(replacement_interval, n_years):
    replacement = np.zeros(n_years)

    for ii in range(replacement_interval - 1,
                    n_years,
                    replacement_interval):
        replacement[ii] = 1
    return replacement


def get_df_this_host_mean(sf3, sprays, cf):
    output = sf3[f'spray_Y{sprays}_host_Y']

    fung_mean = get_dist_mean(
        output['fung_dists'],
        output['k_vec']
    )

    host_mean = get_dist_mean(
        output['host_dists'],
        output['l_vec']
    )

    df = (
        pd.concat([
            pd.DataFrame([dict(host_mean_start=cf.l_mu)]),

            pd.DataFrame([dict(sprays=sprays)]),

            (
                pd.DataFrame(output['yield_vec'])
                .T
                .rename(columns=lambda x: 'yld' + str(x+1))
            ),

            (
                pd.DataFrame(fung_mean)
                .T
                .rename(columns=lambda x: 'fung_mean' + str(x))
            ),

            (
                pd.DataFrame(host_mean)
                .T
                .rename(columns=lambda x: 'host_mean' + str(x))
            ),
        ], axis=1)
    )

    return df


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a run index")

    rep_ind = int(sys.argv[1])

    result = get_data_replacing_N_times(rep_ind)

    result.to_csv(f'../outputs/fig6_replacing_host_{rep_ind}.csv')
