"""For Figure 6"""

import sys
import numpy as np

import pandas as pd
from polymodel.config import Config
from polymodel.consts import DEFAULT_BETA, DEFAULT_I0, DEFAULT_P, FUNG_MUTATION_SCALE, HOST_MUTATION_SCALE, MUTATION_PROP
from polymodel.run import no_joblib_simulations_run

from polymodel.utils import get_dist_mean


def get_data_this_host_mean(
    means,  # host_means_and_vars,
    host_index
):

    N_K = 100
    N_L = 100
    N_YEARS = 50

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
    )

    # cf.l_b = host_means_and_vars.iloc[host_index].b
    # cf.l_mu = host_means_and_vars.iloc[host_index].mu
    cf.l_mu = means[host_index]

    res = no_joblib_simulations_run(
        cf, [DEFAULT_I0]*N_YEARS, [DEFAULT_BETA]*N_YEARS)

    out = pd.concat([
        get_df_this_host_mean(res, 1, cf),
        get_df_this_host_mean(res, 2, cf),
        get_df_this_host_mean(res, 3, cf),
    ])

    return out


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
    # df_in = pd.read_csv(
    # '../data/03_model_inputs/theoretical_host_means_and_variances.csv')

    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a run index")

    host_ind = int(sys.argv[1])

    means = np.linspace(0.01, 0.99, 99)

    result = get_data_this_host_mean(
        means,  # df_in,
        host_ind)

    result.to_csv(f'../outputs/fig6_theoretical_host_{host_ind}.csv')
    # result.to_csv(f'../outputs/theoretical_host_{host_ind}.csv')
