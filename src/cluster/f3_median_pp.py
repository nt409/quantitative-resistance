"""For Figure 3"""

import pandas as pd

N_ITS = 2000
N_K = 500
N_L = 50


def combine_ys_df():

    ys_df = pd.concat([
        pd.read_csv(f'../outputs/fig3_ylds_sevs_df_{ii}_{N_K}_{N_L}.csv')
        for ii in range(N_ITS)
    ], axis=0)

    print(ys_df.shape)

    ys_gpd = (
        ys_df
        .groupby(['sprays', 'year'])
        .median()
        .reset_index()
        .drop('run', axis=1)
    )

    print(ys_gpd.shape)

    fn = '../outputs/combined/fig3_ys_df.csv'
    print(f'saving to {fn}')
    ys_gpd.to_csv(fn, index=False)

    return None


def combine_fung_df():

    f_df = pd.concat([
        pd.read_csv(f'../outputs/fig3_fung_df_{ii}_{N_K}_{N_L}.csv')
        for ii in range(N_ITS)
    ], axis=0)

    print(f_df.shape)

    f_gpd = (
        f_df
        .groupby(['sprays', 'year'])
        .mean()
        .reset_index()
        .drop('run', axis=1)
    )

    print(f_gpd.shape)

    fn = '../outputs/combined/fig3_f_df.csv'
    print(f'saving to {fn}')

    f_gpd.to_csv(fn, index=False)

    return None


if __name__ == "__main__":
    combine_ys_df()
    combine_fung_df()
