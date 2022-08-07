"""For Figure 4"""

import sys

from polymodel.config import Config
from polymodel.consts import (DEFAULT_P, FUNG_MUTATION_SCALE,
                              HOST_MUTATION_SCALE, MUTATION_PROP)
from polymodel.run import no_joblib_multiple_run

from plots.fns import arrange_as_data_frame, dist_means_as_df


def main(run, N_ITS=10, N_YEARS=15):

    cm = Config(
        type='multi',
        sprays=[0, 1, 2, 3],
        host_on=[False],
        n_years=N_YEARS,
        n_k=400,
        n_iterations=N_ITS,
        mutation_proportion=MUTATION_PROP,
        mutation_scale_fung=FUNG_MUTATION_SCALE*DEFAULT_P,
        mutation_scale_host=HOST_MUTATION_SCALE*DEFAULT_P,
    )

    cm.beta_multi = cm.beta_multi[run*N_ITS*N_YEARS:(1+run)*N_ITS*N_YEARS]

    multi_runs = no_joblib_multiple_run(cm)

    sol_0 = multi_runs[f'spray_N_host_N']
    sol_1 = multi_runs[f'spray_Y1_host_N']
    sol_2 = multi_runs[f'spray_Y2_host_N']
    sol_3 = multi_runs[f'spray_Y3_host_N']

    # SEV

    df_ds0 = arrange_as_data_frame(sol_0, 'dis_sev')
    df_ds1 = arrange_as_data_frame(sol_1, 'dis_sev')
    df_ds2 = arrange_as_data_frame(sol_2, 'dis_sev')
    df_ds3 = arrange_as_data_frame(sol_3, 'dis_sev')

    ds_d02 = 100*df_ds0 - 100*df_ds2
    ds_d12 = 100*df_ds1 - 100*df_ds2
    ds_d32 = 100*df_ds3 - 100*df_ds2

    ds_d02.to_csv(f'../outputs/fig4_ds_d02_{run}.csv')
    ds_d12.to_csv(f'../outputs/fig4_ds_d12_{run}.csv')
    ds_d32.to_csv(f'../outputs/fig4_ds_d32_{run}.csv')

    # YIELD

    df_y0 = arrange_as_data_frame(sol_0, 'yield_vec')
    df_y1 = arrange_as_data_frame(sol_1, 'yield_vec')
    df_y2 = arrange_as_data_frame(sol_2, 'yield_vec')
    df_y3 = arrange_as_data_frame(sol_3, 'yield_vec')

    y_d02 = df_y0 - df_y2
    y_d12 = df_y1 - df_y2
    y_d32 = df_y3 - df_y2

    y_d02.to_csv(f'../outputs/fig4_y_d02_{run}.csv')
    y_d12.to_csv(f'../outputs/fig4_y_d12_{run}.csv')
    y_d32.to_csv(f'../outputs/fig4_y_d32_{run}.csv')

    # DIST MEAN

    fdm0 = dist_means_as_df(sol_0, 'fung', cm)
    fdm1 = dist_means_as_df(sol_1, 'fung', cm)
    fdm2 = dist_means_as_df(sol_2, 'fung', cm)
    fdm3 = dist_means_as_df(sol_3, 'fung', cm)

    fd_diff02 = fdm0 - fdm2
    fd_diff12 = fdm1 - fdm2
    fd_diff32 = fdm3 - fdm2

    fd_diff02.to_csv(f'../outputs/fig4_fd_diff02_{run}.csv')
    fd_diff12.to_csv(f'../outputs/fig4_fd_diff12_{run}.csv')
    fd_diff32.to_csv(f'../outputs/fig4_fd_diff32_{run}.csv')

    return None


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a run index")

    run_index = int(sys.argv[1])

    main(run_index)
