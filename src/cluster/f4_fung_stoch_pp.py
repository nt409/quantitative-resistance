"""For Figure 4"""

import numpy as np
import pandas as pd


def main():

    for name in [
        'ds_d02',
        'ds_d12',
        'ds_d32',
        'y_d02',
        'y_d12',
        'y_d32',
        'fd_diff02',
        'fd_diff12',
        'fd_diff32',
    ]:
        print(name)

        (
            pd.concat(
                [
                    pd.read_csv(f'../outputs/fig4_{name}_{int(run)}.csv')
                    for run in np.arange(100)
                ],
                axis=1)
            .to_csv(f'../outputs/combined/fig4_{name}.csv')
        )

    return None


if __name__ == '__main__':
    main()
