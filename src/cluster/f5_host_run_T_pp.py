"""For Figure 5"""

import pandas as pd


def combine(host_on_str, n_k, n_l):
    # host_on_str = T or F

    res = pd.concat([
        pd.read_csv(f'../outputs/fig5_{host_on_str}_{ii}_{n_k}_{n_l}.csv')
        for ii in range(200)
    ], axis=1)

    print(res.shape)

    fn = f'../outputs/combined/fig5_{host_on_str}.csv'

    print(f'saving to {fn}')

    res.to_csv(fn)

    return None


if __name__ == "__main__":
    combine('T', 100, 100)
