"""For Figure 6"""

import pandas as pd


def combine():

    res = pd.concat([
        # pd.read_csv(f'../outputs/theoretical_host_{ii}.csv')
        pd.read_csv(f'../outputs/fig6_theoretical_host_{ii}.csv')
        for ii in range(99)
    ], axis=0)

    print(res.shape)

    filename = '../outputs/combined/fig_6_theoretical.csv'
    print(f'saving to {filename}')

    res.to_csv(filename)

    return None


if __name__ == "__main__":
    combine()
