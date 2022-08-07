"""For Figure 6"""

import pandas as pd


def combine():

    res = pd.concat([
        pd.read_csv(f'../outputs/fig6_replacing_host_{ii}.csv')
        for ii in range(1, 50)
    ], axis=0)

    print(res.shape)

    filename = '../outputs/combined/fig_6_replacing.csv'

    print(f'saving to {filename}')

    res.to_csv(filename)

    return None


if __name__ == "__main__":
    combine()
