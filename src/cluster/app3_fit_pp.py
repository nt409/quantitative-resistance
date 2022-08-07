"""For Figure 3"""

import pandas as pd


def combine():

    p_values = [0.01, 0.1, 0.5, 0.9]

    comb = pd.concat([
        pd.read_csv(f"../outputs/fitting_{str(p).replace('.', ',')}.csv")
        for p in p_values
    ])

    print(comb.shape)

    f = '../outputs/combined/app3_fit.csv'
    print(f'saving to {f}')
    comb.to_csv(f, index=False)

    return None


if __name__ == "__main__":
    combine()
