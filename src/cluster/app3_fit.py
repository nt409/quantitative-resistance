import sys

import optuna
import pandas as pd
from optuna.samplers import TPESampler
from polymodel.config import Config
from polymodel.consts import (FUNG_MUTATION_SCALE, HOST_MUTATION_SCALE,
                              MUTATION_PROP)
from polymodel.fitting import FungicideObjective


def fitting_fung_with_mutation(p, n_trials=600):

    config = Config(
        'single',
        n_k=500,
        n_l=10,
        mutation_proportion=MUTATION_PROP,
        mutation_scale_fung=p * FUNG_MUTATION_SCALE,
        mutation_scale_host=p * HOST_MUTATION_SCALE,
    )

    optuna.logging.set_verbosity(0)

    sampler = TPESampler(seed=0)
    study = optuna.create_study(sampler=sampler)
    obj = FungicideObjective(config)

    study.optimize(obj, n_trials=n_trials)

    out = (
        pd.DataFrame([study.best_params])
        .assign(
            p=p,
            loss=study.best_value,
        )
    )

    p_str = str(p).replace('.', ',')

    out.to_csv(f'../outputs/fitting_{p_str}.csv', index=False)


if __name__ == '__main__':

    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a run index")

    p_values = [0.01, 0.1, 0.5, 0.9]

    P_VAL = p_values[int(sys.argv[1])]

    fitting_fung_with_mutation(P_VAL)
