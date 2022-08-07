import pandas as pd
import numpy as np

FOLDER_PATH = '../data/03_model_inputs'

DEFAULT_I0 = float(
    pd.read_csv(f'{FOLDER_PATH}/I0_value.csv')
    .loc[:, ['I0_value']]
    .iloc[0]
)


LAMBDA_FITTED = float(
    pd.read_csv(f'{FOLDER_PATH}/lambda_fitted.csv')
    .loc[:, ['lambda_fitted']]
    .iloc[0]
)

DEFAULT_BETA = float(
    pd.read_csv(f'{FOLDER_PATH}/beta_value.csv')
    .loc[:, ['beta_median']]
    .iloc[0]
)

ALL_BETAS = np.asarray(
    pd.read_csv(f'{FOLDER_PATH}/beta_sampled.csv').beta
)

FUNG_MUTATION_SCALE = float(
    pd.read_csv(f'{FOLDER_PATH}/mutation_scale.csv')
    .loc[:, ['mutation_scale']]
    .iloc[0]
)

HOST_MUTATION_SCALE = FUNG_MUTATION_SCALE

# From Alexey paper:
MUTATION_PROP = (0.5 * (28 + 130) * 1e6) / (0.5 * (2.3 + 10.5) * 1e12)

DEFAULT_P = 0.1

DEFAULT_MUT_SCALE_HOST = DEFAULT_P * HOST_MUTATION_SCALE
DEFAULT_MUT_SCALE_FUNG = DEFAULT_P * FUNG_MUTATION_SCALE

FUNG_DECAY_RATE = 0.5 * (6.91e-3 + 1.11e-2)
# OLD OPTIONS:
# FUNG_DECAY_RATE = 6.91e-3
# FUNG_DECAY_RATE = 1e-2
# FUNG_DECAY_RATE = 1.11e-2


TRAIN_TEST_SPLIT_PROPORTION = 2/3
# OLD OPTIONS:
# TRAIN_TEST_SPLIT_PROPORTION = 0.75
# TRAIN_TEST_SPLIT_PROPORTION = 0.8
# TRAIN_TEST_SPLIT_PROPORTION = 1
