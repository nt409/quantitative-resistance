from math import log
import scipy.stats
import numpy as np
import pandas as pd

from posterior_fns import Priors


priors = Priors()


y_std = priors.logI0s


chi_square_statistics = []
# 11 equi-distant bins of observed Data

percentile_bins = np.linspace(0, 100, 11)
percentile_cutoffs = np.percentile(y_std, percentile_bins)

observed_frequency, bins = (np.histogram(y_std, bins=percentile_cutoffs))
# cum_observed_frequency = np.cumsum(observed_frequency)

dist_names = ["weibull_min", "norm", "weibull_max", "beta",
              "invgauss", "uniform", "gamma", "expon",
              "lognorm", "pearson3", "triang"]

# Loop through candidate distributions
for distribution in dist_names:
    # Set up distribution and get fitted distribution parameters
    dist = getattr(scipy.stats, distribution)
    param = dist.fit(y_std)
    print(f"{dist} {param}")

    # Get expected counts in percentile bins
    # cdf of fitted sistrinution across bins
    cdf_fitted = dist.cdf(percentile_cutoffs, *param)
    expected_frequency = []
    for bin in range(len(percentile_bins)-1):
        expected_cdf_area = cdf_fitted[bin+1] - cdf_fitted[bin]
        expected_frequency.append(expected_cdf_area)

    # Chi-square Statistics
    expected_frequency = np.array(expected_frequency) * len(y_std)
    # cum_expected_frequency = np.cumsum(expected_frequency)

    print(expected_frequency)
    print(observed_frequency)
    print(((expected_frequency - observed_frequency) ** 2) / expected_frequency)

    ss = sum(((expected_frequency - observed_frequency) ** 2) / expected_frequency)
    chi_square_statistics.append(ss)


# Sort by minimum ch-square statistics
results = pd.DataFrame()
results['Distribution'] = dist_names
results['chi_square'] = chi_square_statistics
results.sort_values(['chi_square'], inplace=True)


print('Distributions listed by Betterment of fit:')
print('............................................')
print(results)
