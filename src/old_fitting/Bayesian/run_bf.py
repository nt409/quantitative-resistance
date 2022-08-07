import pickle

from posterior_fns import Likelihood, Priors
from plotting_fns import (
    plot_output, plot_sevs, contour_plot, histogram,
    likelihood_contour, prior_hist, plot_sampled_params
)

# location = 'Soenderborg'
locations = ['Soenderborg', 'Karise', 'Trige']
N_grid = 400

N_gen = 10000
sd_obs = 2.5
use_uninformative_prior = False


load_saved_lh = False

plot_likelihood_space = True
plot_hist = False
plot_severities = True
plot_samples = True


# for sd_obs in 2.5]:
for location in locations:

    my_lh = Likelihood(N_grid, location, sd_obs,
                       load_saved_lh, use_uninformative_prior)
    my_lh.whole_space()
    my_lh.generate_severities(N_gen)

    if plot_likelihood_space:

        fig = likelihood_contour(my_lh.likelihood, my_lh.edgeI0, my_lh.edgeB)
        fig.show()
        filename = f'../figures/Fitting/Bayes/Contour/' + \
            f'{location}_{N_grid}_{sd_obs}_{use_uninformative_prior}.png'
        fig.write_image(filename)

    if plot_hist:

        priors = Priors()

        for bool in [True, False]:
            fig = prior_hist(priors, bool)
            fig.show()
            if bool:
                string = 'I0'
            else:
                string = 'beta'

            filename = f'../figures/Fitting/Bayes/Prior/hist_{string}.png'
            fig.write_image(filename)

    if plot_severities:
        fig = plot_sevs(my_lh.data, my_lh.generated)
        fig.show()

        filename = f'../figures/Fitting/Bayes/Generated/_' + \
            f'{location}_{N_grid}_{sd_obs}_{use_uninformative_prior}_{N_gen}.png'
        fig.write_image(filename)

    if plot_samples:
        fig = plot_sampled_params(my_lh.logI0s, my_lh.betas)
        fig.show()

        filename = f'../figures/Fitting/Bayes/Sampled/' + \
            f'{location}_{N_grid}_{sd_obs}_{use_uninformative_prior}_{N_gen}.png'
        fig.write_image(filename)
