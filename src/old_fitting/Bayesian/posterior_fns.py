from math import exp, log, floor
import numpy as np
import pandas as pd
import itertools
import scipy.stats
import pickle
import os
from tqdm import tqdm


np.random.seed(1)


def _uncontrolled_epidemic(I0, beta):
    T = 2500-1456
    out = 1/(((1-I0)/I0)*exp(-beta*T) + 1)
    return out


def generate_severities(I0s, betas):
    sevs = []
    for I0, beta in zip(I0s, betas):
        sevs.append(100*_uncontrolled_epidemic(I0, beta))
    return sevs


path = f'../../Data/I0/I0_Calculated.csv'
data_all = pd.read_csv(path)


class Ranges:
    def __init__(self):

        I0s = data_all.I0s
        betas = data_all.beta

        LI0_min = log(min(I0s)/10)
        LI0_max = log(max(I0s)*4)
        # LI0_max = log(0.01)
        # print(max(I0s)*4)

        beta_min = 10**(-13)
        beta_max = 1.25*max(betas)

        self.L_I0_range = [LI0_min, LI0_max]
        self.beta_range = [beta_min, beta_max]


class Priors:
    def __init__(self):

        I0s = data_all.I0s
        self.logI0s = [log(i) for i in I0s]
        self.betas = data_all.beta

        self.L_I_dist = scipy.stats.norm
        fitted_I = self.L_I_dist.fit(self.logI0s)
        self.L_I0_pars = fitted_I

        self.B_dist = scipy.stats.norm
        fitted_B = self.B_dist.fit(self.betas)
        self.Bpars = fitted_B


# likelihood

class Likelihood:
    def __init__(self, N, location, standard_deviation, load_saved, use_uninformative_prior):
        self.N = N
        self.location = location
        self.sd_obs = standard_deviation

        self.load_saved = load_saved

        self.uninformative_prior = use_uninformative_prior

        self.base_filename = f'Bayesian/pickled/{self.location}'

        self.data = self._get_data()

        ranges_ = Ranges()
        self.L_I0_range = ranges_.L_I0_range
        self.beta_range = ranges_.beta_range

        self.prior_dists = Priors()

        self.fitting_config_str = f"UninfPr={str(self.uninformative_prior)[0]}_N={self.N}_SD={self.sd_obs}"

    def _get_data(self):
        path = f'../../Data/Locations/{self.location}/WorstCults_TopProp6.csv'
        data_all = pd.read_csv(path)
        out = np.asarray(data_all.stb)

        out = out[out < 99]

        filename = self.base_filename + '_points.pickle'
        with open(filename, 'wb') as handle:
            pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return out

    def _single_point(self, L_I0, beta):

        I0 = exp(L_I0)

        F = 100*_uncontrolled_epidemic(I0, beta)

        # get likelihood
        likelihood_ = 0

        for dd in self.data:
            likelihood_ += exp(- ((dd - F)**2) / (2*(self.sd_obs**2)))

        # get prior
        if self.uninformative_prior:
            beta_prior = 1
            I0_prior = 1
        else:
            beta_prior = self.prior_dists.B_dist.pdf(beta,
                                                     *self.prior_dists.Bpars)

            I0_prior = self.prior_dists.L_I_dist.pdf(I0,
                                                     *self.prior_dists.L_I0_pars)

        prior = beta_prior * I0_prior

        posterior = likelihood_*prior
        return posterior

    def _setup_grid(self, range_in):

        jump = 0.5*(range_in[1] - range_in[0])/self.N

        out = np.linspace(range_in[0]+jump,
                          range_in[1]-jump,
                          self.N)
        return out

    def whole_space(self):

        self.edgeI0 = np.linspace(self.L_I0_range[0],
                                  self.L_I0_range[1],
                                  self.N+1)

        self.edgeB = np.linspace(self.beta_range[0],
                                 self.beta_range[1],
                                 self.N+1)

        name_use = f"{self.base_filename}_{self.fitting_config_str}_likelihood"
        lh_filename = name_use.replace(".", ",") + ".pickle"

        if self.load_saved:
            if os.path.isfile(lh_filename):
                lh_loaded = pickle.load(open(lh_filename, 'rb'))
                z = lh_loaded['z']
                grid_out = lh_loaded['grid']

                self.likelihood = z

                return None

        # else run and save
        # middleI0, beta
        I0_vals = self._setup_grid(self.L_I0_range)
        beta_vals = self._setup_grid(self.beta_range)

        z = np.zeros((len(beta_vals), len(I0_vals)))

        for bb, ii in tqdm(itertools.product(range(len(beta_vals)), range(len(I0_vals)))):
            beta = beta_vals[bb]
            I0 = I0_vals[ii]

            z[bb, ii] = self._single_point(I0, beta)

        self.likelihood = z

        grid_out = dict(
            edgeI0=self.edgeI0,
            edgeB=self.edgeB
        )

        # save
        lh_save = dict(grid=grid_out, z=z)
        with open(lh_filename, 'wb') as handle:
            pickle.dump(lh_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return None

    def generate_severities(self, N_gen):

        object_names = ['generated', 'logI0s', 'betas']

        base_name = f"{self.base_filename}_Pr={str(self.uninformative_prior)[0]}_N={self.N}_SD={self.sd_obs}_Gen={N_gen}_"
        filenames = {}
        for name in object_names:
            name_use = base_name + name
            filenames[name] = name_use.replace('.', ',') + '.pickle'

        if self.load_saved:
            loaded_all = True

            gen_out = {}
            for name in object_names:

                filename = filenames[name]

                if os.path.isfile(filename):
                    gen_out[name] = pickle.load(open(filename, 'rb'))
                else:
                    loaded_all = False

            if loaded_all:
                self.generated = gen_out['generated']
                self.logI0s = gen_out['logI0s']
                self.betas = gen_out['betas']

                return None

        # else run new and save

        gen_sevs = np.zeros(N_gen)
        logI0s = np.zeros(N_gen)
        betas = np.zeros(N_gen)

        lh = np.asarray(self.likelihood)
        norm_lh = lh/np.sum(lh)

        oneD = norm_lh.ravel()

        twoDindices = np.arange(len(oneD)).reshape(norm_lh.shape)
        ravelled = twoDindices.ravel()

        random_sample = np.random.choice(ravelled, len(gen_sevs), p=oneD)

        for i in range(len(gen_sevs)):
            indices = np.where(twoDindices == random_sample[i])

            sampled = {}

            for vec, index, name in zip([self.edgeB, self.edgeI0],
                                        indices,
                                        ['beta', 'logI0']):
                ind = int(index)
                lower_bd = vec[ind]
                upper_bd = vec[ind+1]

                sampled[name] = np.random.uniform(low=lower_bd, high=upper_bd)
                # sampled[name] = np.mean([lower_bd,upper_bd])

            betas[i] = sampled['beta']
            logI0s[i] = sampled['logI0']

            beta_in = betas[i]
            I0_in = exp(logI0s[i])

            G = 100*_uncontrolled_epidemic(I0_in, beta_in)

            gen_sevs[i] = G

        for obj, name in zip([gen_sevs, logI0s, betas],
                             object_names):

            filename = filenames[name]
            with open(filename, 'wb') as handle:
                pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

        self.generated = gen_sevs
        self.logI0s = logI0s
        self.betas = betas

        I0s = [exp(ii) for ii in logI0s]

        df_I0_out = pd.DataFrame(dict(I0s=I0s))
        df_I0_out.to_csv(
            f"../../Data/Locations/{self.location}/GeneratedI0s_{self.fitting_config_str}_Gen={N_gen}.csv")

        df_betas_out = pd.DataFrame(dict(betas=betas))
        df_betas_out.to_csv(
            f"../../Data/Locations/{self.location}/GeneratedBetas_{self.fitting_config_str}_Gen={N_gen}.csv")

        return None
