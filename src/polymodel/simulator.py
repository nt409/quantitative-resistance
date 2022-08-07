import numpy as np
from scipy.integrate import ode

from polymodel.params import PARAMS
from polymodel.utils import (
    Fungicide,
    economic_yield_function,
    get_dist_mean,
    get_dispersal_kernel,
    get_fung_dist_params_from_config,
    get_host_dist_params_from_config,
    trait_vec,
    host_growth_function,
    initial_host_dist,
    initial_fung_dist,
    normalise,
    yield_function
)


class SimulatorOneTrait:
    """Sets up and runs a single model simulation for fung OR host trait.
    """

    def __init__(
        self,
        config,
        fungicide_on=True,
        host_plant_on=False,
        number_of_sprays=0,
    ):
        """Init method

        Parameters
        ----------
        config : Config
            See docs for Config

        host_plant_on : bool, optional
            Whether host plant protection is on, by default True
            Used instead of config, since might be iterating over multiple values
            a list in the config.

        fungicide_on : bool, optional
            Whether fungicide protection is on, by default True.
            Used instead of config, since might be iterating over multiple values
            a list in the config.

        number_of_sprays : int, optional
            N sprays per year, by default 0
        """

        self.host_plant_on = host_plant_on
        self.fungicide_on = fungicide_on

        if host_plant_on:
            self.number_of_sprays = 0
        else:
            self.number_of_sprays = number_of_sprays

        self.config_o = config
        self.mutation_array = None

        self.k_a, self.k_b = get_fung_dist_params_from_config(self.config_o)
        self.l_a, self.l_b = get_host_dist_params_from_config(self.config_o)

        self.n_k = self.config_o.n_k
        self.n_l = self.config_o.n_l

        self.k_vec = trait_vec(self.n_k)
        self.l_vec = trait_vec(self.n_l)

        self.initial_k_dist = initial_fung_dist(self.n_k, self.k_a, self.k_b)
        self.initial_l_dist = initial_host_dist(self.n_l, self.l_a, self.l_b)

        #
        #

    def run_model(self, I0_vec, beta_vec):
        """_summary_

        Parameters
        ----------
        I0_vec : np.array/list
            A single I0 value per year
        beta_vec : np.array/list
            A single beta value per year

        Returns
        -------
        dict
            keys:
            - fung_dists: np.array, shape (n_k, n_years+1) - includes year 0
            - host_dists: np.array, shape (n_l, n_years+1)

            - n_k: int
            - n_l: int

            - I0_vec: np.array, shape (n_years, )
            - beta_vec: np.array, shape (n_years, )

            - t: np.array, shape (n_timepoints, )
            - y: np.array, shape (1+ [n_k OR n_l OR n_k*n_l], n_t, n_years)

            - total_I: np.array, shape (n_timepoints, n_years)

            - dis_sev: np.array, shape (n_years, )
            - yield_vec: np.array, shape (n_years, )
            - econ: np.array, shape (n_years, )
        """

        replace_cultivar_array = self.config_o.replace_cultivars

        fung_dists = np.zeros((self.n_k, len(beta_vec)+1))
        host_dists = np.zeros((self.n_l, len(beta_vec)+1))

        fung_dists[:, 0] = self.initial_k_dist
        host_dists[:, 0] = self.initial_l_dist

        dis_sev = np.zeros(len(beta_vec))

        sprays_vec_use = self.number_of_sprays*np.ones(len(beta_vec))

        self._get_mutation_kernels()
        #
        # run 0th year
        (
            sol, t, fung_dists[:, 1], host_dists[:, 1], total_infection
        ) = self.calculate_ode_soln(
            fung_dists[:, 0],
            host_dists[:, 0],
            I0_vec[0],
            beta_vec[0],
            sprays_vec_use[0],
        )

        # get shape of solution arrays from the output
        sol_list = np.zeros((sol.shape[0], sol.shape[1], len(beta_vec)))
        total_I = np.zeros((len(t), len(beta_vec)))

        # set the first year
        sol_list[:, :, 0] = sol
        total_I[:, 0] = total_infection
        dis_sev[0] = total_I[-1, 0]

        S_end = sol_list[-1, -1, 0]
        dis_sev[0] = dis_sev[0] / (dis_sev[0] + S_end)

        if replace_cultivar_array is not None and replace_cultivar_array[0]:
            host_dists[:, 1] = self.initial_l_dist

        #
        # calculate the rest of the years
        # for yr in tqdm(range(1, len(beta_vec))):
        for yr in range(1, len(beta_vec)):

            (
                sol_list[:, :, yr], t, fung_dists[:, yr+1],
                host_dists[:, yr+1], total_I[:, yr]
            ) = self.calculate_ode_soln(
                fung_dists[:, yr],
                host_dists[:, yr],
                I0_vec[yr],
                beta_vec[yr],
                sprays_vec_use[yr],
            )

            dis_sev[yr] = total_I[-1, yr]

            # scale so proportion of final leaf size
            S_end = sol_list[-1, -1, yr]
            dis_sev[yr] = dis_sev[yr] / (dis_sev[yr] + S_end)

            if (replace_cultivar_array is not None
                    and replace_cultivar_array[yr]):
                host_dists[:, yr+1] = self.initial_l_dist

        # calculate yield and economic yield
        yield_vec = [yield_function(sev) for sev in dis_sev]

        econ = [
            economic_yield_function([yield_], int(spray_num))[0]
            for yield_, spray_num in zip(yield_vec, sprays_vec_use)
        ]

        fung_means = get_dist_mean(fung_dists, self.k_vec)
        host_means = get_dist_mean(host_dists, self.l_vec)

        return {
            'fung_dists': fung_dists,
            'host_dists': host_dists,
            'fung_means': fung_means,
            'host_means': host_means,
            'n_k': self.n_k,
            'n_l': self.n_l,
            'k_vec': self.k_vec,
            'l_vec': self.l_vec,
            'I0_vec': I0_vec,
            'beta_vec': beta_vec,
            't': t,
            'y': sol_list,
            'total_I': total_I,
            'dis_sev': dis_sev,
            'yield_vec': np.asarray(yield_vec),
            'econ': np.asarray(econ)
        }

    def calculate_ode_soln(self, D0_k, D0_l, I0_in, beta_in, num_sprays):

        self._get_y0(I0_in, D0_l, D0_k)

        solution_tmp, solutiont = self._solve_it(beta_in, num_sprays)

        solution, fung_dist_out, host_dist_out = self._generate_new_dists(
            solution_tmp, D0_l, D0_k)

        total_infection = [sum(solution[:-1, tt])
                           for tt in range(len(solutiont))]
        total_infection = np.asarray(total_infection)

        return solution, solutiont, fung_dist_out, host_dist_out, total_infection

    def _ode_system_with_mutation(
            self,
            t,
            y,
            beta,
            host_growth_fn,
            strains_dict,
            fungicide,
            mutation_array
    ):

        dydt = np.zeros(len(self.y0))

        S = y[-1]

        I_in = y[:-1]

        # host effect is value of strain (host on) or np.ones (host off)
        host_effect_vec = np.array(strains_dict['host'])

        # fung effect is value of strain into fungicide effect function
        fung_effect_vec = np.array([
            fungicide.effect(strain, t) for strain in strains_dict['fung']
        ])

        offspring = host_effect_vec * fung_effect_vec * I_in

        I_after_mutation = np.matmul(mutation_array, offspring)

        dydt[:-1] = beta * S * I_after_mutation

        dydt[-1] = host_growth_fn(t, S, y) - beta * S * sum(I_after_mutation)

        return dydt

    def _get_y0(self, I0_in, D0_l, D0_k):
        """
        y0 different if we have two active traits vs one
        """

        S0_prop = 1 - I0_in

        # set initial condition
        if self.fungicide_on and not self.host_plant_on:
            y0_use = np.zeros(self.n_k+1)
            y0_use[:-1] = I0_in*D0_k

        elif self.host_plant_on and not self.fungicide_on:
            y0_use = np.zeros(self.n_l+1)
            y0_use[:-1] = I0_in*D0_l

        y0_use[-1] = S0_prop

        self.y0 = PARAMS.host_growth_initial_area*y0_use

    def _get_mutation_kernels(self):

        if self.fungicide_on and not self.host_plant_on:
            # FUNG
            self.mutation_array = get_dispersal_kernel(
                self.k_vec,
                self.config_o.mutation_proportion,
                self.config_o.mutation_scale_fung,
            )

        elif not self.fungicide_on and self.host_plant_on:
            # HOST
            self.mutation_array = get_dispersal_kernel(
                self.l_vec,
                self.config_o.mutation_proportion,
                self.config_o.mutation_scale_host,
            )

    def _solve_it(self, beta_in, num_sprays):

        ode_solver = ode(self._ode_system_with_mutation)

        ode_solver.set_integrator('dopri5', max_step=10)

        t_out = np.linspace(PARAMS.T_1, PARAMS.T_end, 100)
        t_out = list(t_out)

        # += PARAMS.T_1, no longer needed!
        t_out += [PARAMS.T_2, PARAMS.T_3]
        t_out = sorted(t_out)
        t_out = np.asarray(t_out)

        y_out = np.zeros((self.y0.shape[0], len(t_out)))

        ode_solver.set_initial_value(self.y0, t_out[0])

        strains_dict = {}
        if self.fungicide_on and not self.host_plant_on:
            # NB looking to match length of active fung trait vector, not host
            strains_dict['host'] = np.ones(self.n_k)
            strains_dict['fung'] = np.asarray(self.k_vec)

        elif not self.fungicide_on and self.host_plant_on:
            # NB looking to match length of active host trait vector, not fung
            strains_dict['host'] = np.asarray(self.l_vec)
            strains_dict['fung'] = np.ones(self.n_l)

        # add other params
        my_fungicide = Fungicide(num_sprays)

        # comes from utils
        host_growth_fn = host_growth_function

        ode_solver.set_f_params(
            beta_in,
            host_growth_fn,
            strains_dict,
            my_fungicide,
            self.mutation_array
        )

        for ind, tt in enumerate(t_out[1:]):
            if ode_solver.successful():
                y_out[:, ind] = ode_solver.y
                ode_solver.integrate(tt)
            else:
                raise RuntimeError('ode solver unsuccessful')

        y_out[:, -1] = ode_solver.y

        return y_out, t_out

    def _generate_new_dists(self, solution, D0_l, D0_k):

        I_end = normalise(solution[:-1, -1])

        n_t_points = solution.shape[1]

        if self.fungicide_on and not self.host_plant_on:
            soln_out = np.zeros((self.n_k+1, n_t_points))
            soln_out[:-1, :] = solution[:-1, :]

            fung_dist_out = np.zeros(self.n_k)
            fung_dist_out = I_end
            host_dist_out = D0_l

        elif self.host_plant_on and not self.fungicide_on:
            soln_out = np.zeros((self.n_l+1, n_t_points))
            soln_out[:-1, :] = solution[:-1, :]

            host_dist_out = np.zeros(self.n_l)
            host_dist_out = I_end
            fung_dist_out = D0_k

        # susceptible tissue
        soln_out[-1, :] = solution[-1, :]

        return soln_out, fung_dist_out, host_dist_out


#
#
#


class SimulatorBothTraits:
    """Sets up and runs a single model simulation in the fung AND host case
    """

    def __init__(
        self,
        config,
        number_of_sprays=0,
    ):
        """Init method

        Parameters
        ----------
        config : Config
            See docs for Config

        number_of_sprays : int, optional
            N sprays per year, by default 0
        """

        self.number_of_sprays = number_of_sprays

        self.config_b = config

        self.k_a, self.k_b = get_fung_dist_params_from_config(self.config_b)
        self.l_a, self.l_b = get_host_dist_params_from_config(self.config_b)

        self.n_k = self.config_b.n_k
        self.n_l = self.config_b.n_l

        self.k_vec = trait_vec(self.n_k)
        self.l_vec = trait_vec(self.n_l)

        self.fung_kernel = None
        self.host_kernel = None

        self.initial_k_dist = initial_fung_dist(self.n_k, self.k_a, self.k_b)
        self.initial_l_dist = initial_host_dist(self.n_l, self.l_a, self.l_b)

        #
        #

    def run_model(self, I0_vec, beta_vec):
        """_summary_

        Parameters
        ----------
        I0_vec : np.array/list
            A vector of I0 values
        beta_vec : np.array/list
            A vector of beta values

        Returns
        -------
        dict
            keys:
            - fung_dists: np.array, shape (n_k, n_years+1) - includes year 0
            - host_dists: np.array, shape (n_l, n_years+1)

            - n_k: int
            - n_l: int

            - I0_vec: np.array, shape (n_years, )
            - beta_vec: np.array, shape (n_years, )

            - t: np.array, shape (n_timepoints, )
            - y: np.array, shape (1+ [n_k OR n_l OR n_k*n_l], n_t, n_years)

            - total_I: np.array, shape (n_timepoints, n_years)

            - dis_sev: np.array, shape (n_years, )
            - yield_vec: np.array, shape (n_years, )
            - econ: np.array, shape (n_years, )
        """

        replace_cultivar_array = self.config_b.replace_cultivars

        fung_dists = np.zeros((self.n_k, len(beta_vec)+1))
        host_dists = np.zeros((self.n_l, len(beta_vec)+1))

        fung_dists[:, 0] = self.initial_k_dist
        host_dists[:, 0] = self.initial_l_dist

        dis_sev = np.zeros(len(beta_vec))

        sprays_vec_use = self.number_of_sprays*np.ones(len(beta_vec))

        self._get_kernels()
        #
        # run 0th year
        (
            sol, t, fung_dists[:, 1], host_dists[:, 1], total_infection
        ) = self.calculate_ode_soln(
            fung_dists[:, 0],
            host_dists[:, 0],
            I0_vec[0],
            beta_vec[0],
            sprays_vec_use[0],
        )

        # get shape of solution arrays from the output
        sol_list = np.zeros((sol.shape[0], sol.shape[1], len(beta_vec)))
        total_I = np.zeros((len(t), len(beta_vec)))

        # set the first year
        sol_list[:, :, 0] = sol
        total_I[:, 0] = total_infection
        dis_sev[0] = total_I[-1, 0]

        S_end = sol_list[-1, -1, 0]
        dis_sev[0] = dis_sev[0] / (dis_sev[0] + S_end)

        if replace_cultivar_array is not None and replace_cultivar_array[0]:
            host_dists[:, 1] = self.initial_l_dist

        #
        # calculate the rest of the years
        # for yr in tqdm(range(1, len(beta_vec))):
        for yr in range(1, len(beta_vec)):

            (
                sol_list[:, :, yr], t, fung_dists[:, yr+1],
                host_dists[:, yr+1], total_I[:, yr]
            ) = self.calculate_ode_soln(
                fung_dists[:, yr],
                host_dists[:, yr],
                I0_vec[yr],
                beta_vec[yr],
                sprays_vec_use[yr],
            )

            dis_sev[yr] = total_I[-1, yr]

            # scale so proportion of final leaf size
            S_end = sol_list[-1, -1, yr]
            dis_sev[yr] = dis_sev[yr] / (dis_sev[yr] + S_end)

            if (replace_cultivar_array is not None
                    and replace_cultivar_array[yr]):
                host_dists[:, yr+1] = self.initial_l_dist

        # calculate yield and economic yield
        yield_vec = [yield_function(sev) for sev in dis_sev]

        econ = [
            economic_yield_function([yield_], int(spray_num))[0]
            for yield_, spray_num in zip(yield_vec, sprays_vec_use)
        ]

        fung_means = get_dist_mean(fung_dists, self.k_vec)
        host_means = get_dist_mean(host_dists, self.l_vec)

        return {
            'fung_dists': fung_dists,
            'host_dists': host_dists,
            'fung_means': fung_means,
            'host_means': host_means,
            'n_k': self.n_k,
            'n_l': self.n_l,
            'k_vec': self.k_vec,
            'l_vec': self.l_vec,
            'I0_vec': I0_vec,
            'beta_vec': beta_vec,
            't': t,
            'y': sol_list,
            'total_I': total_I,
            'dis_sev': dis_sev,
            'yield_vec': np.asarray(yield_vec),
            'econ': np.asarray(econ)
        }

    def calculate_ode_soln(self, D0_k, D0_l, I0_in, beta_in, num_sprays):

        self._get_y0(I0_in, D0_l, D0_k)

        soln_tmp, t_out = self._solve_it(beta_in, num_sprays)

        soln, fung_dist_out, host_dist_out = self._dists_fung_and_host(
            soln_tmp)

        total_infection = [sum(soln[:-1, tt]) for tt in range(len(t_out))]

        total_infection = np.asarray(total_infection)

        return soln, t_out, fung_dist_out, host_dist_out, total_infection

    def _ode_system_host_and_fung(
            self,
            t,
            y,
            beta,
            host_growth_fn,
            fungicide,
            fung_kernel,
            host_kernel,
    ):
        # 1. offspring
        # 2. mutation

        dydt = np.zeros(len(self.y0))

        # host effect is same as value of strain
        host_effect_vec = np.asarray(self.l_vec)

        # FOR NO HOST EFFECT
        # host_effect_vec = np.ones(self.n_l)

        # fung effect is value of strain into fungicide effect function,
        fung_effect_vec = np.array([
            fungicide.effect(strain, t) for strain in self.k_vec
        ])

        S = y[-1]

        I_in = y[:-1]

        I_array = np.reshape(I_in, (self.n_k, self.n_l))

        I_fung = I_array.sum(axis=1)
        I_host = I_array.sum(axis=0)

        # scale both of these as a proportion
        p_fung = I_fung / I_fung.sum()
        p_host = I_host / I_host.sum()

        # 1. rate of offspring production depends on control of parent
        fung_offspring = fung_effect_vec * p_fung
        host_offspring = host_effect_vec * p_host

        # 2. child may not have same trait as parent
        fung_offsp_mut = np.matmul(fung_kernel, fung_offspring)
        host_offsp_mut = np.matmul(host_kernel, host_offspring)

        # convert back into vector with length n_k*n_l, scaled by sum(I_in)
        scale = sum(I_in)

        disease_states_array = np.outer(
            fung_offsp_mut,
            host_offsp_mut
        ) * scale

        disease_states = np.reshape(
            disease_states_array,
            (self.n_k*self.n_l)
        )

        dydt[:-1] = beta * S * disease_states
        dydt[-1] = host_growth_fn(t, S, y) - beta * S * sum(disease_states)

        return dydt

    def _get_y0(self, I0_in, D0_l, D0_k):
        """
        y0 different if we have two active traits vs one
        """

        S0_prop = 1 - I0_in

        # set initial condition
        y0_use = np.zeros(self.n_k*self.n_l+1)

        infection_array = I0_in * np.outer(D0_k, D0_l)

        infection_vector = np.reshape(
            infection_array,
            (self.n_k*self.n_l)
        )

        y0_use[:-1] = infection_vector

        y0_use[-1] = S0_prop

        self.y0 = PARAMS.host_growth_initial_area*y0_use

    def _get_kernels(self):

        self.fung_kernel = get_dispersal_kernel(
            self.k_vec,
            self.config_b.mutation_proportion,
            self.config_b.mutation_scale_fung,
        )

        self.host_kernel = get_dispersal_kernel(
            self.l_vec,
            self.config_b.mutation_proportion,
            self.config_b.mutation_scale_host,
        )

    def _solve_it(self, beta_in, num_sprays):

        ode_solver = ode(self._ode_system_host_and_fung)

        ode_solver.set_integrator('dopri5', max_step=10)

        t_out = np.linspace(PARAMS.T_1, PARAMS.T_end, 100)
        t_out = list(t_out)

        t_out += [PARAMS.T_2, PARAMS.T_3]
        t_out = sorted(t_out)
        t_out = np.asarray(t_out)

        y_out = np.zeros((self.y0.shape[0], len(t_out)))

        ode_solver.set_initial_value(self.y0, t_out[0])

        # add other params
        my_fungicide = Fungicide(num_sprays)

        # comes from utils
        host_growth_fn = host_growth_function

        ode_solver.set_f_params(
            beta_in,
            host_growth_fn,
            my_fungicide,
            self.fung_kernel,
            self.host_kernel,
        )

        for ind, tt in enumerate(t_out[1:]):
            if ode_solver.successful():
                y_out[:, ind] = ode_solver.y
                ode_solver.integrate(tt)
            else:
                raise RuntimeError('ode solver unsuccessful')

        y_out[:, -1] = ode_solver.y

        return y_out, t_out

    def _dists_fung_and_host(self, solution):
        I_end = solution[:-1, -1]

        n_t_points = solution.shape[1]

        I_end_array = np.reshape(I_end, (self.n_k, self.n_l))

        I0_k_end = I_end_array.sum(axis=1)

        I0_l_end = I_end_array.sum(axis=0)

        fung_dist_out = normalise(I0_k_end)
        host_dist_out = normalise(I0_l_end)

        soln_large_array = np.zeros(((self.n_k, self.n_l, n_t_points)))

        soln_small_array = np.reshape(
            solution[:-1, :],
            (I_end_array.shape[0], I_end_array.shape[1], n_t_points),
        )

        soln_large_array = soln_small_array

        solution_out = np.zeros((self.n_k*self.n_l+1, n_t_points))

        solution_out[:-1, :] = np.reshape(
            soln_large_array,
            (self.n_k * self.n_l, n_t_points),
        )

        # susceptible tissue
        solution_out[-1, :] = solution[-1, :]

        return solution_out, fung_dist_out, host_dist_out
