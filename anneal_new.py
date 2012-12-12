import pysb.integrate
import pysb.util
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import os

from earm.lopez_embedded import model


def objective_func(x, solver, rate_mask, exp_data, obs_names, data_names):
    # Simulate model with rates taken from x (which is log transformed)
    param_values = np.array([p.value for p in model.parameters])
    param_values[rate_mask] = 10 ** x
    solver.run(param_values)

    error = 0
    for obs_name, data_name in zip(obs_names, data_names):
        # Get model observable trajectory
        ysim = solver.yobs[obs_name]
        # Normalize it to 0-1
        traj_min = np.nanmin(ysim)
        traj_max = np.nanmax(ysim)
        if traj_max - traj_min > 0:
            ysim_norm = (ysim - traj_min) / (traj_max - traj_min)
        else:
            ysim_norm = 0
        # Get experimental measurement
        ydata = exp_data[data_name]
        # Allow for experimental error
        sigma = 0.1
        # Compute error between simulation and experiment (chi-squared)
        error += np.sum((ydata - ysim_norm) ** 2 / (2 * sigma ** 2))

    print 'error =', error
    return error


if __name__ == '__main__':

    print 'Estimating rates for model:', model.name

    np.random.seed(6)

    # Initial parameter vector is sampled uniformly from a hypercube in log
    # parameter space, centered around the nominal parameter values.
    # ==========
    # Get parameters for rates only
    rate_params = model.parameters_rules()
    # Build a boolean mask for those params against the entire param list
    rate_mask = np.array([p in rate_params for p in model.parameters])
    # Log-transformed parameter values
    xnominal = np.log10([p.value for p in rate_params])
    # Radius of hypercube for sampling
    rand_scale = 2
    # Uniformly sample a vector from the hypercube
    rand_offset = 2 * rand_scale * (np.random.rand(len(rate_params)) - 0.5)
    # Apply offset to nominal values
    x0 = xnominal + rand_offset
    # Displacement size for annealing moves
    dx = 0.02
    # The scipy documentation for 'lower' and 'upper' seems wrong. See
    # http://projects.scipy.org/scipy/ticket/1126 for more information. This is
    # how to get your search to start at x0 and use a displacement of dx:
    lower = x0 - dx / 2
    upper = x0 + dx / 2

    # Load experimental data file
    earm_path = os.path.dirname(__file__)
    data_path = os.path.join(earm_path, 'xpdata', 'forfits',
                             'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
    exp_data = np.genfromtxt(data_path, delimiter=',', names=True)

    # List of model observables and corresponding data file columns
    obs_names = ['mBid', 'aSmac', 'cPARP']
    data_names = ['norm_ICRP', 'IMSRP_step', 'norm_ECRP']

    # Initialize solver object
    solver = pysb.integrate.Solver(model, exp_data['Time'], rtol=1e-3, atol=1e-3)
    # Perform the annealing
    args = [solver, rate_mask, exp_data, obs_names, data_names]
    (xmin, Jmin, Tfinal, feval, iters, accept, retval) = \
        scipy.optimize.anneal(objective_func, x0, full_output=True,
                              lower=lower, upper=upper, feps=1e-2, args=args)
    # Construct vector with resulting parameter values (un-log-transformed)
    params_estimated = np.array([p.value for p in model.parameters])
    params_estimated[rate_mask] = 10 ** xmin

    # Write parameter values to a file
    fit_filename = 'fit_%s.txt' % model.name.replace('.', '_')
    fit_filename = os.path.join(earm_path, fit_filename)
    print 'Saving parameter values to file:', fit_filename
    pysb.util.write_params(model, params_estimated, fit_filename)

    # Construct 3-by-N matrix of experimental data columns of interest
    exp_obs_norm = exp_data[data_names].view(float).reshape(len(exp_data), -1).T

    # Simulate model with new parameters and construct a 3-by-N matrix of the
    # trajectories of the observables of interest, normalized to 0-1.
    solver.run(params_estimated)
    sim_obs = solver.yobs[obs_names].view(float).reshape(len(solver.yobs), -1)
    sim_obs_norm = (sim_obs / sim_obs.max(0)).T

    # Plot experimental data and simulation on the same axes
    colors = ('r', 'g', 'b')
    for exp, sim, c in zip(exp_obs_norm, sim_obs_norm, colors):
        plt.plot(exp_data['Time'], exp, color=c, marker='.', linestyle=':')
        plt.plot(exp_data['Time'], sim, color=c)
    plt.show()
