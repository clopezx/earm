import earm_anneal
import pysb.integrate
import pysb.util
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import os

from earm.lopez_embedded import model


obs_names = ('mBid', 'aSmac', 'cPARP')
data_names = ('norm_ICRP', 'IMSRP_step', 'norm_ECRP')

def objective_func(x, solver, rate_mask, data):
    param_values = np.array([p.value for p in model.parameters])
    param_values[rate_mask] = 10 ** x
    solver.run(param_values)

    #import ipdb; ipdb.set_trace()
    of_value = 0
    for obs_name, data_name in zip(obs_names, data_names):
        ysim = solver.yobs[obs_name]
        traj_min = np.nanmin(ysim)
        traj_max = np.nanmax(ysim)
        ysim_norm = (ysim - traj_min) / (traj_max - traj_min)
        ydata = data[data_name]
        #sigma = np.maximum(data[data_name] * 0.25, 0.01)  # .25 CV, min of 0.01
        sigma = 0.1
        of_value += np.sum((ydata - ysim_norm) ** 2 / (2 * sigma ** 2))

    print of_value, x[0:5]
    return of_value


if __name__ == '__main__':

    params = model.parameters
    rate_params = model.parameters_rules()
    rate_mask = np.array([p in rate_params for p in params])
    x0 = np.log10([p.value for p in rate_params])
    dx = 0.02
    lower = x0 - dx / 2
    upper = x0 + dx / 2

    earm_path = os.path.join(os.path.dirname(__file__))
    data_path = os.path.join(earm_path, 'xpdata', 'forfits',
                             'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
    data = np.genfromtxt(data_path, delimiter=',', names=True)

    solver = pysb.integrate.Solver(model, data['Time'], rtol=1e-3, atol=1e-3)

    np.random.seed(1)
    (xmin, Jmin, Tfinal, feval, iters, accept, retval) = \
        scipy.optimize.anneal(objective_func, x0, full_output=True,
                              lower=lower, upper=upper, disp=True,
                              args=[solver, rate_mask, data])
    params_estimated = np.array([p.value for p in params])
    params_estimated[rate_mask] = 10 ** xmin


    exp_data = np.load('xpdata/forfits/EC-RP_IMS-RP_IC-RP_data_for_models.npz')['arr_0']
    exp_time = exp_data[0]
    exp_obs_norm = exp_data[(2,5,8), :]

    solver.run(params_estimated)
    sim_obs_norm = (solver.yobs_view / solver.yobs_view.max(0)).T

    colors = ('r', 'g', 'b')
    for exp, sim, c in zip(exp_obs_norm, sim_obs_norm, colors):
        plt.plot(exp_time, exp, color=c, marker='.', linestyle=':')
        plt.plot(exp_time, sim, color=c)
    plt.show()
