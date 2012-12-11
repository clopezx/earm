import earm_anneal
import pysb.integrate
import pysb.util
import numpy as np
import matplotlib.pyplot as plt

from earm.lopez_embedded import model

#pysb.util.update_param_vals(model, pysb.util.load_params('test_fit.txt'))
exp_data = np.load('xpdata/forfits/EC-RP_IMS-RP_IC-RP_data_for_models.npz')['arr_0']

exp_time = exp_data[0]
exp_obs_norm = exp_data[(2,5,8), :]

sim_data = pysb.integrate.odesolve(model, exp_time)
sim_obs = sim_data[model.observables.keys()].view(float).reshape(len(sim_data), -1)
sim_obs_norm = (sim_obs / sim_obs.max(0)).T

colors = ('r', 'g', 'b')
for exp, sim, c in zip(exp_obs_norm, sim_obs_norm, colors):
    plt.plot(exp_time, exp, color=c, marker='.', linestyle=':')
    plt.plot(exp_time, sim, color=c)
plt.show()
