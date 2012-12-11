"""
Perform parameter estimation on rate constants using simulated annealing.
"""

import sys
import os
import re
import scipy.optimize
import numpy
import pysb.util
import pysb.integrate


def run(model):
    params = model.parameters
    rate_params = model.parameters_rules()
    rate_mask = numpy.array([p in rate_params for p in params])
    x0 = numpy.log10(numpy.array([p.value for p in rate_params]))

    earm_path = os.path.join(os.path.dirname(__file__), os.path.pardir)
    data_path = os.path.join(earm_path, 'xpdata', 'forfits',
                             'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
    data = numpy.genfromtxt(data_path, delimiter=',', names=True)

    (xmin, Jmin, Tfinal, feval, iters, accept, retval) = \
        scipy.optimize.anneal(objective_func, x0, full_output=True,
                              args=[model, rate_mask, data])
    params_estimated = numpy.array([p.value for p in params])
    params_estimated[rate_mask] = 10 ** xmin
    
    return pysb.util.write_params(model, params_estimated)


def objective_func(x, model, rate_mask, data):
    param_values = numpy.array([p.value for p in model.parameters])
    param_values[rate_mask] = 10 ** x
    y = pysb.integrate.odesolve(model, data['Time'], param_values,
                                rtol=1e-3, atol=1e-3)
    obs_names = ('mBid', 'aSmac', 'cPARP')
    data_names = ('norm_ICRP', 'IMSRP_step', 'norm_ECRP')

    of_value = 0
    for obs_name, data_name in zip(obs_names, data_names):
        ysim = y[obs_name]
        ysim_norm = (ysim - numpy.nanmin(ysim)) / (numpy.nanmax(ysim) - numpy.nanmin(ysim))
        ydata = data[data_name]
        #sigma = numpy.maximum(data[data_name] * 0.25, 0.01)  # .25 CV, min of 0.01
        sigma = 0.1
        #import pylab
        #pylab.ion()
        #pylab.figure()
        #pylab.plot(data['Time'], numpy.c_[ydata, ysim_norm, (ydata - ysim_norm) ** 2 / (2 * sigma ** 2)])
        of_value += numpy.sum((ydata - ysim_norm) ** 2 / (2 * sigma ** 2))

    return of_value


if __name__ == '__main__':
    # sanity checks on filename
    if len(sys.argv) <= 1:
        raise Exception("You must specify the filename of a model script")
    model_filename = sys.argv[1]
    if not os.path.exists(model_filename):
        raise Exception("File '%s' doesn't exist" % model_filename)
    if not re.search(r'\.py$', model_filename):
        raise Exception("File '%s' is not a .py file" % model_filename)
    sys.path.insert(0, os.path.dirname(model_filename))
    model_name = re.sub(r'\.py$', '', os.path.basename(model_filename))
    # import it
    try:
        # FIXME if the model has the same name as some other "real" module which we use,
        # there will be trouble (use the imp package and import as some safe name?)
        model_module = __import__(model_name)
    except StandardError as e:
        print "Error in model script:\n"
        raise
    # grab the 'model' variable from the module
    try:
        model = model_module.__dict__['model']
    except KeyError:
        raise Exception("File '%s' isn't a model file" % model_filename)
    print run(model)
