import scipy.optimize.anneal
import numpy
from earm_1_5biomods import model
import pysb.anneal_sundials
envlist, paramarr = pysb.anneal_sundials.annlinit(model)
xpfile = numpy.load('xp_mod_data_earm10.npz')
xpdata = xpfile['arr_0']

prmfile = numpy.load('earm_1_5annlst.npz')
annprm = prmfile['arr_0']
annnum = prmfile['arr_1']

lb, ub, lower, upper = pysb.anneal_sundials.getgenparambounds(annprm, omag=2, N=10)

annlout = scipy.optimize.anneal(pysb.anneal_sundials.annealfxn, annprm, args=(annnum, 16000, model, envlist, xpdata, 2, 2, lb, ub), lower=lower, upper=upper)
