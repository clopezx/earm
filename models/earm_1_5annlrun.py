import scipy.optimize.anneal
import numpy
from earm_1_5biomods import model
import pysb.anneal_sundials
envlist, paramarr = pysb.anneal_sundials.annlinit(model)
xpfile = numpy.load('xp_mod_data_earm10.npz')
xpdata = xpfile['arr_0']

prmfile = numpy.load('earm_1_5smacannlst.npz')
smacprm = prmfile['arr_0']
smacnum = prmfile['arr_1']
smaclst = prmfile['arr_2']

lb, ub, lower, upper = pysb.anneal_sundials.getgenparambounds(smacprm, omag=3, N=1000)

annlout = scipy.optimize.anneal(pysb.anneal_sundials.annealfxn, smacprm, args=(smacnum, 25000, model, envlist, xpdata, 6, 3, lb, ub), lower=lower, upper=upper, full_output=1)
#started at 1:56a






#===
bidlist = ['kbidCbidMf','kbidCbidMr','kbidbcl2f','kbidbcl2r','kc8bidf','kc8bidr','kc8bidc']
bidnums = []
bidparms = []
for i in range(len(bidlist)):
    for j in range(len(model.parameters)):
        if bidlist[i] == model.parameters[j].name:
            bidparms.append(model.parameters[j].value)
            bidnums.append(j)
annprm = numpy.asarray(bidparms)
annnum = numpy.asarray(bidnums)
#==
annlout = scipy.optimize.anneal(pysb.anneal_sundials.annealfxn, smacprm, args=(smacnum, 16000, model, envlist, xpdata, , , lb, ub), lower=lower, upper=upper)
