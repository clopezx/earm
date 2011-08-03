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

#==
In [16]: annlout
Out[16]: 
(array([  1.26999718e-02,   2.57964264e-02,   4.23375583e-02,
         1.38112174e-02,   9.98337006e-03,   1.70883682e-02,
         3.60034499e-08,   2.02760954e-03,   4.40485206e+00,
         1.67585963e-07,   2.92204241e-03,   6.77694844e-01,
         2.87273716e-04,   2.01920658e-03,   1.76225528e-05,
         7.66998716e-04,   1.13643222e-04,   2.20907044e-04,
         9.97888984e-06,   2.34302137e-03,   8.78400906e-06,
         5.71383079e-04,   4.56248219e-05,   1.94077303e-03,
         1.87039388e-06,   4.19845879e-03,   1.49570400e-06,
         3.69332249e-03,   2.89577036e-05,   1.00060936e-03,
         4.14389316e-05,   3.81191602e-04,   4.47199097e-05,
         5.79887793e-04,   2.29238269e-08,   9.85768999e-04,
         5.85529281e-08,   1.81442831e-03,   1.31363495e-07,
         4.71057338e-03,   7.97728052e-09,   2.21844981e-03]),
 1489.247030020693,
 6.1133390552906945,
 2301,
 45,
 613,
 5)
