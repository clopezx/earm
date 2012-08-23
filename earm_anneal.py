import numpy as np
import pylab as pl
import pysb.anneal_vode
import pysb.anneal as anneal


def earmanneal(model, obj='m1'):
    # get the index for init concentrations, those won't change during run
    initidx = []
    for i,j in enumerate(model.parameters):
        if "_0" in j.name:
            initidx.append(i)

    # now set the parameter bounds, get the param array for annealfxn
    # an omag range of .3 is roughly equivalent to a parameter search in
    # the range [.5x, 2x]. An omag of 1 is equivalent to a parameter search
    # in the range [.1x, 10x] (one omag lower and higher). An omag of 0 does
    # not allow any changes in the initial value
    lb, ub, paramarr = pysb.anneal_vode.logparambounds(model, omag=.3,
                                                       useparams = initidx,
                                                       usemag = 0)
    
    xpfile = np.load('xpdata/forfits/EC-RP_IMS-RP_IC-RP_data_for_models.npz')
    xpnormdata = xpfile['arr_0']

    if obj == 'm1':
        annlout = anneal.anneal(pysb.anneal_vode.annealfxn, paramarr, \
                                args=(20000, model, xpnormdata, \
                                      [(2,1), (5,2), (8,3)], lb, ub), \
                                lower = 0.0, upper = 0.01, full_output=1)
    elif obj == 'm2':
        annlout = anneal.anneal(pysb.anneal_vode.annealfxn, paramarr, \
                                args=(20000, model, xpnormdata, \
                                      [(2,1), (8,3)], lb, ub, \
                                      [2,9810,7245000,180,3600]), \
                                lower = 0.0, upper = 0.01, full_output=1)
    else:
        print "objective function method not implemented"
        raise ValueError("Invalid Value for obj:",obj)
        
    # process updated fitted values
    paramarr = pysb.anneal_vode.mapprms(annlout[0], lb, ub)
    for i,j in enumerate(model.parameters):
        j.value = paramarr[i]
    return paramarr

def graphdata(model):
        xpfile = np.load('xpdata/forfits/EC-RP_IMS-RP_IC-RP_data_for_models.npz')
        xpnormdata = xpfile['arr_0']

        t = np.linspace(0,20000, 500)
        y = pysb.integrate.odesolve(model, t)

        pl.ion()
        pl.figure()
        pl.subplot(311)
        pl.plot(xpnormdata[0], xpnormdata[2])
        pl.plot(t, y.mBid/max(y.mBid))
        pl.subplot(312)
        pl.plot(xpnormdata[0], xpnormdata[5])
        pl.plot(t, y.aSmac/max(y.aSmac))
        pl.subplot(313)
        pl.plot(xpnormdata[0], xpnormdata[8])
        pl.plot(t, y.cPARP/max(y.cPARP))
        
        
        
        return 0
