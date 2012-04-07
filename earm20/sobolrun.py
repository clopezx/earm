from earm_2_0 import model
import pysb.varsens_sundials
import numpy
envlist, paramarr = pysb.varsens_sundials.odeinit(model)
xpfile = numpy.load('EC-RP_IMS-RP_IC-RP_data_for_models.npz')
xpnormdata = xpfile['arr_0']
sobolA = numpy.loadtxt('dim137_seed64_100pts.dat', delimiter=',')
sobolB = numpy.loadtxt('dim137_seed164_100pts.dat', delimiter=',')
sobprmsA = pysb.varsens_sundials.getsobolarr(sobolA, paramarr)
sobprmsB = pysb.varsens_sundials.getsobolarr(sobolB, paramarr)
sobprmsC = pysb.varsens_sundials.genCmtx(sobprmsA, sobprmsB)
yA, yB, yC = pysb.varsens_sundials.sobolfxn(model, sobprmsA, sobprmsB, sobprmsC, 20000, envlist, xpnormdata, [(2,1),(5,2),(8,3)], True, True, True)

