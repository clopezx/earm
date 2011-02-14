import scipy.optimize
import numpy
import csv



# anneal(func, x0, args=(), schedule='fast', full_output=0, T0=None, Tf=9.9999999999999998e-13, 
#        maxeval=None, maxaccept=None, maxiter=400, boltzmann=1.0, learn_rate=0.5, 
#        feps=9.9999999999999995e-07, quench=1.0, m=1.0, n=1.0, lower=-100, upper=100, dwell=50)


def read_cvs_array(fp):
    """return a numpy array from a cvs set of data"""
    reader = csv.reader(fp)
    for row in reader:
        
    converters = (float, str, int)
    for i, item in enumerate(data):
        print converters[i](item)
    
    
