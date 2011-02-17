import scipy.optimize
import numpy
import csv



# anneal(func, x0, args=(), schedule='fast', full_output=0, T0=None, Tf=9.9999999999999998e-13, 
#        maxeval=None, maxaccept=None, maxiter=400, boltzmann=1.0, learn_rate=0.5, 
#        feps=9.9999999999999995e-07, quench=1.0, m=1.0, n=1.0, lower=-100, upper=100, dwell=50)

def read_csv_array(fp):
    """return a numpy array from a cvs set of data"""
    reader = csv.reader(fp)
    templist = []
    #read in the lists
    for row in reader:
        templist.append(row)
    #remove empty spaces
    for i in range(0, len(templist)):
        templist[i] = filter(None, templist[i])
    #Now put these into a numpy array
    converters = (float, float, float, float)
    headstring = templist.pop(0)
    #assume all entries in the list are the same length
    darray = numpy.zeros((len(templist), len(templist[0])))
    for i, item in enumerate(templist):
        darray[i] = numpy.asarray(templist[i], dtype=darray.dtype)    
    return (headstring, darray)
