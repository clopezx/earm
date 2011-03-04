import scipy.optimize
import scipy.interpolate
import numpy
import csv

def read_csv_array(fp):
    """returns the first string and a numpy array from a csv set of data"""
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


def compare_data(array0, array0var=None, array1):
    """Compares two arrays and returns the X^2 between them"""
    # figure out the array shapes
    # this expects arrays of the form array([time, measurements])
    # the time is assumed to be roughly the same for both and the 
    # shortest time will be taken as reference to regrid the data
    # the regridding is done using a b-spline interpolation
    # array0var shuold be the variances at every time point
    #

    # sanity checks
    # make sure we are comparing the right shape arrays
    arr0shape = array0.shape
    arr1shape = array1.shape
    
    if len(arr0shape) != len(arr1shape):
        raise SystemExit("comparing arrays of different dimensions")
    
    # get the time range where the arrays overlap
    rngmin = max(array0[:,0].min(), array1[:,0].min)
    rngmax = min(array0[:,0].max(), array1[:,0].max)
    rngmin = round(rngmin, -1)
    rngmax = round(rngmax, -1)
    
    # use the experimental gridpoints from the reference array as
    # the new gridset for the model array. notice the time range
    # of the experiment has to be within the model range
    #
    iparray = numpy.zeros(array0.shape)
    iparray[:,0] =  array0[:,0]
    
    # now create a b-spline of the data and fit it to desired range
    tck = scipy.interpolate.splrep(array1[:,0], array1[:,1])
    iparray[:,1] = scipy.interpolate.splev(iparray[:,0], tck)
    
    # we now have x and y axis for the points in the model array
    # calculate the objective function
    
    diffarray = array0[:,1] - iparray[:,1]
    diffsqarray = diffarray * diffarray
    
    # assume a default .05 variance
    if array0var is None:
        array0var = numpy.ones(iparray[:,1].shape)
        array0var = array0var*.05
    
    array0var = numpy.multiply(array0var,array0var)
    array0var = array0var*0.5

    objarray = diffsqarray * array0var
    
    return objarray.sum()
    
    
    
    
    
    
