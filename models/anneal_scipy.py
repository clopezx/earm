import scipy.optimize
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


def compare_data(array0, array1):
    """Compares two arrays and returns the X^2 between them"""
    # figure out the array shapes
    # this expects arrays of the form array([time, measurements])
    # the time is assumed to be roughly the same for both and the 
    # shortest time will be taken as reference to regrid the data
    # the regridding is done using a b-spline interpolation
    #

    # sanity checks
    # make sure we are comparing the right shape arrays
    arr0shape = array0.shape
    arr1shape = array1.shape
    
    if len(arr0shape) != len(arr1shape):
        raise SystemExit("comparing arrays of different dimensions")
    
    # get ref array for regridding...
    if arr0shape[0] > arr1shape[0]:
        refarray = array0
        otharray = array1
    else:
        refarray = array1
        otharray = array0

    

    
    
    
    
    
    
    
