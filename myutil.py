import numpy

class AverageData(object):
    def __init__(self, path, datasets):
        self.data = {}
        for d in datasets:
            self.data[d] = numpy.genfromtxt(path + d + ".csv",dtype = [('t',float),('cyca',float),('mobile',float)])
	    
class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
        
        
# Necessary modules


import scipy.interpolate

class ZeroOrderExtrapolation(object):
    def __init__(self,x,y):
        self.__interp1d = scipy.interpolate.interp1d(x,y,axis = 0,kind = 'zero',bounds_error = False)
        self.__x = x
        self.__y = y
    
    def __call__(self,xo):
        ''' xo the ouput range where to interpolate for '''
    
        # the output shape
        os = (len(xo),) + self.__y.shape[1:]
                
        # Generate an empty output array for "y" values
        yo = numpy.empty(shape = os,dtype = self.__y.dtype)

        # Values lower than the minimum "x" are extrapolated at the same time
        low = xo < self.__x[0]
        yo[low] =  self.__y[0]

        # Values higher than the maximum "x" are extrapolated at same time
        high = xo > self.__x[-1]
        yo[high] = self.__y[-1]

        # Values inside the interpolation range are interpolated directly
        inside = numpy.logical_and(xo >= self.__x[0], xo <= self.__x[-1])
        
        indices = self.__x.searchsorted(xo[inside],side = 'right') - 1
        
        yo[inside] = self.__y[indices]
        #~print self.__interp1d(xo[inside]).dtype
        #~print yo.dtype
        #~yo[inside] = self.__interp1d(xo[inside]).astype(yo.dtype)
        return yo