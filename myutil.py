import numpy

class AverageData(object):
    def __init__(self, path, datasets):
        self.data = {}
        for d in datasets:
            self.data[d] = numpy.genfromtxt(path + d + ".csv",dtype = [('t',float),('cyca',float),('mobile',float)])
	    
class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)