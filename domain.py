
import sys
import os
import numpy

from binarydata import SpatialData
from fakedata import FakeSpatialData

#typedict = {'FEM': SpatialData, 'ODE':ODEData}


class ODEData(object):
    def __init__(self, path, domain, component):
        filename  = path + 'concentration.bin'
        nclusters = len(domain.simulation().clusters())
        #print filename
        #print 'filesize', os.stat(filename).st_size
        #print 'doubles in file', os.stat(filename).st_size/8
        #print 'frames in file',os.stat(filename).st_size/8/nclusters
        
        self.data = numpy.memmap(filename,
                                 dtype = [('t',numpy.float64), ('ca',(numpy.float64,(nclusters,)))],
                                 mode= 'r',
                                 shape = (os.stat(filename).st_size/8/(nclusters+1)))

class Domain(object):
    def __init__(self, path, name, simulation):
        self.__name       = name
        self.__simulation = simulation
        
        config = simulation.config
        
        if config.has_option(name,'Type'):
            self.__type = config.get(name,'Type')
            print "found", name, "domain with type", self.__type
        else:
            self.__type = 'FEM'
        
        
        if self.__type == 'ODE':
            DataType = ODEData
        else:
            DataType = SpatialData
        
        if config.has_option(name,'components'):
            for component in config.get(name,'components').split(','):
                #try:
                self.__dict__[component] = DataType(path, self, component)

               # except:
                #    print "Warning: could not load component:", name, component,sys.exc_info()[0]

    def components(self):
        return self.keys();
    
    def name(self):
        return self.__name
        
    def simulation(self):
        return self.__simulation
        
    def __getitem__(self,component):
        return self.__dict__[component]
