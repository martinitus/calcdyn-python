import numpy
import os
import ConfigParser

import types
import scipy.interpolate

from cluster import Cluster

class ODEData(object):
    def __init__(self, path, domain, component):
        path = path if path[-1] == '/' else path + '/'
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
                                 
        
        config = ConfigParser.RawConfigParser()
        config.read(path + "/parameters.txt")
        self.__c1 = config.getfloat('Calcium','c1')
        
        self.__simulation = domain.simulation()
        
        # add calcium methods to channels and clusters
        
        def cluster_calcium(cluster, t = None):
            '''
                get the local calcium concentration of this cluster.
                if no t given, return the tuple containing t,ca as numpy arrays
                if t given, return linear interpolated concentrations
            '''
            if t == None:
                return self.evolution(cluster = cluster)        
            else:
                ip = scipy.interpolate.interp1d(*self.evolution(cluster = cluster), axis = 0,kind = 'linear')
                return ip(t)
                
                
        def cluster_bulk(cluster, t = None):
            '''
                get the bulk calcium concentration of this cluster.
                if no t given, return the tuple containing t,ca as numpy arrays
                if t given, return linear interpolated concentrations
            '''
            if t == None:
                return self.bulk(cluster)        
            else:
                ip = scipy.interpolate.interp1d(*self.bulk(cluster), axis = 0,kind = 'linear')
                return ip(t)
                
        
        for cl in self.__simulation.clusters():
            # bind the method to the instance only. this will avoid problems when different kinds of spatial data are loaded. 
            # and only some of the cluster instances should have the calcium method
            cl.calcium = types.MethodType( cluster_calcium, cl )
            cl.bulk    = types.MethodType( cluster_bulk, cl )
  
        
        def channel_calcium(channel,t = None):
            '''
                get the bulk calcium concentration of this cluster.
                if no t given, return the tuple containing t,ca as numpy arrays
                if t given, return linear interpolated concentrations
            '''
            if t == None:
                return self.evolution(channel.cluster())        
            else:
                ip = scipy.interpolate.interp1d(*self.evolution(channel.cluster()), axis = 0,kind = 'linear')
                return ip(t)
                
        for ch in self.__simulation.channels():
            # bind the method to the instance only. this will avoid problems when different kinds of spatial data are loaded. 
            # and only some of the cluster instances should have the calcium method
            ch.calcium = types.MethodType( channel_calcium, ch )
                                 
    def evolution(self,cluster):
        ''' return a tuplple containing (frames, concentrations) for the given cluster'''        
        idx = cluster if isinstance(cluster, int) else cluster.index()        
        return self.data['t'], self.data['ca'][:,idx]
        
    def bulk(self, cluster):
        ''' calculate the bulk concentration from the global concentration, return tuple containing t and b '''
        b = self.data['ca'][:,cluster.index()]
        b = b - cluster.open(self.data['t']).sum(-1) * self.__c1
        return self.data['t'], b