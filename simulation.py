#import CalciumData
#import ModelData

import ConfigParser
import numpy
import math
import types

from scipy import integrate

#from timeline import TimeLine
from modeldata import EventData
from fakedata import FakeSpatialData
from membrane import Membrane
from channel import Channel
from cluster import Cluster
from domain import Domain
from binarydata import SpatialData

# hold the combination of event data, spatial data and channel data
class Simulation(object):
    def __init__(self, path,**kwargs):
        self.__path = path
        self.config = ConfigParser.RawConfigParser()
        self.config.read(path + "/parameters.txt")
        
        # the amount of contained channels
        self._channelcount = self.config.getint('ChannelSetup','channels')
        self._clustercount = self.config.getint('ChannelSetup','clusters')
        
        # the event data (ModelData)
        modelname = self.config.get('Meta','channelmodel')
        
        self._events    = EventData(sim = self,**kwargs)
        
        # the channels need to read the channel files themselves 
        self._channels = [Channel(self, index) for index in range(self._channelcount)]
            
        # the cluster data 
        self._clusters = [Cluster(self, index, [channel for channel in self._channels if channel.cluster() == index]) for index in range(self._clustercount)]
        
        self.__domains = {}
        for name in self.config.get('Meta','domains').split(','):
            self.__domains[name] = Domain(path, name, self)
            # add the domain as attribute to the simulation object
            self.__dict__[name] = self.__domains[name]
        
        # the membrane object
        self.__membrane = Membrane(self)
        
    def path(self):
        return self.__path
            
    def domains(self):
        return self.__domains
        
    def domain(self,domain):
        return self.__domains[domain]
        
    def property(self, section, item):
        return self.config.get(section,item)
        
    def membrane(self):
        return self.__membrane
    
    def events(self):
        return self._events 
    
    def channelcount(self):
        return len(self._channels)
        
    def channels(self):
        return self._channels
        
    def clusters(self):
        return self._clusters
        
    def channel(self,c):
        return self._channels[c]
        
    def cluster(self,c):
        return self._clusters[c]
        
    def tmin(self):
        return self.domain('cytosol')['calcium'].tmin()
        
    def tmax(self):
        return self.domain('cytosol')['calcium'].tmax()
        
    def totalfluorescence(self, data, fmin = None, fmax = None):
        assert(self.domain('cytosol').has_key('dye'))
        
        if fmin == None:
            fmin   =  self.config.getfloat('dye','fmin');
        if fmax == None:
            fmax   =  self.config.getfloat('dye','fmax');
            
        Bd     = self.config.getfloat('dye','B');

        return data*fmax + (Bd-data)*fmin
        
    def relativefluorescence(self,data, fmin = None, fmax = None):
        if fmin == None:
            fmin   =  self.config.getfloat('dye','fmin');
        if fmax == None:
            fmax   =  self.config.getfloat('dye','fmax');
        
        total  = self.totalfluorescence(data,fmin,fmax)
        Bd     = self.config.getfloat('dye','B');
        
        c0     = self.config.getfloat('cytosol','c0');
        kminus = self.config.getfloat('dye','kminus');
        kplus  = self.config.getfloat('dye','kplus');
        resting =  Bd * c0 / (kminus / kplus + c0);
        print "B",Bd
        print "c0",c0
        print "k-,k+,k-/k+",kminus,kplus,kminus/kplus
        print "b_rest:",resting
        print "(b_rest*fmax + (Bd_b-rest)*fmin)",(resting*fmax + (Bd-resting)*fmin)
        
        return total / (resting*fmax + (Bd-resting)*fmin) - 1
    
    
    # TODO: subclass numpy ndarray for return values that require invervalls functionality and add the intervalls method to this subclass
    def intervalls(self,t,d, condition, frames = False):
        selection = condition(d)
        #print 'selection',selection
        #selection = numpy.append(selection, [False])
        #~ if the condition evaluates to true for the last frame (selection[-1] == True and selection[0] == False) the following roll-xor combination will lead switch_frame[0] == True
        switch_frames = numpy.logical_xor(numpy.roll(selection, 1), selection)
        switch_frames[0] = selection[0] 
        switch_frames[-1] = False # always drop unfinished intervalls
        
        #print 'switchframes',switch_frames
        # detect where the the condition changes from true to false, the roll will directly mark the first and last frame where the condition is true
        start_end = switch_frames.nonzero()[0] # make the returned 0-dimensional tuple an array
        # we condition is true up to the end, we need drop the last transition to condition = true, since we cannot return a closed interval
        if start_end.shape[0] % 2 == 1:
            start_end = numpy.reshape(start_end[:-1],[start_end.size/2,2])                       # reshape the array to contain start-end pairs
        else:
            start_end = numpy.reshape(start_end,[start_end.size/2,2])                       # reshape the array to contain start-end pairs
        
        # always drop intervalls already started at t=0
        if selection[0]:
            start_end = start_end[1:]
            
        if frames:
            return start_end
        else:
            # return the intervalls where the condition is true
            return numpy.array([t[start_end[:,0]],t[start_end[:,1]]]).transpose()
