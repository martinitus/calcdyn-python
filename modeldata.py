import os
import numpy
import csv
import warnings
import itertools
import gzip

#from timeline import TimeLine
from channel import Channel
from cluster import Cluster

import ConfigParser

import dyk
import deterministic
import ryanodine

channelmodels = {}
channelmodels['DYKModel']      = dyk
channelmodels['Deterministic'] = deterministic
channelmodels['RyRModel']      = ryanodine
        
#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | {channel states} | {channel calciumlevels}
import contextlib

@contextlib.contextmanager
def patch_gzip_for_partial():
    """
    Context manager that replaces gzip.GzipFile._read_eof with a no-op.

    This is useful when decompressing partial files, something that won't
    work if GzipFile does it's checksum comparison.

    """
    _read_eof      = gzip.GzipFile._read_eof
    _add_read_data = gzip.GzipFile._add_read_data
    
    gzip.GzipFile._read_eof      = lambda *args, **kwargs: None
    gzip.GzipFile._add_read_data = lambda *args, **kwargs: None
    yield
    gzip.GzipFile._read_eof      = _read_eof
    gzip.GzipFile._add_read_data = _add_read_data

class EventData(object):
    
    def __init__(self, sim = None, path = None, refresh = None, transitions = False):
        if sim == None:
            assert(path)
            self.config = ConfigParser.RawConfigParser()
            self.config.read(path + "/parameters.txt")
        else:
            path = sim.path()
            self.config = sim.config
                
        # the amount of contained channels
        self._channelcount = self.config.getint('ChannelSetup','channels')
        self._clustercount = self.config.getint('ChannelSetup','clusters')
        
        # the event data (ModelData)
        modelname = self.config.get('Meta','channelmodel')
        
        self.__model = channelmodels[modelname]
        self.__transitions = None
    
        #~ self._states = self._model.states   
        # the set of defined states must at least contain a definition of open and closed state
        #~ assert(self._states.has_key('open'))
        #~ assert(self._states.has_key('closed'))
        
        # the channels 
        #~ channeldata = numpy.genfromtxt(path + '/channels.csv', dtype=[('id', int), ('cluster', int), ('location', float, (3)), ('radius', float)])
        # since genfromtxt returns 1d array for single line files we need to reshape
        #~ if channeldata.ndim == 0:
            #~ channeldata = [channeldata]
        
        #~ assert(self._channelcount == len(channeldata))
        
        def load_gzip(fp):
            f = gzip.open(fp)
            dtype  = numpy.dtype(self.__model.types_binary(self._channelcount,self._clustercount))
            #buff = f.read(dtype.itemsize * 100000)
            buff = f.read()
            return numpy.frombuffer(buff, dtype = self.__model.types_binary(self._channelcount,self._clustercount))

        def load_bin(fp):
            f = open(fp, 'rb')
            return numpy.fromfile(f, dtype = self.__model.types_binary(self._channelcount,self._clustercount))

        def load_csv(fp):
            tmp = self.__model.loadtxt(fp,self._channelcount)
            # kill all unnecessary stuff
            return tmp[['t','chid','clid','noch','nocl','states']]
        
        
        file_order = [('transitions.bin.gz', load_gzip),
                      ('transitions.bin',    load_bin),
                      ('transitions.csv',    load_csv),]
        
        for filename, openfunc in file_order:
            if os.path.isfile(path + '/' + filename):
                self._data = openfunc(path + '/' + filename)                
                break
        
        if transitions:
            selection = numpy.roll((self._data['noch'][1:]!=self._data['noch'][:-1]),1)
            selection[0] = True
            selection[-1] = False
            self._data = self._data[selection]
        
    def __repr__(self):
        return "EventData (Channels: %d, Events: %d)" % (self._channelcount, len(self._data))
        
    def _repr_svg_(self):
        return self.open()._repr_svg_()
        
    def channels(self):
        return self._channels
        
    def channel(self,i):
        return self._channels[i]
        
    def clusters(self):
        return self._clusters
        
    def cluster(self, i):
        return self._clusters[i]
        
    def model(self):
        return self.__model
        
    def openchannels(self):
        return self.observe(lambda x: x['noch'],desc = 'open')
        
    def closedchannels(self):
        return self.observe(lambda x: self._channelcount - x['noch'],desc = 'closed')
        
    def openclusters(self):
        return self.observe(lambda x: x['nocl'],desc = 'open')
        
    def closedclusters(self):
        return self.observe(lambda x: x['noch'],desc = 'open')
        
    # get a raw handle on the data
    def data(self):
        return self._data 
     
    def events(self):
        ''' return the events recarray'''
        return self._data 
        
    def transitions(self):
        '''return a recarray containing only the events where the number of open channels has changed''' 
        if self.__transitions == None:
            allevents = self.events()
            nopen = self.events()['noch']
            mask  = numpy.r_[1,numpy.diff(nopen)]
            self.__transitions = allevents[mask != 0]
        return self.__transitions
    
    def tmin(self):
        return self._data[0]['t']
    
    def tmax(self):
        return self._data[-1]['t']
    
    #~ return array with event times
    def times(self):
        return self._data['t'];
        
    # return a timeline object of the observable defined by the provided function
    def observe(self, observable, desc = None):
        return TimeLine(self._data[:]['t'], map(observable,self._data), interpolationorder = 'zero', desc = desc)
    
    #~ Calculate the intervalls the given condition evealuates to true, in either real time, or frames
    #~ a valid condition should do something like this:
    #~ def condition(data):
    #~      return np.logical_and.reduce([data[:,1] == 0, data[:,2] == 1])
    #~ and return an array of booleans indicating wether the condition is true or false for the given frame
    #~ the first and last intervall
    def intervalls(self, condition, frames = False):
        warnings.warn("deprecated use simulation.intervalls insead", DeprecationWarning)
        selection = condition(self._data)
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
            interv = numpy.array([self._data[start_end[:,0],0],self._data[start_end[:,1],0]]).transpose()
            return interv
        
    #~ return the time the given channel switched into the state of given time
    # use channel.state() instead
    '''def stateswitch(self,channel,time):
        f = self.frame(time)
        s = self.data[f,1+channel]
                
        #~ go back in time until we find a transition of the given channel
        while self.data[f,1+channel] == s:
            f = f - 1
            if f < 0:
                return 0        
        #~ at frame f, the channel was still in another state, hence it changed to its current state at frame f+1
        return self.data[f+1,0]'''
    
    #~ return state for given channel and given time
    # use channel.state() instead
    '''def state(self,channel,time):
        f = self.frame(time)
        return self.data[f,1+channel]
        '''
    # get the dictionary defining the given state   
    #def state(self,state):
    #   return self._states[state]
    
    # get the dictionary of predefined states
    #def states(self):
    #   return self._states
    '''def states(self,frame = None, time = None):
        if frame != None:
            assert(time == None)            
        if time != None:
            assert(frame == None)
            frame = self.frame(time)            
        return self._data[frame,1:1+self.channelcount()]'''
    
    def state(self,t):
        ''' return the state of the system for time t'''
        f = self._data['t'].searchsorted(t)
        return self._data[f]


