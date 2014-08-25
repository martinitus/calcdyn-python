import numpy
import scipy
import myutil

import dyk
import deterministic
from eventcollection import EventCollection

        

class Cluster(object):
    def __init__(self, sim, index , channels):
        self._channels = channels
        self.__eventdata = sim.events()
        self.__model = sim.events().model()
        self.__simulation = sim
        self.__events      = None
        self.__transitions = None
        self.__state       = None
        self.__open        = None
        self.__index  = index
        self.__puffs  = {}
    
    # deprecated
    def channels(self):
        return self._channels
    
    # make the cluster iterable over the channels
    def __iter__(self):
        for c in self._channels:
            yield c
    
    def model(self):
        return self.__model
        
    def index(self):
        return self.__index

    #~# return iterable for iteration over puffs
    #~def puffs(self, tolerance = 0.005):
        #~'''
            #~this method needs to be implemented by a specific puff criterion 
        #~'''
        #~raise Warning("Not implemented")
        
        
    def channel(self,i):
        return self._channels[i]
        
    def open(self,t):
        '''
         provide zero order interpolation of the channels open state for time(s) t.
         t can either be a scalar value, or a 1d array.
        '''
        if not self.__open:
            self.__open  = myutil.ZeroOrderExtrapolation(self.events()['t'], self.__model.open(self.events()['states']))
            #~self.__open = scipy.interpolate.interp1d(self.events()['t'], self.__model.open(self.events()['states']),axis = 0,kind = 'zero')
        return self.__open(t).astype(bool)
        
    def state(self,t):
        '''
         provide zero order interpolation of the channels state for time(s) t.
         t can either be a scalar value, or a 1d array.
        '''
        if not self.__state:
            self.__state  = myutil.ZeroOrderExtrapolation(self.events()['t'], self.events()['states'])
            #~self.__state = scipy.interpolate.interp1d(self.events()['t'], self.events()['states'],axis = 0,kind = 'zero')
        return self.__state(t)
    
    # return projection on the eventdata for this cluster
    def events(self):
        if isinstance(self.__events, types.NoneType):
            # the indices of the channels within this cluster
            cidx  = [channel.index() for channel in self]
            #print cidx
            # the eventdata from where to extract
            data  = self.__eventdata._data

            #select all events for this cluster
            eventmask = data['clid'] == self.__index
            #select first and last frame
            eventmask[0]  = True
            eventmask[-1] = True

            #create recarray that stores the events of this cluster
            if self.__model == dyk:                
                self.__events = numpy.recarray(shape = (eventmask.sum()),dtype = [('t', '<f8'), ('noch', '<i2'), ('chid', int), ('states', '|i1', (len(cidx), 8))])
            elif self.__model == deterministic:
                self.__events = numpy.recarray(shape = (eventmask.sum()),dtype = [('t', '<f8'), ('noch', '<i2'), ('chid', int), ('states', bool, (len(cidx),))])

            # copy time chid and subspace of state column to new recarray
            self.__events['t']       = data[eventmask]['t']
            self.__events['chid']    = data[eventmask]['chid']
            self.__events['states']  = data[eventmask]['states'][:,self.__index,...]
            
            # cache the number of open channels         
            model =  self.__eventdata.model()
            self.__events['noch'] = model.open(self.__events).sum(-1)
            
        return self.__events
        
    def transitions(self):
        '''return a recarray containing only the events where this cluster changes the number of open channels''' 
        if isinstance(self.__transitions, types.NoneType):
            allevents = self.events()
            nopen = self.__model.open(allevents).sum(-1)
            mask  = numpy.r_[1,numpy.diff(nopen)]
            self.__transitions = allevents[mask != 0]
        return self.__transitions
        
    def couplingcoefficient(self,available):
        ec = EventCollection(self.puffs(tolerance = 0.0)).filter(lambda x: x.start()>100 and dyk.available(x.events()[0]) == available)
        po = 1.*(ec['peak'] > 1).sum()/len(ec.events)
        c  = 1.-(1.-po)**(1./(available-1))
        #print "#available = ",available,":" , (ec['peak'] > 1).sum(),'/',len(ec.events),'=',round(po,2),'=> c=',round(c,2)
        return c