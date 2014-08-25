import numpy
import types

from cluster import Cluster

class Puff(object):
    def __init__(self, cluster, ts, te):
        self.__cluster = cluster
        
        # prepare buffers for views on transitions and events
        self.__events      = None
        self.__transitions = None
        
        self.__start   = ts
        self.__end     = te
        
        self.__model = cluster.model()
        #print "created puff", self.start(), self.end(), self.duration(), self.peak()
        #~for attribute in self.__model.collective_event_attributes():
            #~self.__dict__[attribute.name] = attribute.value(self.__events)
    
    
    def trajectory(self, dt, start = None, stop = None):
        start = self.start() if start == None else start;
        stop  = self.end()   if stop  == None else stop;
        for t in numpy.arange(start = start,stop = stop,step = dt):
            yield self.__events[self.__events['t'].searchsorted(t)]
    
    def peak(self):
        return self.transitions()['noch'].max()
        
    def events(self):
        if isinstance(self.__events, types.NoneType):
            e = self.__cluster.events()        
            firstframe = e['t'].searchsorted(self.__start)-1
            lastframe  = e['t'].searchsorted(self.__end)+1
            self.__events = e[firstframe:lastframe]
        return self.__events
        
    def transitions(self):
        if isinstance(self.__transitions, types.NoneType):
            e = self.__cluster.transitions()
            firstframe = e['t'].searchsorted(self.__start)-1
            lastframe  = e['t'].searchsorted(self.__end)+1
            self.__transitions = e[firstframe:lastframe]
        return self.__transitions
        
    def start(self):
        return self.__start
        
    def end(self):
        return self.__end
        
    def duration(self):
        return self.end()-self.start()
        
    def model(self):
        return self.__model
        
    # return the time integral over the number of open channels
    def accumulated(self):
        events = self.transitions()
        no = events['noch']
        t =  events['t']
        dts = (t[1:]-t[:-1])
        return numpy.dot(no[:-1],dts)
    
    def average(self):
        return self.accumulated()/self.duration()


#==============================================================
# Now, add the puffs method to the clusters

def remask_open_state(t, os, tolerance):
    assert(os[0] == False)
    assert(os[1] == True )
    
    converged = False
    
    # the counter of real closings
    count     = numpy.logical_not(os).sum()
    
    while not converged:
        # create mask for open and close events
        close_mask = numpy.r_[False,numpy.logical_and(os[:-1],
                                                  numpy.logical_not(os[1:]))]
                                                  
        open_mask  = numpy.r_[False,numpy.logical_and(numpy.logical_not(os[:-1]),
                                                  os[1:])]
        
        # extract the corresponding times
        close_times = t[close_mask]
        open_times  = t[open_mask]
        
        if (len(open_times) == len(close_times) + 1):
            open_times = open_times[:-1]
        
        # for each open time, check wether the difference to the last close time is less then tolerance
        reopen_event   = (open_times[1:] - close_times[:-1]) <= tolerance
        #print close_times
        #print (open_times[1:] - close_times[:1])[0:100]
        os[close_mask] = numpy.r_[reopen_event,False]
        
        converged = count == numpy.logical_not(os).sum()
        
        count = numpy.logical_not(os).sum()
    
    return os, open_mask, close_mask


# return iterable for iteration over puffs
def puffs(self, tolerance = 0.005):
    # add dictionary member if not already present
    if not '_Cluster__puffs' in self.__dict__:
        self.__dict__['_Cluster__puffs'] = {}
        
    # get handle on the member dictionary
    __puffs = self.__dict__['_Cluster__puffs']
    
    
    #if we do not have the data available for the asked tolerance, calculate it and store it in dictionary
    if not __puffs.has_key(tolerance):
        data = self.events()
        
        # get the transitions of this cluster
        tr = self.transitions()
        
        # calculate the openstate of the cluster,
        # a cluster is called open if at leas one channel is open
        # the resulting array is a mask, masking the times the cluster is open
        # with true, and false otherwise. since we only used the transitions to
        # calculate the mask, there will never be 2 consecutive closed states
        os = self.model().open(tr).any(axis = -1)
        
        # now, remask the open state to respect the tolerance 
        os,om,cm = remask_open_state(tr['t'],os,tolerance)
        
        # om and cm will mask the open and close events separatly
        open_indices, = om.nonzero() # needs unpacking, since nonzero returns tuple
        close_indices, = cm.nonzero()
        
        # create list of puffs and assign it to the entry of the dict
        #__puffs[tolerance] = [Puff(self,o,c) for (o,c) in zip(open_indices, close_indices)]
        open_times  = tr['t'][open_indices]
        close_times = tr['t'][close_indices]
        __puffs[tolerance] = [Puff(self,ot,ct) for (ot,ct) in zip(open_times, close_times)]

    
    # at this point we have calculated the new data, and also the puffs event data
    return __puffs[tolerance]
    
#~ TODO: create something like a puff collection to check wether a cluster is active given a certain tolerance
#~class PuffCollection(object):
    #~def __init__(self, puffs):
        #~self.__puffs = puffs
    #~
    #~def __iter__(self):
        #~for p in puffs:
            #~yield p
    #~
    #~def active
        

# add member method to all the puff class,
# all existing instances and all future instances
# will posses this method
Cluster.puffs = puffs
Cluster.collective_events = puffs

