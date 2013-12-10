import dyk
import numpy


''' statistical tool, needs further improvements'''
class EventCollection(object):
    def __init__(self, *args):
        self.events = [];
        for arg in args:
            self.events = self.events+arg
        
        self.data = numpy.ndarray(shape = (len(self.events)),dtype = [('duration',float),('peak',int),('accumulated',float),('available@start',int),('available@end',int),('withip3@start',int),('withip3@end',int)])
        for (e,i) in zip(self.events,range(len(self.events))):
            self.data[i] = (e.duration(),e.peak(),e.accumulated(),dyk.available(e.events())[0],dyk.available(e.events())[-1],dyk.withip3(e.events()[0]),dyk.withip3(e.events()[-1]))
        
    def filter(self,condition):
        return EventCollection([e for e in self.events if condition(e)])
    
    def size(self):
        return len(self.events)
        
    #def hist(self, key, normed = False):
    #    if normed:
    #    else:
    #        return numpy.histogram()
    
    def __getitem__(self,key):
        if type(key) == str:
            return self.data[key]
        elif type(key) == int:
            return self.events[key]
        raise "Invalid key type" + str(type(key))
        
    def select(self, property, condition):
        return [property(event) for event in self.events if condition(event)]
        
    def count(self, condition):
        return sum(1 if condition(element) else 0 for element in self.events)
        
    def average(self, property, condition = lambda x: True):
        s = sum(property(element) for element in self.events if condition(element))
        return s/count(condition)
        
    #def distribution(self, condition)
        
        
    def trajectories(self, progressbar = None):
        H = numpy.ndarray((17,17,17),dtype=int)
        H[:,:,:] = 0
        pc = np.ndarray((17),dtype = int)
        pc[:] = 0
        samples = 200000
        dt = 950./samples
        for puff in self.events:
            z = puff.available()
            pc[z] = pc[z]+1
            for frame in puff.trajectory(dt = 0.001):
                # calculate the number of available subunits (i.e. subunits that have ip3 bound)
                withip3    = frame['states'][:,[dyk.X100,dyk.X110,dyk.X111,dyk.X101]].sum(axis=1)
                inhibited  = frame['states'][:,[dyk.X101,dyk.X111]].sum(axis=1)
                        
                # the number of open channels
                x = (frame['states'][:,dyk.X110]>=3).sum()
                # the number of inhibited channels (inhibited means, the opening is blocked by inhibition and not ip3 unbinding)
                y = (numpy.logical_or(numpy.logical_and(withip3 == 3,inhibited == 1),numpy.logical_and(withip3 == 4,inhibited == 2))).sum()
                
                H[x,y,z] = H[x,y,z]+1
                progressbar.animate(frame['t'])