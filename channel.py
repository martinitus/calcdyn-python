import math
import numpy
import scipy

#from timeline import TimeLine

class Channel(object):
    ''' 
        infos:       a tuple containing channel id, channel location and cluster id
        transitions: a list of channel transitions, the channel will extract its own transitions from the global transition list
    '''
    def __init__(self,sim,index):
        from StringIO import StringIO
        line = open(sim.path() + '/channels.csv', "r").readlines()[index]
        
        self.__simulation = sim
        self.__eventdata  = sim.events()
        self.__model      = sim.events().model()                
        self.__events     = None
        self.__index      = index
        self.__transitions= None
        self.__calcium    = None
        self.__state      = None
        self.__open       = None
        
        if sim.config.has_option('Meta','ChannelType'):
            ct = sim.config.get('Meta','ChannelType')
            if ct == 'SimpleChannel':
                infos = numpy.genfromtxt(StringIO(line),dtype=[('id', int), ('cluster', int),('state',self.__model.state_type())])
                assert(self.__index == int(infos['id']))
        else:
            infos = numpy.genfromtxt(StringIO(line),dtype=[('id', int), ('cluster', int), ('location', float, (3)), ('radius', float)])
            assert(self.__index == int(infos['id']))
            self.__location   = infos['location']      
            self.__radius     = float(infos['radius'])
            
        
        self.__cluster    = int(infos['cluster'])
        
    def model(self):
        return self.__model

    def events(self):
        '''return a rec array with all events of this channel'''
        if self.__events == None:           
            #the event data from where to extract
            data  = self.__eventdata._data
            
            #select all events for this cluster
            eventmask = data['chid'] == self.__index
            
            #select first and last frame
            eventmask[0]  = True
            eventmask[-1] = True
            
            #create recarray that stores the events of this channel
            # this works for DYK
            # self.__events = numpy.recarray(shape = (eventmask.sum()),dtype = [('t', '<f8'), ('states', '|i1', 8)])
            self.__events = numpy.recarray(shape = (eventmask.sum()),dtype = [('t', '<f8'), ('states', self.__model.state_type())])
            
            # copy time chid and subspace of state column to new recarray
            self.__events['t']       = data[eventmask]['t']
            # this works for DYK
            #self.__events['states']  = data[eventmask]['states'][:,self.__index,:]
            self.__events['states']  = data[eventmask]['states'][:,self.clusterindex(),self.__index]
            
        return self.__events

    def transitions(self):
        '''return a recarray containing only the events where this channel either openes or closes''' 
        if self.__transitions == None:
            allevents = self.events()
            nopen = self.__model.open(allevents)
            mask  = numpy.r_[1,numpy.diff(nopen)]
            self.__transitions = allevents[mask != 0]
        return self.__transitions

    def ioi(self, absolute = False):
        '''
            return a list of inter open intervals for this channel, i.e. the intervals this channel is closed
            if absolute is True, return two list, first containing start times, second the durations of the intervalls
            if absolute is False(default) return a 1D list of intervalls
        '''
        # skip the first frame, since it might be the
        # initial state at t=0 and then it can corrupt statistics
        tr = self.transitions()[1:]
        # if the first time corresponds to a channel opening then we need to drop another frame
        if self.__model.open(tr['states'][0]):
            tr = tr[1:]
            
        if not absolute:
            return (tr['t'][1:]-tr['t'][:-1])[::2]
        else:
            return tr['t'][::2], (tr['t'][1:]-tr['t'][:-1])[::2]

    def ici(self):
        '''return a list of inter closed intervals for this channel, i.e. the intervals this channel is open'''
        tr = self.transitions()[1:]
        if not self.__model.open(tr['states'][0]):
            tr = tr[1:]
        return (tr['t'][1:] - tr['t'][:-1])[::2]
    


    def radius(self):
        return self.__radius
        
    def location(self):
        return self.__location
        
    def index(self):
        return self.__index
    
    def cluster(self):
        return self.__simulation.cluster(self.__cluster)
        
    def clusterindex(self):
        return self.__cluster
        
    def open(self,t):
        '''
         provide zero order interpolation of the channels open state for time(s) t.
         t can either be a scalar value, or a 1d array.
        '''
        if not self.__open:
            self.__open = scipy.interpolate.interp1d(self.events()['t'], self.__model.open(self.events()['states']),axis = 0,kind = 'zero')
        return self.__open(t).astype(bool)
        
    def state(self,t):
        '''
         provide zero order interpolation of the channels state for time(s) t.
         t can either be a scalar value, or a 1d array.
        '''
        if not self.__state:
            self.__state = scipy.interpolate.interp1d(self.events()['t'], self.events()['states'],axis = 0,kind = 'zero')
        return self.__state(t)
        
    def fluxcoefficient(self):
        if self.__simulation.config.has_option('membrane','pc'):
            return self.__simulation.config.getfloat('membrane','pc')
        else:
            ic   = self.__simulation.config.getfloat('membrane','ic')
            r    = self.__radius
            erc0 = self.__simulation.config.getfloat('er','c0')
            cyc0 = self.__simulation.config.getfloat('cytosol','c0')
            e    =  1.602176E-19
            Na   =  6.022141E23
            return ic * 1.E18 / (2 * r * r * math.pi * e * Na * (erc0 - cyc0))
        
    '''def current(self,t):
        ersnapshot = self.__simulation.domain('er')['calcium'].snapshot(t)
        cysnapshot = self.__simulation.domain('cytosol')['calcium'].snapshot(t)
        x,y=self.__location[0:2]
        fa = 96485 # faraday number
        foo = integrate.quad(lambda r: 2*math.pi*r*(ersnapshot([[x+r,y]])[0] - cysnapshot([[x+r,y]])[0]),0,self.__radius)
        return foo[0]*(1E-6*1E-24)*2*fa*self.fluxcoefficient()*1E12'''
        
    def current_integrate(self):
        x,y=self.__location[0:2]
        fa = 96485 # faraday number
        n  = len(self.__simulation.domain('cytosol')['calcium'].frames)
        foo = numpy.ndarray(n,dtype=numpy.float32)
        #print 'steps:', n
        for i in range(n):
        #   print i
            t = self.__simulation.domain('cytosol')['calcium'].frames[i]
            ersnapshot = self.__simulation.domain('er')['calcium'].snapshot(t)
            cysnapshot = self.__simulation.domain('cytosol')['calcium'].snapshot(t)
            if self.open(t):
                foo[i] = integrate.quad(lambda r: 2*math.pi*r*(ersnapshot([[x+r,y]])[0] - cysnapshot([[x+r,y]])[0]),0,self.__radius)[0]
            else:
                foo[i] = 0
        foo = foo*(1E-6*1E-24)*2*fa*self.fluxcoefficient()*1E12
        return self.__simulation.domain('cytosol')['calcium'].frames, foo
        
    def current(self):
        x,y= self.__location[0:2]
        r  = self.__radius
        fa = 96485 # faraday number
        t1,cy = self.__simulation.domain('cytosol')['calcium'].evolution(x,y)
        t2,er = self.__simulation.domain('er')['calcium'].evolution(x,y)
        
        f,s = self.state()
        
        import dyk
        mask = self.__model.open(s)
        current = r*r*math.pi*2*fa*self.fluxcoefficient()*(1E-6*1E-24)*1E12 * (er-cy)*mask
        
        return t1,current
        
    def calcium(self, t = None):
        '''
            get the local calcium concentration of this channel.
            if no t given, return the tuple containing t,ca as numpy arrays
            if t given, return linear interpolated concentrations.
            This function needs to be implemented by the caclium data
            that is attached to the simulation.
        '''
        raise Warning("Not implemented")
        
    '''def current_fast(self):  
        x,y    = self.__location[0:2]
        frames = self.__simulation.domain('cytosol')['calcium'].frames
        data   = 
        cur    = numpy.ndarray(len(frames),dtype=numpy.float32)
        for i in range(len(frames)):
            interpolator = scipy.interpolate.LinearNDInterpolator(nodes, data[7221,:],fill_value = numpy.nan)   
        datas         = self.data.swapaxes(0,1)
        
        # create time domain interpolator
        timeinterpol  = scipy.interpolate.interp1d(self.frames,datas,copy = False, fill_value = numpy.nan)
        
        # calculate scattered data values for time t
        datat         = timeinterpol(t)
        
        #create interpolator for spatial coordinates    
        return scipy.interpolate.LinearNDInterpolator(self.nodes, datat,fill_value = numpy.nan)     '''
