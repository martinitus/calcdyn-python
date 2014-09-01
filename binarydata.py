import os
import types
import numpy
import scipy.spatial
import scipy.integrate
import scipy.interpolate
import collections
from progressbar import ConsoleProgressBar as ProgressBar

from scipy.interpolate import LinearNDInterpolator,  interp1d

# get a chunk of frames from a binary file
def chunks(f,framesize,blocksize=1000):
    f.seek(0,os.SEEK_END)
    frames = f.tell()/4/framesize
    blocks = frames/blocksize
    
    f.seek(0,os.SEEK_SET)
    for i in range(blocks):
        f.seek(i*blocksize*framesize*4, os.SEEK_SET)
        yield numpy.fromfile(f,dtype=(numpy.float32,(framesize,)),count = blocksize)
        
    f.seek(blocks*blocksize*framesize*4,os.SEEK_SET)
    yield numpy.fromfile(f,dtype=(numpy.float32,(framesize,)),count = frames-(blocks+1)*blocksize-1)

# downsample binary dataset of new format
def downsample(path,dataset, force = False, verbose = False):
    if not force:
        if os.path.exists(path + dataset + ".downsampled.bin"):
            if os.path.getmtime(path + dataset + ".downsampled.bin") > os.path.getmtime(path + dataset + ".bin"):
                if verbose: print "No need to downsample " +path+dataset, "files are up2date..."
                return
    
    framefile = open(path + dataset + ".frames.downsampled.bin","wb")
    datafile  = open(path + dataset + ".downsampled.bin","wb")
    transfile = open(path + dataset + ".transposed.bin","wb")
    
    framesize = os.path.getsize(path + dataset + ".coordinates.bin")/4/3+1
#    framesize = 5248+1
    print "framesize ", framesize
    frames    = os.path.getsize(path + dataset + ".bin")/4/framesize
    print "frames ", frames
    f         = 0
    if verbose: print "Downsampling " + path + dataset
    progress = ProgressBar(frames)
    frames   = 0
    for chunk in chunks(open(path + dataset + ".bin","rb"),framesize):
        #print (chunk[1:,0])
        #print (chunk[1:,0]>=chunk[:-1,0])
        #assert((chunk[1:,0]>=chunk[:-1,0]).all())
        #mask = chunk[1:,0]-chunk[:-1,0]>=0.0001
        #print "masked:", len(mask) - mask.sum(),"frames"
        #data = chunk[mask]
        data = chunk
        
        frames = frames + len(data)
        #append the frames
        if not (data[1:,0]>=data[:-1,0]).all():
            fail = (data[1:,0]>=data[:-1,0]).all()
            print "Warning, the data seems to be corrupt: the frame information is not a rising list!",
            print "framesize:", framesize
            print "frames:", frames
            print data[fail,0]
            print data[:,0]
        #~assert()
        data[:,0].tofile(framefile)
        #append the data
        data[:,1:].tofile(datafile)
        #allocate space for the transposed data
        data[:,1:].tofile(transfile)
        f = f + len(chunk)
        progress(f)
    
    framefile.close()
    datafile.close()
    
    if verbose: print "Creating transposed dataset " + path + dataset
        
    assert(frames == os.path.getsize(path + dataset + ".frames.downsampled.bin")/4)
    
    progress = ProgressBar(frames)
    f      = 0
    # since we load junks from the downsampled file now, a frame does not contain the time column any more, hence, framesize is decreased by one
    for chunk in chunks(open(path + dataset + ".downsampled.bin","rb"),framesize-1):
        data = chunk.transpose() # no time column any more, use all data
        assert(len(data) == framesize - 1)
        for n in range(framesize-1):
            assert(len(data[n] == len(chunk)))
            transfile.seek((n*frames+f)*4,os.SEEK_SET)
            data[n].tofile(transfile)
        f = f+len(chunk)
        progress(f)

'''
def downsample(path,dataset,dt = 0.0006):
    print "Downsampling: " + path + "/" + dataset           
    
    floats = os.path.getsize(path + dataset + ".bin")/4
    frames = os.path.getsize(path + dataset + ".coordinates.bin")/4/3 + 1

    data = numpy.memmap(path + dataset + ".bin",dtype=numpy.float32,shape=(floats/frames,frames),mode = 'r')
    #extract frame time colums
    frames = data[:,0]

    print "calculating frame selection"
    selection = (frames[1:]-frames[:-1])>=dt
    selection[0] = True

    print "selected:",selection.sum()
    frames[selection].tofile(path + dataset + ".frames.downsampled.bin")

    # drop frame time column
    print "writing frame file"
    data[selection,1:].tofile(path+dataset+".downsampled.bin")

    print "Creating transposed binary data"
    swaped = data[selection,1:].swapaxes(0,1)
    swaped.tofile(path+dataset+".transposed.bin")
    print "Done"
'''


class SpatialData(object):  
    def __init__(self,path,domain,component,refresh=False):
        #~ TODO: make generic for 3 dimensions
        #~~self.nodes = numpy.zeros([0,3], dtype=numpy.float32)
        dataset = domain.name() + '.' + component
        
        #~ self.interpolator  = None
        self.triangulation = None
        
        # read the coordinates              
        self.__nodes = numpy.fromfile(open(path+dataset+".coordinates.bin"), dtype=(numpy.float32,(3,)))
        
        # if all z-coordinates are equal, assume a 2D dataset
        self.__2D = (self.__nodes[:,2] == self.__nodes[0,2]).all()
        #if self.__2D:
            # drop z-coordinate
            #self.__nodes = self.__nodes[:,0:2]
            
        print self.__2D
        
        if refresh:
            downsample(path,dataset)
        
        self.__frames = numpy.fromfile(path+dataset+".frames.downsampled.bin", dtype=numpy.float32)
        
        f = self.__frames.shape[0]
        n = self.__nodes.shape[0]
        
        self.__data       = numpy.memmap(path+dataset+".downsampled.bin",dtype = numpy.float32, shape = (f,n),mode = 'r')
        self.__transposed = numpy.memmap(path+dataset+".transposed.bin", dtype = numpy.float32, shape = (n,f),mode = 'r')
        
        # add calcium methods to channels
        self.__simulation = domain.simulation()
        
        attrib_name = '__' + component
        method_name = component
        
        def channel_concentration(channel,t = None):
            '''
                get the local calcium concentration of this channel.
                if no t given, return the tuple containing t,ca as numpy arrays
                if t given, return linear interpolated concentrations.
            '''
            if not attrib_name in channel.__dict__.keys():
                channel.__dict__[attrib_name] = self.evolution(channel.location()[0:2])
            
            if t != None:
                ip = scipy.interpolate.interp1d(channel.__dict__[attrib_name][0], channel.__dict__[attrib_name][1],axis = 0,kind = 'linear')
                return ip(t)
            else:
                return channel.__dict__[attrib_name]
            
        for ch in self.__simulation.channels():
            # bind the method to the instance only. this will avoid problems when different kinds of spatial data are loaded. 
            # and only some of the cluster instances should have the calcium method
            ch.__dict__[method_name] = types.MethodType( channel_concentration, ch )
    
    def select(self,xmin,xmax,ymin,ymax):
        selection = numpy.logical_and(numpy.logical_and(xmin<=self.nodes[:,0],self.nodes[:,0] <= xmax),numpy.logical_and(ymin<=self.nodes[:,1],self.nodes[:,1] <= ymax))
        self.nodes = self.nodes[selection]
        self.data  = self.data[:,selection]
        # reset triangulation
        self.triangulation = None
        
    #~ t and x can either be scalar or array types.
    #~ if x is None, the interpolation is done for all nodes
    #~ if t is None, the interpolation is done for all frames
    #~ if both, x and t are not None, the interpolation is done for all x and all t
    #~ if t is scalar, and either, x or y are arrays, the return type will be onedimensional
    def __call__(self,t = None, x = None):
        '''
         perform linear interpolation in space and time 
         If x is None, interpolate data for all nodes in time for the given times stored in t.
         If t is None, interpolate data for all frames in space for the given coordinates stored in x.
         
         If both, t and x are arrays, return an array with shape  (len(t), len(x))
         If either t or x is scalar, return an array with shape
        '''
    
        assert(t != None or x != None)
        
        if t == None:
            # TODO make this more efficient
            t = self.__frames
        if x == None:
            # TODO make this more efficient
            x = self.__nodes
        
        assert(isinstance( x, collections.Iterable ))
        
        # if t or x a scalar values, make them iterable
        x = [x] if not isinstance( x[0], collections.Iterable ) else x
        t = [t] if not isinstance( t, collections.Iterable ) else t
        
        # define the roi we need to load from the memmap
        
        lower_left  = numpy.min(x, axis = 0)
        upper_right = numpy.max(x, axis = 0)
        
        print x
        print lower_left, upper_right
        
        print (self.__nodes >= lower_left).all(axis = 1).sum()
        print (self.__nodes <= upper_right).all(axis = 1).sum()
        
        delauney = scipy.spatial.Delaunay(nodes, qhull_options= "Qbb Qc Qz QJ")
        
        simplices = delauney.find_simplex(x)
        
        # create a mask selecting all nodes within the roi
        mask = numpy.logical_and( (self.__nodes >= lower_left).all(axis = 1),
                                  (self.__nodes <= upper_right).all(axis = 1))
                                  
        mask = numpy.ones(shape= (len(mask),), dtype = bool)
        
        print "Found", mask.sum(), "coordinates within range"
            
        # find min and max frame we need the data for
        fmin  = self.__frames.searchsorted(min(t))-1
        fmax  = self.__frames.searchsorted(max(t))+1
        
        # select the nodes we need to interpolate from
        nodes = self.__nodes[mask]
        
        # depending of which dimesnion is larged, choose either frame ordered memmmap
        # or coordinate ordered memmmap
        if len(x) >= len(t):
            data = self.__transposed[mask,fmin:fmax]
        else:
            data = self.__data[fmin:fmax,mask]
        
        # TODO make this more efficient
        space_interpolator = LinearNDInterpolator(self.__nodes[mask], self.__data[mask,fmin:fmax])
        
        tmp = space_interpolator(x)
        
        time_interpolator  =  scipy.interpolate.interp1d(self.frames,tmp,copy = False)
        
        return time_interpolator(t)
        
        
        #interpolator  = LinearNDInterpolator(nodes[:,[0,1]],data)
        #interpolation = interpolator([[x,y]])[0]
        #return self.__frames[fmin:fmax], interpolation
        '''
        #~ print "x:",x
        #~ print "y:",y
        #~ print "t:",t
        #~ 
        #~ print "hasattr(x, 'shape'):",hasattr(x, 'shape')
        #~ print "hasattr(x, 'size'):",hasattr(x, 'size')
        
        datas = self.data.swapaxes(0,1)
        if hasattr(t, 'shape'):
            firstframe = self.frame(t[0])
            lastframe  = min(self.frames.size-1,self.frame(t[-1])+1)
        else:
            firstframe = self.frame(t)
            lastframe  = min(self.frames.size-1,firstframe+1)
        
        #create interpolator for spatial coordinates    
        interpolator  = scipy.interpolate.LinearNDInterpolator(self.nodes, datas[:,firstframe:lastframe+1],fill_value = fill_value)     
        
        # interpolate all nodes for given time and return scattered data
        if x == None and y == None and not hasattr(t, 'size'):
            interpolator = scipy.interpolate.interp1d(self.frames,datas,copy = False, fill_value = fill_value)
            return interpolator(t)
            
            
        elif hasattr(x, 'size') and hasattr(y, 'size'):
            assert(x.size == y.size)
            coordinates = numpy.zeros([x.size,2],dtype = numpy.float32)
            coordinates[:,0] = x[:]
            coordinates[:,1] = y[:]
        elif hasattr(x, 'size'):
            coordinates = numpy.zeros([x.size,2],dtype = numpy.float32)
            coordinates[:,0] = x[:]
            coordinates[:,1] = y
        elif hasattr(y, 'size'):
            coordinates = numpy.zeros([y.size,2],dtype = numpy.float32)
            coordinates[:,0] = x
            coordinates[:,1] = y[:]
        else:
            coordinates = [[x,y]]
        
        interpolation = interpolator(coordinates)
        print "interpolation:",interpolation
        
        timedomaininterpolator = scipy.interpolate.interp1d(self.frames[firstframe:lastframe+1],interpolation,copy = False, fill_value = fill_value)
        
        result = timedomaininterpolator(t)
        # if we had scalar x and y value, we want to return a onedimensional array
        if result.ndim == 2:
            return result[0,:]
        #~ print repr(result)
        return result'''
        
    # return a callable object representing a time interpolation for the given coordinates
    #~ def timeline(self,*args):
        #~ if len(args) == 1:
            #~ x=args[0][0]
            #~ y=args[0][1]
        #~ elif len(args)==2:
            #~ x=args[0]
            #~ y=args[1]
        #~ else:
            #~ raise Exception("dondt know what to do with:" + str(args))
        #~ datas = self.data.swapaxes(0,1)
        #~ interpolator  = scipy.interpolate.LinearNDInterpolator(self.nodes, datas[:,:])               
        #~ interpolation = interpolator([[x,y]])[0]
        #~ return timeline.TimeLine(self.frames,interpolation,ylabel=self.dataset,yunit='$\mu$M')
        
    def nodes(self):
        return self.__nodes
        
    def frames(self):
        return self.__frames
    
    def data(self, transposed = False):
        if transposed:
            return self.__transposed
        else:
            return self.__data
    
    
    def evolution(self,*args,**kwargs):
        ''' return a tuplple containing (frames, concentrations) for given coordinate'''
       
        if len(args) == 1:
            x=args[0][0]
            y=args[0][1]
        elif len(args)==2:
            x=args[0]
            y=args[1]
        else:
            raise Exception("dondt know what to do with:" + str(args))
            
        xmin = xmax = x
        ymin = ymax = y
        selection = numpy.ndarray([False])
        # create subset of data to avoid loading the whole memory map
        while selection.sum() < 5:
            xmin = xmin - 10
            xmax = xmax + 10
            ymin = ymin - 10 
            ymax = ymax + 10
            selectionx = numpy.logical_and(xmin<=self.__nodes[:,0],self.__nodes[:,0] <= xmax)
            selectiony = numpy.logical_and(ymin<=self.__nodes[:,1],self.__nodes[:,1] <= ymax)
            if not self.__2D:
                selectionz = self.__nodes[:,2] == 500.
                selection = numpy.logical_and(numpy.logical_and(selectionx,selectiony),selectionz)
            else:
                selection = numpy.logical_and(selectionx,selectiony)
            
        # take respect to tmin and tmax
        fmin  = 0 if not kwargs.has_key('tmin') else self.__frames.searchsorted(kwargs['tmin'])-1
        fmax  = len(self.__frames)-1 if not kwargs.has_key('tmax') else self.__frames.searchsorted(kwargs['tmax'])+1
        
        # print self.__nodes
        # print selection.sum()
        
        nodes = self.__nodes[selection]
        data  = self.__transposed[selection,fmin:fmax]
            
        interpolator  = LinearNDInterpolator(nodes[:,[0,1]],data)
        interpolation = interpolator([[x,y]])[0]
        return self.__frames[fmin:fmax], interpolation
        
    def linescan(self,tmin,tmax,xmin,xmax,y,dt=0.001,dx=1):
        sinterpolator  = LinearNDInterpolator(self.nodes(),self.data(transposed = True))
        #interpolate in space
        sinterpolation = sinterpolator([[x,y] for x in numpy.linspace(xmin,xmax, (xmax-xmin)/dx)])
        # regularize grid in time
        tinterpolator = interp1d(self.frames(),sinterpolation,axis = 1)
        return tinterpolator(numpy.linspace(tmin,tmax,(tmax-tmin) / dt)),numpy.linspace(tmin,tmax,(tmax-tmin) / dt),numpy.linspace(xmin,xmax, (xmax-xmin)/dx)
    
    
    # retrief a interpolator object for 2D slice at given time
    def snapshot(self,t):
        f0 = self.frame(t)
        f1 = f0 + 1
        
        t0 = self.__frames[f0]
        t1 = self.__frames[f1]
        
        # do manual interpolation to avoid loading the whole memmap
        return self.__data[f0] + (self.__data[f1]-self.__data[f0])*(t-t0)/(t1-t0)
    
    def write(self, path):
        # write coordinate file
        coordfile = open(path + '/coordinates.bin','w')
        self.nodes.tofile(coordfile)
        
        # write data file
        datafile = open(path + '/calcium.bin','w')
        self.data.tofile(datafile)
        
        # write frame file
        framefile = open(path + '/frames.bin','w')
        self.frames.tofile(framefile)
        
    def shift_coordinates(self,x,y):
        self.nodes[:,0] = self.nodes[:,0] + x
        self.nodes[:,1] = self.nodes[:,1] + y
        
    #return the frame number closest to the given time
    def frame(self, time):
        return numpy.searchsorted(self.__frames,time)
    
    def tmin(self):
        return self.__frames[0]
    
    def tmax(self):
        return self.__frames[-1]    
    
    # return index of node closest to [x,y]
    def node(self,*args):
        if len(args) == 1:
            x=args[0][0]
            y=args[0][1]
        elif len(args)==2:
            x=args[0]
            y=args[1]
        else:
            raise Exception("dondt know what to do with:" + str(args))
        n = 0
        mn = 0
        mindist = 9E30
        for node in self.__nodes:
            #~ print node
            dist = (node[0] - x)*(node[0] - x) + (node[1] -y)*(node[1] -y)
            if dist < mindist:              
                mindist = dist
                mn = n
                #~ print dist, node, mn
            n = n + 1
        return mn
        
        
    def grid(self, time, xmin,xmax,ymin,ymax,resolution):
        xi = numpy.linspace(xmin,xmax,resolution)
        yi = numpy.linspace(ymin,ymax,resolution)
        dim3 = True
        if not self.__2D:            
            zslize = 500.
            selection = self.__nodes[:,2] == zslize
            nodes = self.__nodes[selection,0:2]
            snapshot = self.snapshot(time)
            snapshot = snapshot[selection]
        else:
            nodes = self.__nodes[:,0:2]
            snapshot = self.snapshot(time)
            
        zi = scipy.interpolate.griddata(nodes, snapshot, (xi[None,:], yi[:,None]), method='cubic')
        return zi
        
    def extend(self):
        return numpy.array([self.__nodes[:,0].min(),self.__nodes[:,0].max(),self.__nodes[:,1].min(),self.__nodes[:,1].max()])
    
    def center(self):
        e  = self.extend()      
        return numpy.array([e[1]/2,e[3]/2])