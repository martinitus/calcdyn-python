import os
import numpy
import scipy.spatial
import scipy.integrate
import scipy.interpolate
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
        
    f.seek(blocks*blocksize*4,os.SEEK_SET)
    yield numpy.fromfile(f,dtype=(numpy.float32,(framesize,)),count = frames-(blocks+1)*blocksize-1)

# downsample binary dataset of new format
def downsample(path,dataset, force = False, verbose = True):
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
        assert((data[1:,0]>=data[:-1,0]).all())
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
	def __init__(self,path,dataset,refresh=False):
		#~ TODO: make generic for 3 dimensions
		#~~self.nodes = numpy.zeros([0,3], dtype=numpy.float32)
		self.dataset = dataset
		
		#~ self.interpolator  = None
		self.triangulation = None
		
		# read the coordinates				
		self.__nodes = numpy.fromfile(open(path+dataset+".coordinates.bin"), dtype=(numpy.float32,(3,)))
		# drop z-coordinate
		self.__nodes = self.__nodes[:,0:2]
		
		if refresh:
			downsample(path,dataset)
		
		self.__frames = numpy.fromfile(path+dataset+".frames.downsampled.bin", dtype=numpy.float32)
		
		f = self.__frames.shape[0]
		n = self.__nodes.shape[0]
		
		self.__data       = numpy.memmap(path+dataset+".downsampled.bin",dtype = numpy.float32, shape = (f,n),mode = 'r')
		self.__transposed = numpy.memmap(path+dataset+".transposed.bin", dtype = numpy.float32, shape = (n,f),mode = 'r')
	
	def select(self,xmin,xmax,ymin,ymax):
		selection = numpy.logical_and(numpy.logical_and(xmin<=self.nodes[:,0],self.nodes[:,0] <= xmax),numpy.logical_and(ymin<=self.nodes[:,1],self.nodes[:,1] <= ymax))
		self.nodes = self.nodes[selection]
		self.data  = self.data[:,selection]
        # reset triangulation
		self.triangulation = None
		
	#~ t, x, and y can all be either scalar or array types, if x and y are vectors, they both need to be the same size
	#~ if all three are vector types, the return value will be a two dimensional array containing frames in first, and spatial values in second dimension
	#~ if t is scalar, and either, x or y are arrays, the return type will be onedimensional
	'''def __call__(self,t = None, x = None,y = None,fill_value = numpy.nan):
		if self.triangulation == None:
			self.triangulation = scipy.spatial.Delaunay(self.nodes) 
		
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
	
	# return a frames and concentrations for given coordinate
	def evolution(self,*args,**kwargs):
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
			selection = numpy.logical_and(numpy.logical_and(xmin<=self.__nodes[:,0],self.__nodes[:,0] <= xmax),numpy.logical_and(ymin<=self.__nodes[:,1],self.__nodes[:,1] <= ymax))
			
		# take respect to tmin and tmax
		fmin  = 0 if not kwargs.has_key('tmin') else self.__frames.searchsorted(kwargs['tmin'])-1
		fmax  = len(self.__frames)-1 if not kwargs.has_key('tmax') else self.__frames.searchsorted(kwargs['tmax'])+1
		
		nodes = self.__nodes[selection]
		data  = self.__transposed[selection,fmin:fmax]
			
		interpolator  = LinearNDInterpolator(nodes,data)
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
		zi = scipy.interpolate.griddata(self.__nodes, self.snapshot(time), (xi[None,:], yi[:,None]), method='cubic')
		return zi
		
	def extend(self):
		return numpy.array([self.__nodes[:,0].min(),self.__nodes[:,0].max(),self.__nodes[:,1].min(),self.__nodes[:,1].max()])
	
	def center(self):
		e  = self.extend()		
		return numpy.array([e[1]/2,e[3]/2])
		
class RankData(SpatialData):
	two_dimensions = False
	
	def __init__(self, path, dataset, rank, verbose = False):
		super(RankData, self).__init__(dataset)
		
		# read the coordinates				
		coordfile = open(RankData.coordfile(path,dataset,rank))					
		self.nodes = numpy.fromfile(coordfile, dtype=numpy.float32)
		
		if(RankData.two_dimensions):
			#print "using 2D dataset", self.nodes.size
			self.nodes = numpy.reshape(self.nodes,(self.nodes.size / 2,2))			
		else:
			#print "using 3D dataset", self.nodes.size
			self.nodes = numpy.reshape(self.nodes,(self.nodes.size / 3,3))
			self.nodes = self.nodes[:,0:2]
			

		# read the data file in 10mb junks
		datafile = open(RankData.datafile(path,dataset,rank))
		floats_per_block = 10*1024*1024/4
		floats_in_file   = os.path.getsize(RankData.datafile(path,dataset,rank))/4
		self.data        = numpy.zeros(floats_in_file,dtype = numpy.float32)
		block            = 0
		blocks           = floats_in_file / floats_per_block
		
		# print 'fpb=%i, fif=%i, blocks=%i' % (floats_per_block,floats_in_file,blocks)	
		
		# read blocks
		while block < blocks:
			if verbose:
				print 'Reading {path}/calcium{rank}.bin {read:>11}/{work:>11} {progress:<.0%}'.format(path=os.path.relpath(path), rank = rank, read = block*floats_per_block, work = floats_in_file, progress = block*floats_per_block*1./floats_in_file) + '\r',
			self.data[block*floats_per_block:(block+1)*floats_per_block] = numpy.fromfile(datafile, dtype=numpy.float32, count=floats_per_block)
			block = block + 1		
			
		# read remaining floats
		self.data[blocks*floats_per_block:] = numpy.fromfile(datafile, dtype=numpy.float32, count = floats_in_file - blocks*floats_per_block)
		if verbose:
			print 'Reading {path}/calcium{rank}.bin {read:>11}/{work:>11} {progress:<.0%}'.format(path=os.path.relpath(path), rank = rank, read = floats_in_file, work = floats_in_file, progress = 1.),
		
		# determine number of frames
		f = self.data.size / (self.nodes.shape[0]+1)
		
		foo = ' '
		# reshape to 2D array  containing (frames x time+nodes)		
		if not self.data.size % (self.nodes.shape[0] + 1) == 0:
			foo = '[corrected corrupted file]'
			self.data.resize(f*(self.nodes.shape[0] + 1))
			
		self.data = numpy.reshape(self.data,(f, self.nodes.shape[0] + 1))		
		# extract frames from array
		self.frames = self.data[:,0]
		
		# drop frame time column
		self.data = self.data[:,1:]
		if verbose:
			print '   t=[{start:.8f}, {end:.8f}], {frames} frames {flag}'.format(start=float(self.frames[0]), end = float(self.frames[-1]), frames = self.frames.size, flag = foo)
			
			
	@staticmethod
	def coordfile(path, dataset, rank):
		if (os.path.exists(os.path.join(path,  dataset + '.coordinates.' + str(rank) + '.bin'))):
			return os.path.join(path,  dataset + '.coordinates.' + str(rank) + '.bin')
		print 'Warning: Could not locate: ' + os.path.join(path,  dataset + '.coordinates.' + str(rank) + '.bin') + ' (old nameing sheme?)'
		return None
	
	@staticmethod
	def datafile(path, dataset, rank):
		if (os.path.exists(os.path.join(path,  dataset + '.' + str(rank) + '.bin'))):
			return os.path.join(path,  dataset + '.' + str(rank) + '.bin')
		print 'Warning: Could not locate data file for ' + dataset + '... (old nameing sheme?)'
		return None
		
	@staticmethod		
	def exists(path, dataset, rank):
		return os.path.exists(os.path.join(path,  dataset + '.coordinates.' + str(rank) + '.bin')) and os.path.exists(os.path.join(path,  dataset + '.' + str(rank) + '.bin'))

class ParallelData(SpatialData):
	new_file_format = True
	def __init__(self,path,dataset):
		super(ParallelData, self).__init__(dataset)	
		
		if ParallelData.new_file_format:
			floats_in_file   = os.path.getsize(RankData.datafile(path,dataset,rank))/4
			self.data = numpy.memmap(RankData.datafile(path,dataset,rank),dtype=numpy.float32,shape=(floats_in_file/self.nodes.shape[0],self.nodes.shape[0]))
			#extract frame time colums
			self.frames = self.data[:,0]
			# drop frame time column
			self.data = self.data[:,1:]
		else:
			#print 'Loading: ', path, dataset
			ranks = []		
			while RankData.exists(path,dataset,len(ranks)):			
				ranks.append(RankData(path, dataset, len(ranks)))	
			
					
			if len(ranks) == 0:
				raise Exception('Cannot find ' + dataset + ' files in ' + os.path.relpath(path))
			
			# combine the positions
			self.nodes = numpy.concatenate(list(rank.nodes for rank in ranks))		
			
			# get lowest frame number from rank files
			f = min(list(rank.frames.size for rank in ranks))
			self.frames = ranks[0].frames[0:f]
				
			self.data = numpy.zeros((f,self.nodes.shape[0]),dtype=numpy.float32)
			
			# merge all data in one large array
			startnode = 0
			for rank in ranks:
				self.data[:,startnode:startnode+rank.nodes.shape[0]] = rank.data[0:self.frames.size,:]
				startnode = startnode + rank.nodes.shape[0]
				
			#print 'Found', self.nodes, 'nodes and', self.frames, 'frames (t=[',self.time(0),',',self.time(self.frames-1),'])'		

class CalciumData:
	def __init__(self, path):
		super(CalciumData, self).__init__()
		print 'Mergeing files from', os.path.abspath(path)
								
		# read the data file in 10mb junks
		self.data = numpy.zeros(0, dtype = numpy.float32)
		self.frames = numpy.zeros(0, dtype = numpy.float32)
		self.nodes = numpy.zeros((0,2), dtype = numpy.float32)
				
		paths = self.resumePaths(path)
		for subpath in paths:
			subdata = ParallelData(os.path.join(path, subpath))
			
			# take nodes of first subdata and sort them lexicographically first for x, then for y 
			if self.nodes.shape[0] == 0:
				permutation = numpy.lexsort((subdata.nodes[:,1],subdata.nodes[:,0]))
				self.nodes  = subdata.nodes[permutation]				
				
			# get the insertion index
			index = self.frames.searchsorted(subdata.frames[0])
			
			# allocate additional memory
			print 'Allocating memory for %i additional frames...' % (index+subdata.frames.size - self.frames.size),
			self.frames = numpy.resize(self.frames,index+subdata.frames.size)			
			self.data = numpy.resize(self.data, (self.frames.size,self.nodes.shape[0]))
			print 'done'			
			
			# calculate permutiation
			print 'Calculating node permutation...',
			permutation = numpy.lexsort((subdata.nodes[:,1],subdata.nodes[:,0]))
			# check if we have a indentical set of nodes
			assert((self.nodes == subdata.nodes[permutation]).all())			
			print 'done'
					
			# merge data	
			self.frames[index:] = subdata.frames
			for i in range(subdata.frames.size):
				print 'Merging datasets... {progress:.2%} \r'.format(progress = i*1./(subdata.frames.size+1)),
				self.data[index + i,:] = subdata.data[i,permutation]
				#~ self.data[index:,:] = subdata.data[:,permutation]
			print '\n',
		
	def writeFiles(self, path):
		# write coordinate file
		coordfile = open(path + '/coordinates.bin','w')
		self.nodes.tofile(coordfile)
		
		# write frame file
		framefile = open(path + '/frames.bin','w')
		self.frames.tofile(framefile)
		
		# write data file
		datafile = open(os.path.join(path, 'calcium.bin'),'w')
		blocksize = 10*1024*1024/4
		blocks    = self.data.size / blocksize	
		block            = 0
		# write blocks
		while block < blocks:
			print 'Writing {path}/calcium.bin {write:>11}/{work:>11} {progress:<.2%} \r'.format(path = os.path.relpath(path), write = block*blocksize, work = self.data.size, progress = block*blocksize*1./self.data.size),
			self.data.flat[block*blocksize:(block+1)*blocksize].tofile(datafile)
			block = block + 1		
			
		# write remaining data
		self.data.flat[blocks*blocksize:].tofile(datafile)
		print 'Writing {path}/calcium.bin {write:>11}/{work:>11} {progress:<.2%}'.format(path = os.path.relpath(path), write = blocks*blocksize, work = self.data.size, progress = 1.)	
		
	# return a list of all paths
	def resumePaths(self, path):
		paths = []
		
		# check for start folder	
		if os.path.exists(os.path.join(path, 'start')):
			paths.append('start')	
			
			# check for resume folders
			resume = 1
			while os.path.exists(os.path.join(path, 'resume' + str(resume))):	
				paths.append('resume'+str(resume))
				resume += 1
				
		# check for running simulation
		if os.path.exists(os.path.join(path,'1')):
			paths.append('1')		
		
		# did we find any usefull data?
		if len(paths) == 0:
			raise Exception('Cannot find any simulation directories')
		return paths
		

