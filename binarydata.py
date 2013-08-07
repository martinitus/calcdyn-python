import os
import numpy
import scipy.spatial
import scipy.interpolate.interpnd
import scipy.integrate
import timeline

class Domain(dict):
	def __init__(self, path, name, components = ['calcium']):
		super(Domain, self).__init__()
		# the name of the domain
		self.name       = name		
		# create empty dictionary
		#self.components = {}
		# fill the dictionary with the corresponding data sets
		for component in components:
			self[component] = NewDataFormat(path,self.name + '.' + component)
			
	def components(self):
		return self.keys();
			
	#def __repr__(self):
	#	return self.components.__repr__()

class DataSet(object):
	def __init__(self,path):
		print 'Loading spatial data from', path
		self.erdomain = ERDomain(path)	
		self.cydomain = CYDomain(path)
		print self.erdomain, self.cydomain
		#print self.erdomain.frames
		#print self.cydomain.frames
		#assert((self.erdomain.frames == self.cydomain.frames).all())
		if(self.erdomain.frames.size<self.cydomain.frames.size):
			self.frames = self.erdomain.frames;		
		else:
			self.frames = self.cydomain.frames;	
		
		self.nodes = self.cydomain.nodes
		# permute data of erdomain to match cydomain node distribution
		'''permutation = numpy.zeros(self.erdomain.nodes.shape[0],dtype = numpy.int)
		i=0
		for node in self.erdomain.nodes:    
			permutation[i] = numpy.where(numpy.logical_and(self.cydomain.nodes[:,0] == node[0],self.cydomain.nodes[:,1] == node[1]))[0][0]
			i = i+1
		self.erdomain.nodes         = self.cydomain.nodes
		self.erdomain.calcium.nodes = self.cydomain.nodes
		self.erdomain.debug.nodes   = self.cydomain.nodes
		self.erdomain.calcium.data  = self.erdomain.calcium.data[:,permutation]
		self.erdomain.debug.data    = self.erdomain.debug.data[:,permutation]
	'''
	def select(self,xmin,xmax,ymin,ymax):
		self.erdomain.select(xmin,xmax,ymin,ymax)
		self.cydomain.select(xmin,xmax,ymin,ymax)
	
	def tmax(self):
		return self.frames[-1];
		
	def tmin(self):
		return self.frames[0];
		
	'''
	interpolate all contained data for x y coordinate
	'''
	def timeline(self, x , y ):
		datas = self.data.swapaxes(0,1)
		interpolator  = scipy.interpolate.LinearNDInterpolator(self.nodes, datas[:,:])				
		interpolation = interpolator([[x,y]])[0]
		return TimeLine.TimeLine(self.frames,interpolation,self.tmin(),self.tmax())
		

class SpatialData(object):
	
	def __init__(self,dataset,path=None):
		#~ TODO: make generic for 3 dimensions
		#~~self.nodes = numpy.zeros([0,3], dtype=numpy.float32)
		self.dataset = dataset
		self.nodes = numpy.zeros([0,2], dtype=numpy.float32)
		self.data  = numpy.zeros([0,0], dtype=numpy.float32)
		self.frames= numpy.zeros([0,0], dtype=numpy.float32)		
		#~ self.interpolator  = None
		self.triangulation = None
		
		if path!=None:
			print 'using new loading interface!'
			load(path, dataset, self)
	
	def select(self,xmin,xmax,ymin,ymax):
		selection = numpy.logical_and(numpy.logical_and(xmin<=self.nodes[:,0],self.nodes[:,0] <= xmax),numpy.logical_and(ymin<=self.nodes[:,1],self.nodes[:,1] <= ymax))
		self.nodes = self.nodes[selection]
		self.data  = self.data[:,selection]
        # reset triangulation
		self.triangulation = None
		
	#~ t, x, and y can all be either scalar or array types, if x and y are vectors, they both need to be the same size
	#~ if all three are vector types, the return value will be a two dimensional array containing frames in first, and spatial values in second dimension
	#~ if t is scalar, and either, x or y are arrays, the return type will be onedimensional
	def __call__(self,t = None, x = None,y = None,fill_value = numpy.nan):
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
		return result
		
	# return a callable object representing a time interpolation for the given coordinates
	def timeline(self,*args):
		if len(args) == 1:
			x=args[0][0]
			y=args[0][1]
		elif len(args)==2:
			x=args[0]
			y=args[1]
		else:
			raise Exception("dondt know what to do with:" + str(args))
		datas = self.data.swapaxes(0,1)
		interpolator  = scipy.interpolate.LinearNDInterpolator(self.nodes, datas[:,:])				
		interpolation = interpolator([[x,y]])[0]
		return timeline.TimeLine(self.frames,interpolation,ylabel=self.dataset,yunit='$\mu$M')
	
	def snapshot(self,t):
		#TODO return spatial data interpolator for time t
		return None
	
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
	def frame(self, time, fmin = 0, fmax = -1):
		if fmax == -1:
			fmax = self.frames.size -1
			
		#~ print 'time',time,'fmin',fmin,'fmax',fmax
		#abort recursion
		if fmax - fmin == 1:
			#~ if time - self.frames[fmin] < self.frames[fmax] - time :
			return fmin
			#~ else: 
				#~ return fmax
			
		f = fmin + (fmax-fmin)/2
		#search in left hand side
		if self.frames[f] > time:
			#~ print 'lhs'
			return self.frame(time,fmin,f)
		#search in right hand side
		if self.frames[f] <= time:
			#~ print 'rhs'
			return self.frame(time,f,fmax)
		
		raise Exception("should never be reached...")
	
	def tmin(self):
		return self.frames[0]
	
	def tmax(self):
		return self.frames[-1]	
	
	# return index of node closest to [x,y]
	def node(self,x,y):
		n = 0
		mn = 0
		mindist = 9E30
		for node in self.nodes:
			#~ print node
			dist = (node[0] - x)*(node[0] - x) + (node[1] -y)*(node[1] -y)
			if dist < mindist:				
				mindist = dist
				mn = n
				#~ print dist, node, mn
			n = n + 1
		return mn
		
		
	def griddata(self, time, xmin,xmax,ymin,ymax,resolution):
		xi = numpy.linspace(xmin,xmax,resolution)
		yi = numpy.linspace(ymin,ymax,resolution)
		zi = scipy.interpolate.griddata((self.nodes[:,0],self.nodes[:,1]), self(time), (xi[None,:], yi[:,None]), method='cubic')
		return zi
		
	def spatialextend(self):
		return [[self.nodes[:,0].min(),self.nodes[:,0].max()],[self.nodes[:,1].min(),self.nodes[:,1].max()]]
		
		
class NewDataFormat(SpatialData):
	def __init__(self,path,dataset):
		super(NewDataFormat, self).__init__(dataset)	
		
		filename = path + dataset + ".bin"
		# read the coordinates				
		coordfile = open(path+dataset+".coordinates.bin")					
		self.nodes = numpy.fromfile(coordfile, dtype=numpy.float32)
		print self.nodes.shape
		#print "using 3D dataset", self.nodes.size
		self.nodes = numpy.reshape(self.nodes,(self.nodes.size / 3,3))
		self.nodes = self.nodes[:,0:2]
		
		floats_in_file   = os.path.getsize(filename)/4
		
		framesize = self.nodes.shape[0] +1
		
		self.data = numpy.memmap(filename,dtype=numpy.float32,shape=(floats_in_file/framesize,framesize))
		#extract frame time colums
		self.frames = self.data[:,0]
		# drop frame time column
		self.data = self.data[:,1:]

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
		

