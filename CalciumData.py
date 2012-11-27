import os
import numpy
import scipy.spatial
import scipy.interpolate.interpnd
import scipy.integrate

# provide a continuous time evolution of a discrete variable
class TimeLine(object):
	def __init__(self,frames,data,t0,tend):
		self.frames = frames
		self.data   = data
		self.t0     = t0
		self.tend   = tend
		assert(frames[0]<=t0)
		assert(tend     <=frames[-1])
		self.interp = scipy.interpolate.interp1d(frames,data,copy = True)
		
	def __call__(self,t):
		assert(self.t0 <= t)
		assert(t  <= self.tend)
		return self.interp(t)
	
	def tmin(self):
		return self.t0
	
	def tmax(self):
		return self.tend	

class SpatialData(object):
	
	# TODO: check if this class is used somewhere...
	class My2DInterpolator(scipy.interpolate.interpnd.LinearNDInterpolator): 
		def __init__(self, tri, values, fill_value = numpy.nan, tol=1e-6, maxiter=400): 
			scipy.interpolate.interpnd.NDInterpolatorBase.__init__(self, tri.points, values, ndim=2, fill_value=fill_value) 
			self.tri = tri			
			self.grad = scipy.interpolate.interpnd.estimate_gradients_2d_global(self.tri, self.values, tol=tol, maxiter=maxiter) 
	
	def __init__(self):
		self.nodes = numpy.zeros([0,2], dtype=numpy.float32)
		self.data  = numpy.zeros([0,0], dtype=numpy.float32)
		self.frames= numpy.zeros([0,0], dtype=numpy.float32)		
		#~ self.interpolator  = None
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
		
		
	def spatialinterpolation(self,x,y,firstframe,lastframe):
		warnings.warn("spatialinterpolation is derprecated, use __call__ instead!", DeprecationWarning )
		datas = self.data.swapaxes(0,1)
		interpolator  = scipy.interpolate.LinearNDInterpolator(self.nodes,datas[:,firstframe:lastframe+1],fill_value = fill_value)
		
	# return a callable object representing a time interpolation for the given coordinates
	def timeline(self, x,y):
		datas = self.data.swapaxes(0,1)
		interpolator  = scipy.interpolate.LinearNDInterpolator(self.nodes, datas[:,:])		
		
		interpolation = interpolator([[x,y]])[0]
		return TimeLine(self.frames,interpolation,self.tmin(),self.tmax())
		#print "interpolation:",interpolation
		
		#return scipy.interpolate.interp1d(self.frames,interpolation, copy = True)
	
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
	
	#~ t0 and t1 can either be scalars or arrays given start end endtimes of separated integration intervals
	def timeintegration(self,x,y,t0,t1):		
		if t0.__class__.__name__ != 'ndarray':
			assert(t0.__class__.__name__ != 'ndarray')
			t0 = numpy.array([t0])
			t1 = numpy.array([t1])		
		
		firstframe = self.frame(t0[0])
		lastframe  = self.frame(t1[-1])+1
				
		#~ interpolate all required frame values for the given coordinates
		datas = self.data.swapaxes(0,1)
		spatial_interpolation  = scipy.interpolate.LinearNDInterpolator(self.nodes, datas[:,firstframe:lastframe+1])([[x,y]])[0]
		
		#~ print "requested range: [",t0[0],',',t1[-1],']'
		#~ print "interpolatorrange: [",self.frames[firstframe],',',self.frames[lastframe],']'
		#~ create time domain interpolator object
		interpolator = scipy.interpolate.interp1d(self.frames[firstframe:lastframe+1],spatial_interpolation,copy = False)
		
		result = 0
		#~ since we cannot set custom ranges for scipy.integrate.trapz, we need to excplicitly handle the limits
		for i in range(t0.shape[0]):
			f0     = self.frame(t0[i]) # the frame before the start of the present intervall
			f1     = self.frame(t1[i]) # the frame before the end   of the present intervall
			df     = (self.frames[f0+1] - t0[i]          ) # the delta for front
			de     = (t1[i]             - self.frames[f1]) # the delta for end
		
			front  = df * interpolator(t0[i]           + df / 2.)
			end    = de * interpolator(self.frames[f1] + de / 2.)
			middle = scipy.integrate.trapz(spatial_interpolation[f0-firstframe+1:f1-firstframe+1],self.frames[f0+1:f1+1])
		
			result = result + front + end + middle
		
		return result
		
	def timeaverage(self,x,y,t0,t1):
		#~ print "hasattr(t0, '__class__')",hasattr(t0, '__class__')		
		#~ if hasattr(t0, '__class__'):
			#~ print "t0.__class__.__name__", t0.__class__.__name__
		#~ print "hasattr(t0, '__get__')",hasattr(t0, '__get__')		
		#~ print "hasattr(t0, 'shape')",hasattr(t0, 'shape')		
		#~ print "hasattr(t0, 'size')",hasattr(t0, 'size')		
		
		if t0.__class__.__name__ != 'ndarray':
			assert(t0.__class__.__name__ != 'ndarray')
			t0 = numpy.array([t0])
			t1 = numpy.array([t1])		
		
		duration = numpy.sum(t1-t0)
		#~ print "duration:",duration
		
		return self.timeintegration(x,y,t0,t1) / duration
	

class RankData(SpatialData):
	def __init__(self, path, rank = '',verbose = False):
		super(RankData, self).__init__()
		# read the coordinates				
		coordfile = open(os.path.join(path, 'coordinates' + rank + '.bin'))		
		self.nodes = numpy.fromfile(coordfile, dtype=numpy.float32)
		self.nodes = numpy.reshape(self.nodes,(self.nodes.size / 2,2))					
				
		# read the data file in 10mb junks
		datafile = open(os.path.join(path, 'calcium' + rank + '.bin'))
		floats_per_block = 10*1024*1024/4
		floats_in_file   = os.path.getsize(os.path.join(path, 'calcium' + rank + '.bin'))/4
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

class ParallelData(SpatialData):
	def __init__(self, path):
		super(ParallelData, self).__init__()	
		#print 'Trying to merge files from', path
		ranks = self.loadRankFiles(path)
		
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
		
	def loadRankFiles(self,path):
		ranks = []
		rank = 0
		while os.path.exists(os.path.join(path, 'calcium_rank_' + str(rank) + '.bin')) and os.path.exists(os.path.join(path, 'coordinates_rank_' + str(rank) + '.bin')):	
			ranks.append(RankData(path, '_rank_' + str(rank)))	
			rank += 1
			
		if len(ranks) == 0:
			# check for sequential output file
			if os.path.exists(os.path.join(path, 'calcium.bin')) and os.path.exists(os.path.join(path, 'coordinates.bin')):	
				ranks.append(RankData(path))
				
		if len(ranks) == 0:
			raise Exception('Cannot find any data files in ' + os.path.relpath(path))
		return ranks
		

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
		

