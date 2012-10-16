import os
import numpy

class RankData:
	def __init__(self, path, rank = ''):
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
			print 'Reading {path}/calcium{rank}.bin {read:>11}/{work:>11} {progress:<.0%}'.format(path=os.path.relpath(path), rank = rank, read = block*floats_per_block, work = floats_in_file, progress = block*floats_per_block*1./floats_in_file) + '\r',
			self.data[block*floats_per_block:(block+1)*floats_per_block] = numpy.fromfile(datafile, dtype=numpy.float32, count=floats_per_block)
			block = block + 1		
			
		# read remaining floats
		self.data[blocks*floats_per_block:] = numpy.fromfile(datafile, dtype=numpy.float32, count = floats_in_file - blocks*floats_per_block)
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
		print '   t=[{start:.8f}, {end:.8f}], {frames} frames {flag}'.format(start=float(self.frames[0]), end = float(self.frames[-1]), frames = self.frames.size, flag = foo)

class ParallelData:
	def __init__(self, path):		
		#print 'Trying to merge files from', path
		ranks = self.loadRankFiles(path)
		
		# combine the positions
		self.nodes = numpy.concatenate(list(rank.nodes for rank in ranks))		
		
		# get lowest frame number from rank files
		f = min(list(rank.frames.size for rank in ranks))
		self.frames = ranks[0].frames[0:f]
			
		self.data = numpy.zeros((f,self.nodes.size),dtype=numpy.float32)
		
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

class CalciumData:
	def __init__(self, path):
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
		

