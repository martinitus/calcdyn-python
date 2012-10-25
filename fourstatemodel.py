import os
import numpy
import csv

class FourStateData:
	def __init__(self, path):
		self.data = numpy.genfromtxt(os.path.join(path, 'transitions.csv'), delimiter=' ')
		
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
