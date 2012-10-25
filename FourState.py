import os
import numpy
import csv

class Data:
	def __init__(self, path):
		self.data = numpy.genfromtxt(os.path.join(path, 'transitions.csv'), delimiter=' ')
		self.states = {}
		self.states['resting']   = 0
		self.states['active']    = 1
		self.states['open']      = 2
		self.states['inhibited'] = 3
		#~ print self.data
		#~ self.data = numpy.recfromcsv(os.path.join(path, 'transitions.csv'), delimiter=' ')
	
	def tmin(self):
		return self.data[0,0];
		
	def tmax(self):
		return self.data[-1,0];
		
	def times(self):
		return self.data[:,0];
		
	def open(self):		
		return self.data[:,7];
		
	def resting(self):
		return (self.data[:,9:19] == 0).sum(1)
		
	def active(self):
		return (self.data[:,9:19] == 1).sum(1)
		
	def inhibited(self):
		return (self.data[:,9:19] == 3).sum(1)
	
	# return an array of puffs
	# puffs[p,0] puff start
	# puffs[p,1] puff end
	# puffs[p,2] puff duration
	# puffs[p,3] puff peak open count
	def puffs(self,tolerance = 0):
		active    = False
		puffs = numpy.zeros([0,4],dtype = numpy.float32)
		puff     = 0

		for line in self.data:
			# a puff just started
			if not active and line[8] >= 1:				
				# we really have a new puff
				if puff == 0 or line[0] - puffs[-1,1] > tolerance:					
					puffs.resize([puff+1,4]) #resize output array
					puff = puff + 1          # increase puff counter
					puffs[-1,0] = line[0]    # store open time
					puffs[-1,3] = 1          # set puff peak level
					print 'puff started: t=',puffs[-1,0]
				# the last puff just continued within tolerance
				else:
					pass
					
			# a puff might terminate 
			elif active and line[8] == 0:
				puffs[-1,1] = line[0]                  # store or overwrite possible puff end time
				puffs[-1,2] = line[0] - puffs[-1,1]    # store or overwrite possible puff duration
				print 'puff terminated: t=', puffs[-1,1],'duration=',puffs[-1,1]-puffs[-1,0], 'peak=',puffs[-1,2]
					
			# save present active state for next event
			active = line[8] >= 1
			
			# update puff peak level
			if active:
				puffs[-1,3]  = max(puffs[-1,3], line[7])
		
		return puffs
		
	def statedurationdistribution(self, statename):
		state = self.states[statename]
		#~ print state, statename
		
		durations  = numpy.zeros(0,dtype = numpy.float32)
		starttimes = -1*numpy.ones(10,dtype = numpy.float32)
		count = 1
		#~ print starttimes
		#~ starttimes[:] = -1
		#~ print starttimes
		#~ endtimes   = numpy.zeros(10,dtype = numpy.float32)
				
		for line in self.data:
			# the channel id of the event
			channel = line[4]			
			# channel just switched to the desired state				
			if line[channel+9] == state:
				starttimes[channel] = line[0]
				
			# channel just left state
			elif line[channel+9] != state and starttimes[channel] != -1:
				durations.resize(count);
				count = count+1
				durations[-1] = (line[0] - starttimes[channel])*1000 # return in ms
				#~ duration = 
				starttimes[channel] = -1
		#~ print durations	
		return durations

		#~ 
		#~ # read blocks
		#~ while block < blocks:
			#~ print 'Reading {path}/calcium{rank}.bin {read:>11}/{work:>11} {progress:<.0%}'.format(path=os.path.relpath(path), rank = rank, read = block*floats_per_block, work = floats_in_file, progress = block*floats_per_block*1./floats_in_file) + '\r',
			#~ self.data[block*floats_per_block:(block+1)*floats_per_block] = numpy.fromfile(datafile, dtype=numpy.float32, count=floats_per_block)
			#~ block = block + 1		
			#~ 
		#~ # read remaining floats
		#~ self.data[blocks*floats_per_block:] = numpy.fromfile(datafile, dtype=numpy.float32, count = floats_in_file - blocks*floats_per_block)
		#~ print 'Reading {path}/calcium{rank}.bin {read:>11}/{work:>11} {progress:<.0%}'.format(path=os.path.relpath(path), rank = rank, read = floats_in_file, work = floats_in_file, progress = 1.),
		#~ 
		#~ # determine number of frames
		#~ f = self.data.size / (self.nodes.shape[0]+1)
		#~ 
		#~ foo = ' '
		#~ # reshape to 2D array  containing (frames x time+nodes)		
		#~ if not self.data.size % (self.nodes.shape[0] + 1) == 0:
			#~ foo = '[corrected corrupted file]'
			#~ self.data.resize(f*(self.nodes.shape[0] + 1))
			#~ 
		#~ self.data = numpy.reshape(self.data,(f, self.nodes.shape[0] + 1))		
		#~ # extract frames from array
		#~ self.frames = self.data[:,0]
		#~ 
		#~ # drop frame time column
		#~ self.data = self.data[:,1:]	
		#~ print '   t=[{start:.8f}, {end:.8f}], {frames} frames {flag}'.format(start=float(self.frames[0]), end = float(self.frames[-1]), frames = self.frames.size, flag = foo)
