import os
import numpy
import csv



import Utility
		
class Channel(object):
	def __init__(self,data,index):
		self.data  = data
		self.index = index
	
	#~ kwargs are forwarded to the TimeLine constructor
	def stateevolution(self,**kwargs):
		return Utility.TimeLine(self.data.data[:,0],self.data.data[:,self.index+1],**kwargs)
		
	#~ kwargs are forwarded to the TimeLine constructor	
	def calciumevolution(self,**kwargs):
		return Utility.TimeLine(self.data.data[:,0],self.data.data[:,self.index+1+self.data.channelcount()],**kwargs)
		
	def x(self):
		return self.data.channels()[self.index]['x']
	def y(self):
		return self.data.channels()[self.index]['y']
	def location(self):
		return self.data.channels()[self.index][['x','y']]


#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | {channel states} | {channel calciumlevels}

#~ internally we drop some columns...
#~ time | {channel states} | {channel calciumlevels}

class Data(object):
	def __init__(self, path):
		self.__channels = numpy.genfromtxt(os.path.join(path, 'channels.csv'),dtype=[('id', int), ('x', float), ('y', float), ('r', float)])
				
		#~ if there is only a single channel reshape the 1D array to a 1 x n 2D array
		if self.__channels.ndim == 0:
			self.__channels = self.__channels.reshape([1, self.__channels.size])
			
		self.__channelcount = self.__channels.shape[0]
		#print "found",self.__channelcount, "channels"
		
		#~ read csv file
		self.data = numpy.genfromtxt(os.path.join(path, 'transitions.csv'))
		
		#~ drop columns
		if self.data.shape[1] == 6+3+self.channelcount()*2:
			#old format with coordinates x y z
			print "Warning, still contains columns xyz, dropping!"
			self.data = self.data[:,[0] + range(9,9+2*self.channelcount())]
		elif self.data.shape[1] == 6+self.channelcount()*2:
			#new format
			self.data = self.data[:,[0] + range(6,6+2*self.channelcount())]
		else:
			raise "Unknown CSV Format!"
		
		#~ insert one line for the initial states
		self.data = numpy.insert(self.data, 0, numpy.zeros(self.data.shape[1]),axis = 0)
	
	def __repr__(self):
		return "EventData (Channels: %d, Events: %d)" % (self.__channelcount, self.data.shape[0])
		
	def _repr_svg_(self):
		return Utility.TimeLine(self.times(),self.open())
     
	#~ return the state for given time (or frame) and given channel, if neither time (frame) nor channel are give, all states for all frames are returned
	#~ def state(self, time = None, frame = None, channel = None):	
		#~ if time == None and frame != None:
			#~ pass
		#~ elif time != None and frame == None:
			#~ pass
		#~ elif time == None and frame == None:
			#~ pass
		#~ else 
			#~ raise Exception("cannot specify time and frame")
	
	def shift_coordinates(self,x,y):
		self.__channels[:]['x'] = self.__channels[:]['x'] + x
		self.__channels[:]['y'] = self.__channels[:]['y'] + y
	
	#~ return number of channels
	def channels(self):
		return self.__channels
		
	def channel(self,i):
		return Channel(self,i)
		#return self.__channels[i] #x = x.view(np.recarray) 
	
	def channelcount(self):
		return self.__channelcount
		
	#~ return number of events
	def events(self):
		return self.data.shape[0]
	
	def tmin(self):
		return self.data[0,0];
		
	def tmax(self):
		return self.data[-1,0];
	
	#~ return array with event times
	def times(self):
		return self.data[:,0];
	
	#~ return array with calcium level for given channel for all event times
	def calcium(self, channel):
		return self.data[:, 6 + self.channelcount() + channel]
		
	#~ Calculate the intervalls the given condition evealuates to true, in either real time, or frames
	#~ a valid condition should do something like this:
	#~ def condition(data):
    #~ 		return np.logical_and.reduce([data[:,1] == 0, data[:,2] == 1])
    #~ and return an array of booleans indicating wether the condition is true or false for the given frame
    #~ the first and last intervall
	def intervalls(self, condition, frames = False):
		selection = condition(self.data)
		#print 'selection',selection
		#selection = numpy.append(selection, [False])
		#~ if the condition evaluates to true for the last frame (selection[-1] == True and selection[0] == False) the following roll-xor combination will lead switch_frame[0] == True
		switch_frames = numpy.logical_xor(numpy.roll(selection, 1), selection)
		switch_frames[0] = selection[0] 
		switch_frames[-1] = False # always drop unfinished intervalls
		#print 'switchframes',switch_frames
		# detect where the the condition changes from true to false, the roll will directly mark the first and last frame where the condition is true
		start_end = switch_frames.nonzero()[0] # make the returned 0-dimensional tuple an array
		# we condition is true up to the end, we need drop the last transition to condition = true, since we cannot return a closed interval
		if start_end.shape[0] % 2 == 1:
			start_end = numpy.reshape(start_end[:-1],[start_end.size/2,2])                       # reshape the array to contain start-end pairs
		else:
			start_end = numpy.reshape(start_end,[start_end.size/2,2])                       # reshape the array to contain start-end pairs
		
		# always drop intervalls already started at t=0
		if selection[0]:
			start_end = start_end[1:]
			
		if frames:
			return start_end
		else:
			# return the intervalls where the condition is true
			interv = numpy.array([self.data[start_end[:,0],0],self.data[start_end[:,1],0]]).transpose()
			return interv
		
	#~ return the time the given channel switched into the state of given time
	def stateswitch(self,channel,time):
		f = self.frame(time)
		s = self.data[f,1+channel]
				
		#~ go back in time until we find a transition of the given channel
		while self.data[f,1+channel] == s:
			f = f - 1
			if f < 0:
				return 0		
		#~ at frame f, the channel was still in another state, hence it changed to its current state at frame f+1
		return self.data[f+1,0]
	
	#~ return state for given channel and given time
	def state(self,channel,time):
		f = self.frame(time)
		return self.data[f,1+channel]
	
	def states(self,frame = None, time = None):
		if frame != None:
			assert(time == None)			
		if time != None:
			assert(frame == None)
			frame = self.frame(time)			
		return self.data[frame,1:1+self.channelcount()]
		
	#~ return the smallest frame number past the given time	
	def frame(self, time, fmin = 0, fmax = -1):
		if fmax == -1:
			fmax = self.events()-1
			
		#~ print 'time',time,'fmin',fmin,'fmax',fmax
		#abort recursion
		if fmax - fmin == 1:
			return fmin
			
		f = fmin + (fmax-fmin)/2
		#search in left hand side
		if self.data[f,0] > time:
			#~ print 'lhs'
			return self.frame(time,fmin,f)
		#search in right hand side
		if self.data[f,0] <= time:
			#~ print 'rhs'
			return self.frame(time,f,fmax)
		
		raise Exception("should never be reached...")
	


