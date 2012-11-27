import os
import numpy
import csv
import ModelData

#~ At the moment the data in the csv is arranged as follows:
#~ time | channel location (x y z) | channel id | cluster id | transition | open channels | open clusters | channels x state               | channels x calciumlevel
#~ old format :-/
#~ 0.25   2000 2000 0 					0 			0 			open 			1 				1 				1(opencount) 0(closecount) 
#~ 1 	  2000 2000 0 					0 			0 			close 			0 				0 				0 1 
class Data(ModelData.Data):
	def __init__(self, path):
		#~  load csv file and set up channel locations in base class
		super(Data, self).__init__(path)
		
		#~ self.data = numpy.recfromcsv(os.path.join(path, 'transitions.csv'), delimiter=' ')
		
		#~ check if everything went right
		#~ print self.data
		
		#~ define states for this model
		self.states = {}		
		self.states['closed'] = 0
		self.states['open']   = 1
		
	#~ def calcium(self, channel):
		#~ return self.data[:,9+self.channels.shape[0]+channel]
	
	def open(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['open']).sum(1)	
	
	def closed(self):
		return self.channelcount() - self.open();
	
	#~ return coordinates of channels in given state for given time
	def locations(self, time, state = 'open'):
		s = self.states[state]
		f = self.frame(time)
				
		#~ get the indices of channels in the given state
		indices = numpy.where(super(Data,self).states(time = time) == s)		
		return self.channels()[indices]
