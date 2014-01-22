import os
import numpy
import csv
import ModelData
import scipy.integrate
import scipy.optimize
import math
import matplotlib.figure
import StringIO
import matplotlib.pyplot as plt

from IPython.core.pylabtools import print_figure

#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel

class Data(ModelData.Data):
	def __init__(self, path):
		
		# define the set of data types to import from the csv
		types = [('t', float), ('chid', int), ('clid', int), ('tr','>S18'), ('noch',int),('nocl',int,1), ('states',int,(10)), ('calcium',float,(10))]
		
		#~  load csv file and set up channel locations in base class
		super(Data, self).__init__(path, types)
		
		# define a set of conditions making up the state
		self.states = {}
		self.states['resting']   = 0
		self.states['active']    = 1
		self.states['open']      = 2
		self.states['inhibited'] = 3
		#~ print self.channels()

	def resting(self):
		return self.observe(lambda x: (x['states'] == self.states['resting']).sum(),desc = 'resting')
		
	def active(self):
		return self.observe(lambda x: (x['states'] == self.states['active']).sum(),desc = 'active')
		
	def inhibited(self):
		return self.observe(lambda x: (x['states'] == self.states['inhibited']).sum(),desc = 'inhibited')
				
	def open(self):
		return self.observe(lambda x: (x['states'] == self.states['open']).sum(),desc = 'open')
		
	def _repr_svg_(self):
		fig,ax=plt.subplots(figsize=(10,3))
		
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)
		ax.set_ylabel('channels')
		ax.set_xlabel('time [s]')
	
		#~ Plot the state evolution to axes object	
		self.active().plot(   ax,  c='green',lw=2)
		self.open().plot(     ax,  c='red',  lw=2)
		self.inhibited().plot(ax,  c='cyan', lw=2)
		self.resting().plot(  ax,  c='blue', lw=2)
		ax.legend(loc=2)
		
		data = print_figure(fig,'svg')
		plt.close(fig)
		return data.decode('utf-8')

	
	#~ return coordinates of channels in given state for given time
	def locations(self, time, state = 'open'):
		s = self.states[state]
		f = self.frame(time)
				
		#~ get the indices of channels in the given state
		indices = numpy.where(super(Data,self).states(time = time) == s)		
		return self.channels()[indices]
			
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
	
	
	def transitionlist(self, channel):
		states = self.data[:,[0,1+channel]]
		#~ calculate the indices of the channels transitions
		transitions = numpy.where((numpy.roll(states[:,1], 1) - states[:,1]) != 0)
		print states[transitions]
		#~ array([0, 0, 1, 0, 0, 1, 0])
	
	
	def stateintervals(self, channel, statename,frames = False):
		#~ state = self.states[statename]
		
		#~ def cond_resting(data):
			#~ return data[:,1+channel] == 0
		#~ def cond_active(data):
			#~ return data[:,1+channel] == 1
		#~ def cond_open(data):
			#~ return data[:,1+channel] == 2
		#~ def cond_inhibited(data):
			#~ return data[:,1+channel] == 3
		#~ print "requested state",statename,"for channel",channel
		def cond(data, c, s):
			#~ print "requested state",s,"for channel",c
			return data[:,1+c] == s
		
		return self.intervalls(lambda data: cond(data,c=channel,s=self.states[statename]), frames)
		
		
		#~ if statename == 'resting':
			#~ return self.intervalls(cond_resting,frames)
		#~ elif statename == 'active':
			#~ return self.intervalls(cond_active,frames)
		#~ elif statename == 'open':
			#~ return self.intervalls(cond_open,frames)
		#~ elif statename == 'inhibited':
			#~ return self.intervalls(cond_inhibited,frames)			
		#~ else:
			#~ raise "Invalid State:" + statename
		
		#times  = self.data[:,0]
		#states = self.data[:,1+channel]
				
		#~ calculate the indices of the channels transitions
		#transitions = numpy.where((numpy.roll(states, 1) - states) != 0)[0]
		
		#~ filter out all transitions different than the requested state name,
		#~ and ignore state[0] and state[-1] because they will have corrupt time.
		#~ we need to add 1 since state[0] is ignored and all indices hence are shifted to the left
		#requested   = numpy.where(states[transitions][1:-1] == state)[0] + 1
		
		#~ create array with start end end times
		#result = numpy.array([times[transitions][requested], times[numpy.roll(transitions, -1)][requested]]).transpose()
		
		#return result
		
	def statedurationdistribution(self, statename):
		result = numpy.zeros(0,dtype = numpy.float32)
		
		for channel in range(self.channelcount()):			
			intervalls = self.stateintervals(channel,statename)
			#~ append this channels durations to the return array
			result = numpy.append(result,intervalls[:,1]-intervalls[:,0])
				
		#~ return the distribution in miliseconds		
		return result * 1000		
