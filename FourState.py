import os
import numpy
import csv
import ModelData
import scipy.integrate
import scipy.optimize
import math
import matplotlib.figure
import StringIO

# Default plot for StateTimeLine
class StateTimeLine(matplotlib.figure.Figure):	
	def __init__(self, data, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True):
		super(StateTimeLine, self).__init__(figsize, dpi, facecolor, edgecolor, frameon)
		self.data  = data
			
		ax = self.add_subplot(1,1,1)
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)
		ax.axis([1, self.data.tmax(), 0, 11])
		ax.set_ylabel('channels')
		ax.set_xlabel('time [s]')
	
		#~ Plot the state evolution to axes object	
		a,= ax.plot(data.times(),data.active(),   c='green',drawstyle='steps-post',lw=2)
		o,= ax.plot(data.times(),data.open(),     c='red',  drawstyle='steps-post',lw=2)
		i,= ax.plot(data.times(),data.inhibited(),c='cyan', drawstyle='steps-post',lw=2)	
		r,= ax.plot(data.times(),data.resting(),  c='blue', drawstyle='steps-post',lw=2)

		ax.legend([r, a, o, i], ["resting","active","open","inhibited"], loc=2)
		matplotlib.backends.backend_agg.FigureCanvasAgg(self)
		


#~ At the moment the data in the csv is arranged as follows:
#~ time | channel id | cluster id | transition | open channels | open clusters | channels x state | channels x calciumlevel

#~ internally we drop some columns...
#~ time | {channel states} | {channel calciumlevels}



class Data(ModelData.Data):
	def __init__(self, path):
		#~  load csv file and set up channel locations in base class
		super(Data, self).__init__(path)
		
		self.states = {}
		self.states['resting']   = 0
		self.states['active']    = 1
		self.states['open']      = 2
		self.states['inhibited'] = 3
		#~ print self.channels()

	def resting(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['resting']).sum(1)
		
	def active(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['active']).sum(1)
		
	def inhibited(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['inhibited']).sum(1)
				
	def open(self):
		return (self.data[:,1:1+self.channelcount()] == self.states['open']).sum(1)	
	
	def closed(self):
		return self.channelcount() - self.open();
		
	def _repr_svg_(self):
		foo = StateTimeLine(self)
		imgdata = StringIO.StringIO()
		foo.savefig(imgdata, format='svg')
		imgdata.seek(0)  # rewind the data
		return imgdata.getvalue()
	
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
