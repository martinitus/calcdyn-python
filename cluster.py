import numpy
import dyk
from eventcollection import EventCollection

class Puff(object):
	def __init__(self, cluster, firstframe, lastframe):
		self.__cluster = cluster
		self.__start   = firstframe
		self.__end     = lastframe
		self.__events  = cluster.events()[firstframe:lastframe]
		
		self.__model = cluster.model()
		#print "created puff", self.start(), self.end(), self.duration(), self.peak()
		for attribute in self.__model.collective_event_attributes():
			self.__dict__[attribute.name] = attribute.value(self.__events)
	
	
	def trajectory(self, dt, start = None, stop = None):
		start = self.start() if start == None else start;
		stop  = self.end()   if stop  == None else stop;
		for t in numpy.arange(start = start,stop = stop,step = dt):
			yield self.__events[self.__events['t'].searchsorted(t)]
	
	def peak(self):
		return self.__events['noch'].max()
		
	def events(self):
		return self.__events
		
	def start(self):
		return self.__events['t'][0]
		
	def end(self):
		return self.__events['t'][-1]
		
	def duration(self):
		return self.end()-self.start()
		
	def model(self):
		return self.__model
		
	# return the time integral over the number of open channels
	def accumulated(self):
		no = self.__events['noch']
		t =  self.__events['t']
		dts = (t[1:]-t[:-1])
		return numpy.dot(no[:-1],dts)
	
	def average(self):
		return self.accumulated()/self.duration()
		

class Cluster(object):
	def __init__(self, index , channels, eventdata):
		self._channels = channels
		self.__eventdata = eventdata
		self.__model = eventdata.model()
		#self.__simulation = simulation
		self.__events = None
		self.__index  = index
		self.__puffs  = {}
	
	# deprecated
	def channels(self):
		return self._channels
	
	# make the cluster iterable over the channels
	def __iter__(self):
		for c in self._channels:
			yield c
	
	def model(self):
		return self.__model

	# return iterable for iteration over puffs
	def puffs(self, tolerance = 0.005):
		
		#if we do not have the data available for the asked tolerance, calculate it and store it in dictionary
		if not self.__puffs.has_key(tolerance):
			data = self.events()
			
			wasopen = data[0]['noch'] > 0
			self.__puffs[tolerance] = []
			
			puffs = self.__puffs[tolerance]
			
			active    = False
			lastclose = -2*tolerance

			for i in range(len(data)):
				isopen = data[i]['noch'] > 0
					
				#first channel opened
				if not wasopen and isopen:
					withintolerance = data[i]['t']-lastclose <= tolerance
				
					if not withintolerance:
						self.__puffs[tolerance] = self.__puffs[tolerance] + [[i,-1]]
			
				#last channel closed
				if not isopen and wasopen:
					lastclose = data[i]['t']    
					self.__puffs[tolerance][-1][1] = i
				if isopen:
					self.__puffs[tolerance][-1][1] = i
				
				wasopen = isopen
				
			# finished calculation, transform raw data to list of puff objects, NOTE: add +1 to include the last close event
			self.__puffs[tolerance] = [Puff(self,puff[0],puff[1]+1) for puff in self.__puffs[tolerance]]
		
		
		# at this point we have calculated the new data, and also the puffs event data
		#~ puffs = self.__puffs[tolerance]	
		
		#~ for puff in puffs:
		#	#~ yield Puff(self,puff[0],puff[1])
		return self.__puffs[tolerance]
        
        
        
	def channel(self,i):
		return self._channels[i]
	
	# return projection on the eventdata for this cluster
	def events(self):
		if self.__events == None:
			# the indices of the channels within this cluster
			cidx  = [channel.index() for channel in self]
			#print cidx
			# the eventdata from where to extract
			data  = self.__eventdata._data

			#select all events for this cluster
			eventmask = data['clid'] == self.__index
			#select first and last frame
			eventmask[0]  = True
			eventmask[-1] = True

			#create recarray that stores the events of this cluster
			self.__events = numpy.recarray(shape = (eventmask.sum()),dtype = [('t', '<f8'), ('noch', '<i2'), ('chid', int), ('states', '|i1', (len(cidx), 8))])

			# copy time chid and subspace of state column to new recarray
			self.__events['t']       = data[eventmask]['t']
			self.__events['chid']    = data[eventmask]['chid']
			self.__events['states']  = data[eventmask]['states'][:,cidx,:]
			
			# cache the number of open channels			
			model =  self.__eventdata.model()
			self.__events['noch'] = model.noch(self.__events)
			
		return self.__events
		
	def couplingcoefficient(self,available):
		ec = EventCollection(self.puffs(tolerance = 0.0)).filter(lambda x: x.start()>100 and dyk.available(x.events()[0]) == available)
		po = 1.*(ec['peak'] > 1).sum()/len(ec.events)
		c  = 1.-(1.-po)**(1./(available-1))
		#print "#available = ",available,":" , (ec['peak'] > 1).sum(),'/',len(ec.events),'=',round(po,2),'=> c=',round(c,2)
		return c