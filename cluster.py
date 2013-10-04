import numpy

class Cluster(object):
	def __init__(self, index , channels, eventdata,simulation):
		self._channels = channels
		self.__eventdata = eventdata
		self.__simulation = simulation
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
		
		# at this point we have calculated the new data, and also the puffs event data
		puffs = self.__puffs[tolerance]	
		
		for puff in puffs:
			yield self.__events[puff[0]:puff[1]]
        
        
        
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
			self.__events = numpy.recarray(shape = (eventmask.sum()),dtype = [('t', '<f8'), ('noch', '<i2'), ('states', '|i1', (len(cidx), 8))])

			# copy time chid and subspace of state column to new recarray
			self.__events['t']       = data[eventmask]['t']			
			self.__events['states']  = data[eventmask]['states'][:,cidx,:]
			
			# cache the number of open channels
			opencondition = self.__eventdata._states['open']['condition']
			self.__events['noch']    = numpy.apply_along_axis(opencondition,2,self.__events['states']).sum(axis = 1)
		return self.__events