


class Cluster(object):
	def __init__(self, index , channels, eventdata,simulation):
		self._channels = channels
		self.__eventdata = eventdata
		self.__simulation = simulation
		self.__events = None
		self.__index  = index
		
	def channels(self):
		return self._channels
	
	def channel(self,i):
		return self._channels[i]
	
	# return projection on the eventdata for this cluster
	def events(self):
		if self.__events == None:
			# the id of this cluster
			i     = self.__index

			# the indices of the channels within this cluster
			cidx  = [channel.index() for channel in self.channels()]
			#print cidx
			# the eventdata from where to extract
			data  = self.__eventdata._data

			initialstate = data['states'][0,cidx]
			#print initialstate

			finalstate   = data['states'][-1,cidx]
			#print finalstate.shape

			initialtime  = self.__simulation.tmin()
			finaltime    = self.__simulation.tmax()

			frames  = data[data['clid'] == i]['t']
			states  = data[data['clid'] == i]['states'][:,cidx]

			#print states.shape

			#insert initial and final state
			states  = numpy.insert(states,0, initialstate,axis = 0)
			states  = numpy.append(states,   [finalstate],axis = 0)

			frames  = numpy.insert(frames,0,initialtime,axis = 0)
			frames  = numpy.append(frames,  finaltime)
			
			#self.__state = TimeLine(frames,states,interpolationorder = 'zero')
			self.__events = (frames,states)
		return self.__events