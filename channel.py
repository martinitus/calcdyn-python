import math


class Channel(object):
	''' 
		infos:       a tuple containing channel id, channel location and cluster id
		transitions: a list of channel transitions, the channel will extract its own transitions from the global transition list
	'''
	def __init__(self,infos,eventdata,simulation):
		#print infos
		self.__location = infos['location']
		self.__index = infos['id']
		self.__radius = infos['radius']
		self.__cluster = infos['cluster']
		self.__eventdata = eventdata
		self.__simulation = simulation
		
		self.__state   = None;
	
	
	#~ return the event array for this channel
	def events(self):
		return self.state();
	
	# deprecated
	def state(self):
		if self.__state == None:
			i     = self.__index
			data  = self.__eventdata._data
			
			initialstate = data['states'][0][i]
			finalstate   = [data['states'][-1][i]]
			
			initialtime  = self.__simulation.tmin()
			finaltime    = self.__simulation.tmax()
			
			frames  = data[data['chid'] == i]['t']
			states  = data[data['chid'] == i]['states'][:,i]
			
			#insert initial and final state
			states  = numpy.insert(states,0, initialstate,axis = 0)
			states  = numpy.append(states,   finalstate,  axis = 0)
			
			frames  = numpy.insert(frames,0,initialtime,axis = 0)
			frames  = numpy.append(frames,  finaltime)
			
			#self.__state = TimeLine(frames,states,interpolationorder = 'zero')
			self.__state = (frames,states)
		return self.__state

	def radius(self):
		return self.__radius
		
	def location(self):
		return self.__location
		
	def index(self):
		return self.__index
		
	def cluster(self):
		return self.__cluster
		
	def open(self,t):
		tl = TimeLine(*self.state(),interpolationorder= 'zero')
		return self.__eventdata.state('open')['condition'](tl.at(t))
		
	def fluxcoefficient(self):
		if self.__simulation.config.has_option('membrane','pc'):
			return self.__simulation.config.getfloat('membrane','pc')
		else:
			ic   = self.__simulation.config.getfloat('membrane','ic')
			r    = self.__radius
			erc0 = self.__simulation.config.getfloat('er','c0')
			cyc0 = self.__simulation.config.getfloat('cytosol','c0')
			e    =  1.602176E-19
			Na   =  6.022141E23
			return ic * 1.E18 / (2 * r * r * math.pi * e * Na * (erc0 - cyc0))
		
	'''def current(self,t):
		ersnapshot = self.__simulation.domain('er')['calcium'].snapshot(t)
		cysnapshot = self.__simulation.domain('cytosol')['calcium'].snapshot(t)
		x,y=self.__location[0:2]
		fa = 96485 # faraday number
		foo = integrate.quad(lambda r: 2*math.pi*r*(ersnapshot([[x+r,y]])[0] - cysnapshot([[x+r,y]])[0]),0,self.__radius)
		return foo[0]*(1E-6*1E-24)*2*fa*self.fluxcoefficient()*1E12'''
		
	def current_integrate(self):
		x,y=self.__location[0:2]
		fa = 96485 # faraday number
		n  = len(self.__simulation.domain('cytosol')['calcium'].frames)
		foo = numpy.ndarray(n,dtype=numpy.float32)
		#print 'steps:', n
		for i in range(n):
		#	print i
			t = self.__simulation.domain('cytosol')['calcium'].frames[i]
			ersnapshot = self.__simulation.domain('er')['calcium'].snapshot(t)
			cysnapshot = self.__simulation.domain('cytosol')['calcium'].snapshot(t)
			if self.open(t):
				foo[i] = integrate.quad(lambda r: 2*math.pi*r*(ersnapshot([[x+r,y]])[0] - cysnapshot([[x+r,y]])[0]),0,self.__radius)[0]
			else:
				foo[i] = 0
		foo = foo*(1E-6*1E-24)*2*fa*self.fluxcoefficient()*1E12
		return self.__simulation.domain('cytosol')['calcium'].frames, foo
		
	def current(self):
		x,y= self.__location[0:2]
		r  = self.__radius
		fa = 96485 # faraday number
		t1,cy = self.__simulation.domain('cytosol')['calcium'].evolution(x,y)
		t2,er = self.__simulation.domain('er')['calcium'].evolution(x,y)
		
		state = self.state()
		
		mask = [self.__simulation._events.state('open')['condition'](state(ti)) for ti in t1]
		current = r*r*math.pi*2*fa*self.fluxcoefficient()*(1E-6*1E-24)*1E12 * (er-cy)*mask
		
		return t1,current
		
	def calcium(self):
		return self.__simulation.domain('cytosol')['calcium'].evolution(self.__location[0:2])
		
		
	'''def current_fast(self):	
		x,y    = self.__location[0:2]
		frames = self.__simulation.domain('cytosol')['calcium'].frames
		data   = 
		cur    = numpy.ndarray(len(frames),dtype=numpy.float32)
		for i in range(len(frames)):
			interpolator = scipy.interpolate.LinearNDInterpolator(nodes, data[7221,:],fill_value = numpy.nan)	
		datas         = self.data.swapaxes(0,1)
		
		# create time domain interpolator
		timeinterpol  = scipy.interpolate.interp1d(self.frames,datas,copy = False, fill_value = numpy.nan)
		
		# calculate scattered data values for time t
		datat         = timeinterpol(t)
		
		#create interpolator for spatial coordinates	
		return scipy.interpolate.LinearNDInterpolator(self.nodes, datat,fill_value = numpy.nan)		'''