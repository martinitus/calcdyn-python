class FakeSpatialData(object):
	
	def __init__(self,realdata):
		self._realdata = realdata
		
	# return a callable object representing a time interpolation for the given coordinates
	def timeline(self,*args):
		return lambda x: 700
	
	# retrief a interpolator object for 2D slice at given time
	def snapshot(self,t):
		#create interpolator for spatial coordinates	
		return (lambda x: [700])
	
	#return the frame number closest to the given time
	def frame(self, time, fmin = 0, fmax = -1):
		return self._realdata.frame(time,fmin,fmax)
	
	def tmin(self):
		return self._realdata.tmin()
	
	def tmax(self):
		return self._realdata.tmax()
	
	# return index of node closest to [x,y]
	def node(self,*args):
		return self._realdata.nodes(args)		
		
	def griddata(self, time, xmin,xmax,ymin,ymax,resolution):
		xi = numpy.linspace(xmin,xmax,resolution)
		yi = numpy.linspace(ymin,ymax,resolution)
		zi = scipy.interpolate.griddata((self.nodes[:,0],self.nodes[:,1]), self(time), (xi[None,:], yi[:,None]), method='cubic')
		return zi
		
	def spatialextend(self):
		return self._realdata.spatialextend()
	
	def center(self):
		self._realdata.center()
		
class FakeDomain(dict):
	def __init__(self, name, simulation, components = ['calcium']):
		super(FakeDomain, self).__init__()
		# the name of the domain
		self.name       = name	
		self.simulation = simulation	
		# fill the dictionary with the corresponding data sets
		for component in components:
			self[component] = FakeSpatialData(simulation.domain('cytosol')['calcium'])
			
	def components(self):
		return self.keys();		