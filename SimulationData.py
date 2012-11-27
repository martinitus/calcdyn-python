import CalciumData
import ModelData

# hold the combination of event data, spatial data and channel data
class SimulationData(object):
	def __init__(self, path):
		raise "Not implemetnted"
	
	def __init__(self,spatial,events,channels):
		# the spatial data (CalciumData)
		self.spatial  = spatial
		# the event data (ModelData)
		self.events   = events
		# the channels 
		self.channels = channels


import matplotlib.figure

# Default plot for spatial calcium profile
class SpatialFigure(matplotlib.figure.Figure):	
	def __init__(self, data, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True, aoi = None, resolution = 100):
		super(SpatialFigure, self).__init__(figsize, dpi, facecolor, edgecolor, frameon)
		
		if aoi == None:
			xmin = min(data.channels[:]['x'])
			xmax = max(data.channels[:]['x'])
			width = xmax-xmin
			xmin = xmin - width * 0.2
			xmax = xmax + width * 0.2
		
			ymin = min(data.channels[:]['y'])
			ymax = max(data.channels[:]['y'])
			height = ymax-ymin
			ymin = ymin - height * 0.2
			ymax = ymax + height * 0.2
		else:
			xmin,xmax,ymin,ymax = aoi[:]
		
		xi = numpy.linspace(xmin,xmax,resolution)
		yi = numpy.linspace(ymin,ymax,resolution)
		
		# grid the data.	
		frame = data.spatial.frame(self.frame_time)
		
		data = self.calcium_data.data[frame,:]

		zi = scipy.interpolate.griddata((self.calcium_data.nodes[:,0],self.calcium_data.nodes[:,1]), data, (xi[None,:], yi[:,None]), method='cubic')

		# plot countour lines
		#~ cont2 = self.spatial.contourf(y1, z1, data, norm=plt.colors.LogNorm(lev[0],lev[len(lev)-1])) 
		from matplotlib import colors 
		from matplotlib.colors import LogNorm
		
		#~ self.contour = self.spatial.contourf(xi,yi,zi,15,norm=LogNorm())
		self.contour = self.spatial.contourf(xi,yi,zi,self.contourlevels,norm=LogNorm())
				
		#~ plot open and closed channels as X and O 
		oc = self.state_data.locations(self.frame_time, 'open')
		ic = self.state_data.locations(self.frame_time, 'inhibited')
		rc = self.state_data.locations(self.frame_time, 'resting')
		ac = self.state_data.locations(self.frame_time, 'active')
		scato = self.spatial.scatter(oc[:]['x'],oc[:]['y'],marker='o',c='black',s=100)
		scati = self.spatial.scatter(ic[:]['x'],ic[:]['y'],marker='x',c='black',s=100)
		scatr = self.spatial.scatter(rc[:]['x'],rc[:]['y'],marker='s',c='black',s=100)
		scata = self.spatial.scatter(ac[:]['x'],ac[:]['y'],marker='^',c='black',s=100)	
		n     = self.spatial.scatter([self.calcium_data.nodes[self.node,0]],[self.calcium_data.nodes[self.node,1]],marker='h',c='black',s=100)	
