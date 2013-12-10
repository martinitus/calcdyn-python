import matplotlib.pyplot as plt
import numpy
import math
import scipy.interpolate
import sys, traceback
import binarydata

from matplotlib.patches import Circle

class Overview(object):	
	def __init__(self, dataset,downsample = 200):
		
		self.data   = dataset
		
		self.spatial_data = dataset.cytosol.calcium
		
		self.fig = plt.figure(self.data.path)
		self.fig.subplots_adjust(left = 0.06, bottom = 0.04, right = 0.97, top = 0.97, wspace = 0.05, hspace = 0.05)

		#~ Create for subplots with shared axes for the right side
		self.states        = self.fig.add_subplot(2,2,1)
		self.model         = self.fig.add_subplot(2,2,2,sharex=self.states)
		self.concentration = self.fig.add_subplot(2,2,3,sharex=self.states)
		self.spatial       = self.fig.add_subplot(2,2,4)
		
		# set up grids
		self.states.grid(True)
		self.concentration.grid(True)
		self.model.grid(True)
		
		# set up thick horizontal base line
		self.states.axhline(0, color='black', lw=2)
		self.concentration.axhline(0, color='black', lw=2)
		
		# set up x and y labels
		self.states.set_xlabel('t [s]')
		self.concentration.set_xlabel('t [s]')
		
		self.states.set_ylabel('channels')
		self.concentration.set_ylabel('[muM]')
		
		# set up the lines marking the currently selected frame
		self.statevax         = self.states.axvline(x=0,ymin=-0.02,ymax=10,c="red",linewidth=3,zorder=0)
		self.concentrationvax = self.concentration.axvline(x=0,ymin=-0.02,ymax=10,c="black",linewidth=3,zorder=0)

		# plot number of open channels
		self.noch = self.states.plot(self.data.events().data()['t'],self.data.events().data()['noch'])
		
		# plot model specific information
		self.data.events().model().overview(self.model,self.data.events())
		
		# plot calcium evolution at center
		self.xy    = self.spatial_data.center()
		self.ca,   = self.concentration.plot(*self.spatial_data.evolution(self.xy),label = 'cy-calcium',color = 'black')
		
		# add second axes for er concentration
		#self.concentrationer = self.concentration.twinx()
		#self.erca, = self.concentrationer.plot(*dataset.domain('er')['calcium'].evolution(self.xy),label = 'er-calcium',color = 'green')
		
		#draw legend
		self.concentration.legend()
		
		# set x plotrange of concentration and state frame
		self.concentration.set_xlim(self.data.tmin(),self.data.tmax())
		
		# plot spatial distribution
		xmin,xmax,ymin,ymax = self.spatial_data.extend()
		zi = self.spatial_data.grid(10,xmin,xmax,ymin,ymax,100)
		from matplotlib import colors 
		from matplotlib.colors import LogNorm
		self.contour = self.spatial.imshow(zi,norm=LogNorm(0.02,250,clip=True),origin='lower',extent=self.spatial_data.extend(),aspect='auto')
				
		from matplotlib.ticker import LogLocator, LogFormatter 
		l_f = LogFormatter(10, labelOnlyBase=False) 
		self.cbar = plt.colorbar(self.contour, ax = self.spatial) 
		
		for channel in self.data.channels():
			self.spatial.add_artist(Circle(channel.location(),radius = channel.radius(),fill = channel.open(0)))
		
		# connect callbacks for interactivity
		self.fig.canvas.mpl_connect('key_press_event', self.key_pressed)
		

		#~ the state statistics upper right
		
		#self.stats.grid(True)
		#self.stats.axhline(0, color='black', lw=2)
		#self.stats.axis([0,500,0,0.3])
		#self.stats.set_ylabel('state duration distribution')
		#self.stats.set_xlabel('time[ms]')

		#~ calcium levels time line for given coordinate lower right
		
		
		#self.concentration.axhline(0, color='black', lw=2)
		#self.timeline.axis([self.data.tmin(), self.data.tmax(),0,200])
		#self.concentration.set_ylabel()
		#self.concentration.set_xlabel('time [s]')
		#self.concentration.set_yscale('log')
		#self.concentration = None
		#self.concentration = None
		
		#~ the spatial calcium distribution lower left
		
		#self.spatial.axis(self.spatial_data.extend())
		#self.spatial.grid(True)
		#self.spatial.set_title('spatial calcium')
		
	#	self.state_time = 0
	#	self.x,self.y = self.spatial_data.center()
		#self.contourlevels=[0.001, 0.01, 0.1, 1, 2, 4, 8,16, 32, 64, 128, 256]
		#self.contourlevels=[400,500,550,600,650,675]
		
		#self.update_states()
		#self.update_stats()
		#self.update_timeline()
		#self.update_spatial()

		
		
		#
		#self.timeline.callbacks.connect('xlim_changed', self.timerescale)
		
		#self.spatial.callbacks.connect('xlim_changed', self.update_spatial)
		#self.spatial.callbacks.connect('ylim_changed', self.update_spatial)

	def set_time(self,t):
		'''a,= self.states.plot(self.state_data.times(),self.state_data.active(),c='green',drawstyle='steps-post',lw=2)
		o,= self.states.plot(self.state_data.times(),self.state_data.open(),c='red',drawstyle='steps-post',lw=2)
		i,= self.states.plot(self.state_data.times(),self.state_data.inhibited(),c='cyan',drawstyle='steps-post',lw=2)	
		r,= self.states.plot(self.state_data.times(),self.state_data.resting(),c='blue',drawstyle='steps-post',lw=2)'''
		ti = self.spatial_data.frames().searchsorted(t);
		t  = self.spatial_data.frames()[ti-1]
		
		self.statevax.set_xdata([t,t])
		self.concentrationvax.set_xdata([t,t])
		
		xmin,xmax = self.spatial.get_xlim()
		ymin,ymax = self.spatial.get_ylim()
		
		self.contour.set_data(self.spatial_data.grid(t,xmin,xmax,ymin,ymax,1000))
		self.contour.set_extent((xmin,xmax,ymin,ymax))
		#self.spatial.axis([xmin,xmax,ymin,ymax])
		self.time = t
		
	def set_location(self,x,y):
		self.xy = (x,y)
		
		self.ca.set_data(*self.spatial_data.evolution(x,y))
		#self.erca.set_data(*self.data.domain('er')['calcium'].evolution(x,y))
		#self.states.legend([r, a, o, i], ["resting","active","open","inhibited"], loc=2)

	def update_stats(self):
		pass
		'''rsd = self.state_data.statedurationdistribution('resting')
		asd = self.state_data.statedurationdistribution('active')
		osd = self.state_data.statedurationdistribution('open')
		isd = self.state_data.statedurationdistribution('inhibited')
		#~ foo = self.stats.hist([rsd,asd,osd,isd], 50, normed=1, facecolor='blue', alpha=0.5)
		rh = self.stats.hist(rsd, 50, normed=1, facecolor='blue', alpha=0.5)
		ah = self.stats.hist(asd, 50, normed=1, facecolor='green', alpha=0.5)
		oh = self.stats.hist(osd, 50, normed=1, facecolor='red', alpha=0.5)
		ih = self.stats.hist(isd, 50, normed=1, facecolor='cyan', alpha=0.5)
				
		proxy = [plt.Rectangle((0,0),1,1,fc = color) for color in ['blue','green','red','cyan']]
		self.stats.legend(proxy, [
			"rest =({0:3.2f} +- {1:3.2f})ms".format(numpy.mean(rsd),numpy.std(rsd)),\
			"active = ({0:3.2f} +- {1:3.2f})ms".format(numpy.mean(asd),numpy.std(asd)),\
			"open = ({0:3.2f} +- {1:3.2f})ms".format(numpy.mean(osd),numpy.std(osd)),\
			"inhibited = ({0:3.2f} +- {1:3.2f})ms".format(numpy.mean(isd),numpy.std(isd))])
		#~ 
		#~ text = ' \n active = {active:3.2f}ms\n open = {open:3.2f}ms\n inhibited = {inhibited:3.2f}ms'.\
			#~ format(rest = , active = numpy.mean(asd), open = numpy.mean(osd), inhibited = numpy.mean(isd))		
		#~ self.stats.text(200, 0.2, text, fontsize=18, ha='center', va='top')'''

	# spatial calcium plot
	def update_spatial(self, axes = None):
		
		xmin,xmax = self.spatial.get_xlim()
		ymin,ymax = self.spatial.get_ylim()
		
		self.spatial.clear()
		
		zi = self.spatial_data.grid(self.state_time,xmin,xmax,ymin,ymax,100/(xmax-xmin))

		# plot countour lines
		#~ cont2 = self.spatial.contourf(y1, z1, data, norm=plt.colors.LogNorm(lev[0],lev[len(lev)-1])) 
		from matplotlib import colors 
		from matplotlib.colors import LogNorm
		
		#~ self.contour = self.spatial.contourf(xi,yi,zi,15,norm=LogNorm())
		#self.contour = self.spatial.contourf(xi,yi,zi,self.contourlevels,norm=LogNorm())
		norm = None
		
		if self.spatial_data == self.data.domain('er')['calcium']:
			norm = LogNorm(1,700,clip=True)
		else:
			norm = LogNorm(0.02,250,clip=True)
		
		
		'''if hasattr(self,"contour"):
			self.contour.set_data(zi)
			self.contour.set_extent([xmin,xmax,ymin,ymax])
			if self.spatial_data == self.data.domain('er')['calcium']:
				self.contour.set_norm(LogNorm(700,100,clip=True))
			else:
				self.contour.set_norm(LogNorm(0.02,250,clip=True))
		else:
			self.contour = self.spatial.imshow(zi,norm=LogNorm(0.02,250,clip=True),origin='lower',extent=[xmin,xmax,ymin,ymax],aspect='auto')
			'''
			
		self.contour = self.spatial.imshow(zi,norm=norm,origin='lower',extent=[xmin,xmax,ymin,ymax],aspect='auto')
		
		self.spatial.axis([xmin,xmax,ymin,ymax])
		
		n     = self.spatial.scatter([self.x],[self.y],marker='h',c='black',s=100)
		
		for channel in self.data.channels():
			self.spatial.add_artist(Circle(channel.location(),radius = channel.radius(),fill = channel.open(self.state_time)))
						
		'''for state in self.data.events()._states.values():
			satename  = state['name']
			condition = state['condition']
			marker    = state['marker']		
			channels  = [c for c in self.data.channels() if condition(c.state().at(self.frame_time))]
			locations = numpy.array([c.location()[0:2] for c in channels])
			if locations.shape[0] > 0:
				scato = self.spatial.scatter(locations[:,0],locations[:,1],marker=marker,c='black',s=100)'''
		
		#~ plot open and closed channels as X and O 
		'''oc = self.state_data.locations(self.frame_time, 'open')
		ic = self.state_data.locations(self.frame_time, 'inhibited')
		rc = self.state_data.locations(self.frame_time, 'resting')
		ac = self.state_data.locations(self.frame_time, 'active')
		scato = self.spatial.scatter(oc[:]['x'],oc[:]['y'],marker='o',c='black',s=100)
		scati = self.spatial.scatter(ic[:]['x'],ic[:]['y'],marker='x',c='black',s=100)
		scatr = self.spatial.scatter(rc[:]['x'],rc[:]['y'],marker='s',c='black',s=100)
		scata = self.spatial.scatter(ac[:]['x'],ac[:]['y'],marker='^',c='black',s=100)	'''
		

	def key_pressed(self, event):		
		try:
			#print 'you pressed', event.key, event.xdata, event.ydata, event.inaxes
			# button 2 = middle
			if event.key == ' ':
				if event.inaxes == self.states or event.inaxes == self.concentration:
					
					self.set_time(event.xdata)
					#self.update_states()
					#self.update_spatial()
					self.fig.canvas.draw()
				elif event.inaxes == self.spatial:
					self.set_location(event.xdata,event.ydata)
					#self.update_timeline()
					#self.update_spatial()
					self.fig.canvas.draw()
			if event.key == 'n':
				if event.inaxes == self.spatial or event.inaxes == self.concentration:
					#self.select_next_spatial_data()
					#self.update_timeline()
					#self.update_spatial()
					self.fig.canvas.draw()
			if event.key == 'c':
				if event.inaxes == self.spatial:
					#self.spatial_data = self.data.domain('cytosol')['calcium']
					#self.update_timeline()
					#self.update_spatial()
					self.fig.canvas.draw()
		except Exception as detail:
			tb = traceback.format_exc()
			print tb
		finally:
			pass
		self.fig.canvas.draw()
			
	def select_next_spatial_data(self):
		if not hasattr(self,'spatial_data_cycle'):
			datasets = []
			for domain in self.data.domains():
				for comp in domain.items():
					datasets = datasets + []
		self.spatial_data = self.data.domain('cytosol')['calcium']

	def timerescale(self,axes):
		if self.time_average_ax != None:
			self.time_average_ax.remove()
			self.time_average_ax = None
		x,y   = self.spatial_data.nodes()[self.node]
		t0,t1 = self.timeline.get_xlim()
		#~ print x,y
		#~ print t0,t1
		#average = self.spatial_data.timeaverage(x,y,t0,t1)
		#self.time_average_ax = self.timeline.axhline(y=average,c="blue",linewidth=3,zorder=0)

#plot = Plot('./2/')
#~ formc = "channel #{id:1d}: {state:>9s}: ca = {avg:3.2f} T = {duration:1.4f}"
#~ forma = "average for {state:>9s} state: ca = {avg:3.2f} T = {duration:1.4f}"
#~ for state in ['resting','active','inhibited']:
	#~ channelaverage = 0
	#~ totalduration  = 0
	#~ totalstates    = 0
	#~ for channel in plot.state_data.channels():
		#~ intervalls = plot.state_data.stateintervals(channel['id'],state)		
		#~ calciumavg = plot.calcium_data.timeaverage(channel['x'],channel['y'], intervalls[:,0],intervalls[:,1])		
		#~ duration   =  numpy.sum(intervalls[:,1] - intervalls[:,0])		
		#~ print formc.format(id = channel['id'],state = state, avg = calciumavg, duration = duration / intervalls.shape[0])
		#~ 
		#~ channelaverage = channelaverage + duration * calciumavg
		#~ totalduration  = totalduration + duration
		#~ totalstates    = totalstates + intervalls.shape[0]
	#~ 
	#~ print forma.format(state = state, avg = channelaverage / totalduration, duration = totalduration / totalstates)
	#~ 
#plt.show()
