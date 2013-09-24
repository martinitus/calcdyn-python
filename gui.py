import matplotlib.pyplot as plt
import numpy
import math
import scipy.interpolate
import sys, traceback
import binarydata

from matplotlib.patches import Circle

class Overview(object):	
	def __init__(self, dataset,downsample = 200):
		
		self.downsample = downsample
		self.data   = dataset
		
		self.spatial_data = dataset.domain('cytosol')['calcium']
		
		self.fig = plt.figure(self.data.path)
		self.fig.subplots_adjust(left = 0.06, bottom = 0.04, right = 0.97, top = 0.97, wspace = 0.05, hspace = 0.05)

		#~ the timeline upper left
		self.states = self.fig.add_subplot(2,2,1)
		self.states.grid(True)
		self.states.axhline(0, color='black', lw=2)
		self.states.set_xlim([self.data.tmin(), self.data.tmax()])
		self.states.set_ylabel('channels')
		self.states.set_xlabel('time [s]')
		self.statevax = self.states.axvline(x=0,ymin=-0.02,ymax=10,c="black",linewidth=3,zorder=0)
		self.framevax = self.states.axvline(x=0,ymin=-0.02,ymax=10,c="black",linewidth=3,zorder=0)

		#~ the state statistics upper right
		self.stats = self.fig.add_subplot(2,2,2)
		self.stats.grid(True)
		#self.stats.axhline(0, color='black', lw=2)
		#self.stats.axis([0,500,0,0.3])
		#self.stats.set_ylabel('state duration distribution')
		#self.stats.set_xlabel('time[ms]')

		#~ calcium levels time line for given coordinate lower right
		self.timeline = self.fig.add_subplot(2,2,3,sharex=self.states)
		self.timeline.grid(True)
		self.timeline.axhline(0, color='black', lw=2)
		self.timeline.axis([self.data.tmin(), self.data.tmax(),0,200])
		self.timeline.set_ylabel('[muM]')
		self.timeline.set_xlabel('time [s]')
		self.timeline.set_yscale('log')
		self.timelineplot = None
		self.time_average_ax = None
		
		#~ the spatial calcium distribution lower left
		self.spatial = self.fig.add_subplot(2,2,4)
		self.spatial.axis(self.spatial_data.spatialextend())
		self.spatial.grid(True)
		self.spatial.set_title('spatial calcium')
		
		self.state_time = 0
		self.frame_time = 0
		self.node = self.spatial_data.node(self.spatial_data.center())
		#self.contourlevels=[0.001, 0.01, 0.1, 1, 2, 4, 8,16, 32, 64, 128, 256]
		#self.contourlevels=[400,500,550,600,650,675]
		
		self.update_states()
		self.update_stats()
		self.update_timeline()
		self.update_spatial()
		
#~ self.fig.colorbar(self.contour, ax = self.spatial)
		from matplotlib.ticker import LogLocator, LogFormatter 
		l_f = LogFormatter(10, labelOnlyBase=False) 
		self.cbar = plt.colorbar(self.contour, ax = self.spatial) 
		
		
		self.fig.canvas.mpl_connect('key_press_event', self.key_pressed)
		#self.timeline.callbacks.connect('xlim_changed', self.timerescale)
		
		#self.spatial.callbacks.connect('xlim_changed', self.update_spatial)
		#self.spatial.callbacks.connect('ylim_changed', self.update_spatial)

	def update_states(self):
		#self.statevax.remove()		
		#self.framevax.remove()
		self.states.clear()
		xmin,xmax = self.states.get_xlim()
		
		'''a,= self.states.plot(self.state_data.times(),self.state_data.active(),c='green',drawstyle='steps-post',lw=2)
		o,= self.states.plot(self.state_data.times(),self.state_data.open(),c='red',drawstyle='steps-post',lw=2)
		i,= self.states.plot(self.state_data.times(),self.state_data.inhibited(),c='cyan',drawstyle='steps-post',lw=2)	
		r,= self.states.plot(self.state_data.times(),self.state_data.resting(),c='blue',drawstyle='steps-post',lw=2)'''
		self.data.events().openchannels().plot(self.states)
		self.statevax = self.states.axvline(x=self.state_time,ymin=-0.02,ymax=10,c="yellow",linewidth=3,zorder=0)
		self.framevax = self.states.axvline(x=self.frame_time,ymin=-0.02,ymax=10,c="black",linewidth=3,zorder=0)
		self.states.grid(True)
		self.states.set_xlim([xmin,xmax])
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

	def update_timeline(self):
		if self.timelineplot != None:
			self.timelineplot.pop(0).remove()
			
		#self.timeline.clear()
		self.timelineplot = self.timeline.plot(self.spatial_data.frames[::self.downsample],self.spatial_data.data[::self.downsample,self.node],c='black',lw=2)

	# spatial calcium plot
	def update_spatial(self, axes = None):
		
		xmin,xmax = self.spatial.get_xlim()
		ymin,ymax = self.spatial.get_ylim()
		
		self.spatial.clear()
		#self.spatial.axis([xmin,xmax,ymin,ymax])
		
		xi = numpy.linspace(xmin,xmax,100)
		yi = numpy.linspace(ymin,ymax,100)
		
		# grid the data.	
		frame = self.spatial_data.frame(self.frame_time)
		#~ data = numpy.log(self.calcium_data.data[frame,:])
		data = self.spatial_data.data[frame,:]

		zi = scipy.interpolate.griddata((self.spatial_data.nodes[:,0],self.spatial_data.nodes[:,1]), data, (xi[None,:], yi[:,None]), method='cubic')

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
		
		n     = self.spatial.scatter([self.spatial_data.nodes[self.node,0]],[self.spatial_data.nodes[self.node,1]],marker='h',c='black',s=100)
		
		for channel in self.data.channels():
			self.spatial.add_artist(Circle(channel.location(),radius = channel.radius(),fill = channel.open(self.frame_time)))
						
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
				if event.inaxes == self.states or event.inaxes == self.timeline:
					self.state_time = self.data.events()._data[self.data.events().frame(event.xdata)]['t']
					self.frame_time = self.spatial_data.frames[self.spatial_data.frame(event.xdata)]				
					self.update_states()
					self.update_spatial()
					self.fig.canvas.draw()
				elif event.inaxes == self.spatial:
					self.node = self.spatial_data.node(event.xdata,event.ydata)
					self.update_timeline()
					self.update_spatial()
					self.fig.canvas.draw()
			if event.key == 'n':
				if event.inaxes == self.spatial or event.inaxes == self.timeline:
					self.select_next_spatial_data()
					self.update_timeline()
					self.update_spatial()
					self.fig.canvas.draw()
			if event.key == 'c':
				if event.inaxes == self.spatial:
					self.spatial_data = self.data.domain('cytosol')['calcium']
					self.update_timeline()
					self.update_spatial()
					self.fig.canvas.draw()
					
		except Exception as detail:
			tb = traceback.format_exc()
			print tb
		finally:
			pass
			
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
		x,y   = self.spatial_data.nodes[self.node]
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
