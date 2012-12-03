# provide a continuous time evolution of a discrete variable
class TimeLine(object):
	def __init__(self,frames,data,t0,tend):
		self.frames = frames
		self.data   = data
		self.t0     = t0
		self.tend   = tend
		assert(frames[0]<=t0)
		assert(tend     <=frames[-1])
		self.interp = scipy.interpolate.interp1d(frames,data,copy = True)
		
	def __call__(self,t):
		assert(self.t0 <= t)
		assert(t  <= self.tend)
		return self.interp(t)
	
	def tmin(self):
		return self.t0
	
	def tmax(self):
		return self.tend	
		
	def _repr_svg_(self):
		plt.ioff() # turn off interactive mode
		fig=plt.figure(figsize=(10,3))
		ax = fig.add_subplot(111)
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)		
		ax.set_ylabel('ca [microMolar]')
		ax.set_xlabel('time [s]')		
		ax.plot(self.frames,self.data, c='black',lw=2)
		ax.axis([self.tmin(), self.tmax(), 0, 11])
		output = StringIO.StringIO()
		fig.savefig(output,format='svg')
		output.seek(0)
		plt.ion() # turn on interactive mode
		return output.getvalue()
