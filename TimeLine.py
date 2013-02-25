import scipy.interpolate
import StringIO
import matplotlib.figure
import matplotlib.pyplot as plt

# provide a continuous time evolution of a discrete variable
class TimeLine(object):
	def __init__(self,frames,data,t0 = None,tend = None, desc="", ylabel="",yunit=""):
		self.frames = frames
		self.data   = data
		self.ylabel = desc
		self.yunit  = yunit
		self.desc   = desc
		
		if t0 == None:
			self.t0 = frames[0]
		else:
			self.t0 = t0
			
		if tend == None:
			self.tend = frames[-1]
		else:
			self.tend = tend			
		
		self.minframe = self.frame(self.tmin())
		self.maxframe = self.frame(self.tmax())
		assert(frames[0]<=self.t0)
		assert(self.tend     <=frames[-1])
		self.interp = scipy.interpolate.interp1d(frames,data,copy = True)
		
		
	def __call__(self,t = None, tmin = None, tmax = None):
		if t != None:
			assert(self.t0 <= t)
			assert(t  <= self.tend)
			return self.interp(t)
		else:
			assert(tmin != None)
			assert(tmax != None)
			return TimeLine(self.frames,self.data,tmin,tmax)	
	
	def tmin(self):
		return self.t0
	
	def tmax(self):
		return self.tend	
		
	def _repr_svg_(self):
		fig=self.plot()
		output = StringIO.StringIO()		
		fig.savefig(output,format='svg')
		output.seek(0)
		return output.getvalue()
		
	def plot(self, tmin = None, tmax = None, ymin = None, ymax = None):
		tmin = self.tmin() if tmin == None else tmin
		tmax = self.tmax() if tmax == None else tmax

		
		fig=plt.figure(figsize=(10,3))
		ax = fig.add_subplot(111)
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)		
		ax.set_ylabel(self.yunit + '[' + self.yunit +']')
		ax.set_xlabel('time [s]')		
		
		ax.plot(self.frames,self.data, c='black',lw=2)
		ax.set_xlim(tmin, tmax)
		
		if(ymin != None):
			ax.set_ylim(ymin,ymax)		
		
		#matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
		
		return fig
		
			#~ return the smallest frame number past the given time	
	def frame(self, time, fmin = 0, fmax = -1):
		if fmax == -1:
			fmax = self.frames.size-1
			
		#~ print 'time',time,'fmin',fmin,'fmax',fmax
		#abort recursion
		if fmax - fmin == 1:
			return fmin
			
		f = fmin + (fmax-fmin)/2
		#search in left hand side
		if self.frames[f] > time:
			#~ print 'lhs'
			return self.frame(time,fmin,f)
		#search in right hand side
		if self.frames[f] <= time:
			#~ print 'rhs'
			return self.frame(time,f,fmax)
		
		raise Exception("should never be reached...")
	
