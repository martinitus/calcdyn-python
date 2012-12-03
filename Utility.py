import scipy.interpolate.interpnd
import matplotlib.figure
import StringIO
import matplotlib.backends.backend_agg

# provide a continuous time evolution of a discrete variable
class TimeLine(object):
	def __init__(self,frames,data, tmin = None, tmax = None, ylabel='y', yscale = 'linear', ymin = None, ymax = None):
		self.frames = frames
		self.data   = data		
		if(tmin != None):		
			assert(frames[0]<=tmin)
			self.__tmin = tmin
		else:
			self.__tmin = frames[0]
		if(tmax != None):			
			assert(frames[-1]>=tmax)		
			self.__tmax = tmax
		else:
			self.__tmax = frames[-1]
		self.interp = scipy.interpolate.interp1d(frames,data,copy = False)
		self.ylabel = ylabel
		self.yscale  = yscale
		self.__ymin = ymin
		self.__ymax = ymax
		
	def __call__(self,t):
		assert(self.frames[0] <= t)
		assert(t  <= self.frames[-1])
		return self.interp(t)
	
	def tmin(self):
		return self.__tmin
	
	def tmax(self):
		return self.__tmax	
	
	def tlimits(self, tmin = None, tmax = None):
		if tmin!=None:
			assert(self.frames[0]<=tmin)
			self.__tmin=tmin
		if tmax!=None:
			assert(tmax<=self.frames[-1])
			self.__tmax=tmax
		return self.__tmin, self.__tmax
	
	def ylabel(self,string):
		self.ylabel = string
	
	def yscale(self,scale):
		assert(scale == 'log' or scale == 'linear')
		self.yscale = scale
	
	def ylimits(ymin=None,ymax=None):
		self.__ymin = ymin
		self.__ymax = ymax	
		
	def _repr_svg_(self):
		fig=matplotlib.figure.Figure(figsize=(10,3))
		ax = fig.add_subplot(111)
		ax.grid(True)
		ax.axhline(0, color='black', lw=2)		
		ax.set_ylabel(self.ylabel)
		ax.set_xlabel('time [s]')
		ax.plot(self.frames,self.data, c='black',lw=2)
		ax.set_yscale(self.yscale)
		ax.set_xlim(self.__tmin, self.__tmax)
		ax.set_ylim(self.__ymin, self.__ymax)		
		matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
		output = StringIO.StringIO()
		fig.savefig(output,format='svg')
		output.seek(0)
		return output.getvalue()
		
	#~ t0 and t1 mark start end endtimes of integration interval
	def integrate(self,t0,t1):	
		assert(t0 >= self.frames[0])
		assert(t1 <= self.frames[-1])
				
		firstframe = self.frame(t0)
		lastframe  = self.frame(t1)+1
	
		f0     = self.frame(t0) # the frame before the start of the intervall
		f1     = self.frame(t1) # the frame before the end   of the intervall
		df     = (self.frames[f0+1] - t0             ) # the delta for front
		de     = (t1                - self.frames[f1]) # the delta for end
	
		front  = df * self(t0              + df / 2.)
		end    = de * self(self.frames[f1] + de / 2.)
		middle = scipy.integrate.trapz(self.data[f0+1:f1+1],self.frames[f0+1:f1+1])
	
		return front + end + middle	
		
	#~ return the smallest frame number past the given time	
	def frame(self, time, fmin = 0, fmax = -1):
		if fmax == -1:
			fmax = self.frames.shape[0]-1
			
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
	
	def average(self,t0,t1):		
		return self.integrate(t0,t1) / (t1-t0)
