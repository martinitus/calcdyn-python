#import StochasticModel
import numpy
import scipy
import math
#import timeline
#import deterministic

def types_binary(channels):
    return [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32), ('noch',numpy.int32),('nocl',numpy.int32), ('states',numpy.bool,(channels))] #''
    

def loadtxt(filename,channels):
    tmp = numpy.loadtxt(filename, dtype = [('t', numpy.double), ('chid', numpy.int32), ('clid', numpy.int32),
                                           ('tr','>S8'), ('noch',numpy.int32),('nocl',numpy.int32),
                                           ('states','>S6',(channels))])
                                   
    return numpy.fromiter(itertools.izip(tmp['t'],tmp['chid'],tmp['clid'],tmp['noch'],tmp['nocl'],[[s == 'open' for s in f] for f in tmp['states']]),
                                     dtype = types_binary(channels))

def state_type():
    return ('>S6')

def open(state):
    return state == 'open'
#states = deterministic.states

'''RyRModel::usage = "
Definition of propensities a[t], CDF P[t] and probability density \[Rho][t]:
P(\[Tau]) =  1-Exp[-\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Tau]\)]\)a[s]\[DifferentialD]s] 
\[Rho](\[Tau]) =  a[\[Tau]] Exp[-\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Tau]\)]\)a[s]\[DifferentialD]s]"
'''
alpha = 0.15;
km    = 720.;
rhop  = 30000.;
rhom  = 480.;

'''
open::usage     = "calculate open rate in [1/s]"
close::usage    = "calculate close rate  in [1/s]"
a::usage        = "calculate open propensity = n*openrate"
P::usage        = "calculate CDF"
\[Rho]::usage        = "calculate probability density"
T::usage        = "calculate expectation value of first passage time"
e::usage        = "defines the stickyness energy, set to Zero for no stickyness"
Popen::usage    = "calculate open probability for given cytosolic and luminal calcium concentration"
coupling::usage = "the coupling term"

Kdopen::usage     = "return the \!\(\*SubscriptBox[\(K\), \(d\)]\) for the open rate given luminal calcium concentration for given \[Alpha] and \!\(\*SubscriptBox[\(K\), \(max\)]\)"
Kd::usage         = "return the \!\(\*SubscriptBox[\(K\), \(d\)]\) for the open probability"
'''


def Kd(er):
    return math.pow(rhom/(rhop+rhom), 1./4.)*(alpha*(km-er));

def kclose(cy,er):
    return rhom;
    
def kopen(cy,er):
    return rhop * (cy**4)/(cy**4+(alpha*(km-er))**4);
    
def tclose(cy,er):
    return 1./kopen(cy,er);
    
def topen(cy,er):
    return 1./kclose(cy,er);

def popen(cy,er):
    return topen(cy,er)/(topen(cy,er)+tclose(cy,er));
    
def a(cy,er,N=1.):
    return lambda t: N*kopen(cy(t),er(t))
    
def rho(N,cy,er,frames):
    aa   = a(cy,er,N)
    rho = lambda t: aa(t)*np.exp(-scipy.integrate.quadrature(aa,0,t,maxiter = 200)[0])
    return frames,map(rho,frames)
    
def cdf(N,cy,er,frames):
    aa   = a(cy,er,N)
    P = lambda t: 1. - np.exp(-scipy.integrate.quadrature(aa,0,t,maxiter = 200)[0])
    return frames,map(P,frames)
#def rho():
    
'''
Off[NIntegrate::nlim]
a[cy_,er_]:=Function[t,n open[cy[t],er[t]]]
P[cy_,er_]:=Function[t,1- Exp[-NIntegrate[a[cy,er][x],{x,0,t},PrecisionGoal->2,AccuracyGoal->2]]]
\[Rho][cy_,er_]:=Function[t,a[cy,er][t]*Exp[-NIntegrate[a[cy,er][x],{x,0,t}]]]
T[cy_,er_,tend_]:=NIntegrate[t \[Rho][cy,er][t],{t,0,tend}]

End[]
EndPackage[]
$ContextPath = DeleteCases[$ContextPath, "RyRModel`"]
'''

def overview(ax1,data):
    t=data.events()['t']
    noch=data.events()['noch']
    ax1.plot(t,noch,lw=1,c='black',label = '# open',drawstyle = 'steps-post')
    ax1.legend(loc='upper left')
