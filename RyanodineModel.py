import StochasticModel
import numpy
import math


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
