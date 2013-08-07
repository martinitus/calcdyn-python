from scipy.stats import nanmean
from scipy.stats import nanstd
import numpy as np

class Spine(object):
    def __init__(self, data):
        self.data = data;
        
    def __repr__(self):
        return repr(self.data)
    
    def all(self, normed = False):
        if normed:
            return self.data / self.premean()
        else:
            return self.data;
    
    def pre(self, normed = False):
        if normed:
            return self.data[0:6] / self.premean()
        else:
            return self.data[0:6]
    
    def post(self, normed = False):
        if normed:
            return self.data[7:] / self.premean()
        else:
            return self.data[7:]
        
    
    def premean(self):
        return nanmean(self.data[0:6])
    
    def postmean(self):
        return nanmean(self.data[7:])
    
    def boosting(self):
        return nanmean(self.post())/nanmean(self.pre())


class SpineCollection(object):
    def __init__(self, path):        
        self.raw = np.genfromtxt(path,delimiter=",")
        self.data = [Spine(self.raw[:,i]) for i in range(self.raw.shape[1])]
    
    def spine(self,i):
        return self.data[i]
    
    def spines(self, condition = lambda x: True):
        return [s for s in self.data if condition(s)]
    
    def pre(self, condition = lambda x: True, normed = False):
       return np.array([val for spine in self.spines(condition) for val in spine.pre(normed) if not np.isnan(val)])
    
    def post(self, condition = lambda x: True, normed = False):
        return  np.array([val for spine in self.spines(condition) for val in spine.post(normed) if not np.isnan(val)])
    
    def all(self, condition = lambda x: True, normed = False):
        return  np.array([val for spine in self.spines(condition) for val in spine.all(normed) if not np.isnan(val)])
        