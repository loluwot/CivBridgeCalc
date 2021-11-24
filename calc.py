import numpy as np
from itertools import accumulate
from bisect import bisect_left, bisect_right
from matplotlib import pyplot as plt
from numpy.core import numeric
def get_ybar(bases, lengths):
    endings = list(accumulate(lengths, lambda x, y: x + y))
    centroids = list(lambda i: endings[i] - lengths[i]/2), map(range(len(endings)))
    areas = list(map(lambda x, y: x*y, lengths, bases))
    ybar = list(map(lambda x, y: x*y, areas, centroids))/sum(areas)
    return ybar

class CrossSection: #assuming can be simplified to horizontally symmetric group of positive and negative area rectangles
    def __init__(self, lengths, bases) -> None:
        self.lengths = lengths
        self.endings = list(accumulate(lengths, lambda x, y: x + y))
        self.bases = bases
        self.centroids = list(map(lambda i: self.endings[i] - self.lengths[i]/2, range(len(self.endings))))
        self.areas = list(map(lambda x, y: x*y, self.lengths, self.bases))
        self.ybar = sum(list(map(lambda x, y: x*y, self.areas, self.centroids)))/sum(self.areas)
        Isects = list(map(lambda x, y: x*y**3/12, self.bases, self.lengths))
        Iadj = list(map(lambda a, c: a*(self.ybar - c)**2, self.areas, self.centroids))
        self.I = sum(Isects) + sum(Iadj)
        self.Q = list(map(lambda a, c: a*c, self.areas, self.centroids))

    def get_Qy(self, y1): #y0 assumed to be bottom 
        idx = bisect_left(self.endings, y1)
        new_section = CrossSection(self.lengths[:idx] + [self.lengths[idx] - (self.endings[idx] - y1)], self.bases[:idx+1])
        # print(new_section.lengths)
        # print(new_section.Q)
        return sum(new_section.Q)
    
    def get_shear_per_force(self, y1):
        idx = bisect_left(self.endings, y1)
        return self.get_Qy(y1)/self.bases[idx]/self.I
        pass
class Bridge: #Constant x thickness, non constant y thickness hollow member
    def __init__(self, cross, L, applied_loads, reaction_locs) -> None: #assume 2 reaction locations
        self.L = L
        self.cross = cross
        self.applied_loads = applied_loads #location, forces
        self.reaction_locs = reaction_locs
        moments = sum(list(map(lambda x: x[1]*(x[0] - self.reaction_locs[0]), self.applied_loads)))
        rxn2 = moments/(self.reaction_locs[1] - self.reaction_locs[0])
        self.reaction_forces = [sum([x[1] for x in self.applied_loads]) - rxn2, rxn2]
        self.reaction_loads = list(zip(self.reaction_locs, self.reaction_forces))
    def shear_stress(self, y, L):
        all_loads = self.applied_loads + self.reaction_loads
        all_loads = sorted(all_loads, key=lambda x: x[0])
        locs, forces = zip(*all_loads)
        location = bisect_left(locs, L)
        nforces = list(accumulate(forces, lambda x, y: x + y))
        V = nforces[location]
        return self.cross.get_shear_per_force(y)*V

A = CrossSection([100-1.27,1.27, 1.27*2], [2*1.27, 20,120])
print(A.I)
print(A.ybar)
bridge_test = Bridge(A, 550+510+190, [(550, 100), (550+510+190, 100)], [0, 550+510])
xdata = np.linspace(0, 75 - 1.27, num=1000)
ydata = []
for X in xdata:
    # print(bridge_test.shear_stress(75 - 1.27, L))
    ydata.append(bridge_test.shear_stress(X, (550+510+190)/3))
plt.plot(xdata, ydata)
plt.show()