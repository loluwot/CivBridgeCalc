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
        # self.Q = list(map(lambda a, c: a*c, self.areas, self.centroids))

    def get_Qy(self, y1): #y0 assumed to be bottom 
        idx = bisect_left(self.endings, y1)
        # print(idx)
        new_section = CrossSection(self.lengths[:idx] + [self.lengths[idx] - (self.endings[idx] - y1)], self.bases[:idx+1])
        # print(new_section.lengths)
        # print(list(map(lambda a, c: a*(c - self.ybar), new_section.areas, new_section.centroids)))
        return sum(list(map(lambda a, c: a*(c - self.ybar), new_section.areas, new_section.centroids)))
        # return sum(new_section.Q)
    
    def get_shear_per_force(self, y1):
        idx = bisect_left(self.endings, y1)
        return self.get_Qy(y1)/self.bases[idx]/self.I
class Bridge: #Constant x thickness, non constant y thickness hollow member
    def __init__(self, cross, L, applied_loads, reaction_locs) -> None: #assume 2 reaction locations
        self.L = L
        self.cross = cross #list of (locs_end, crosssections)
        self.applied_loads = applied_loads #location, forces
        self.reaction_locs = reaction_locs
        moments = sum(list(map(lambda x: x[1]*(x[0] - self.reaction_locs[0]), self.applied_loads)))
        rxn2 = moments/(self.reaction_locs[1] - self.reaction_locs[0])
        self.reaction_forces = [-(sum([x[1] for x in self.applied_loads]) - rxn2), -rxn2]
        self.reaction_loads = list(zip(self.reaction_locs, self.reaction_forces))
        all_loads = self.applied_loads + self.reaction_loads
        self.all_loads = sorted(all_loads, key=lambda x: x[0])

    def sfd(self):
        return [(0,0)] + list(accumulate(self.all_loads, lambda x, y: (y[0], x[1] + y[1])))+ [(self.L, 0)] 

    def plot_sfd(self):
        net_loads = self.sfd()    
        new_plot = [(net_loads[i][0], net_loads[i-1][1]) for i in range(1, len(net_loads))]
        tot = [net_loads, new_plot]
        interlaced = [tot[i % 2][i//2] for i in range(len(new_plot) + len(net_loads))]
        print(interlaced)
        xdata, ydata = zip(*interlaced)
        plt.plot(xdata, ydata)
        plt.show()
        
    def bmd(self):
        sf = self.sfd()
        bmd = [(sf[i][0], sf[i-1][1]*(sf[i][0] - sf[i-1][0])) for i in range(1, len(sf))]
        bmd = [(0,0)] + list(accumulate(bmd, lambda x, y: (y[0], y[1] + x[1])))
        return bmd

    def plot_bmd(self):
        xdata, ydata = zip(*self.bmd())
        plt.plot(xdata, ydata)
        plt.show()

    def shear_stress(self, y, L):
        locs, forces = zip(*self.all_loads)
        location = bisect_left(locs, L)
        nforces = list(accumulate(forces, lambda x, y: x + y))
        V = nforces[location]
        return self.cross.get_shear_per_force(y)*V

    

A = CrossSection([10, 100, 10], [80, 2*1.27, 120])
print(A.I)
print(A.ybar)
bridge_test = Bridge(A, (50+150+400)*2, [(50, 100/2), (50+150+400, 100), ((50+150+400)*2 - 50, 100/2)], [200, (50+150+400)*2 - 200])
print(bridge_test.reaction_loads)
bridge_test.plot_sfd()
bridge_test.plot_bmd()

# xdata = np.linspace(0, 120, num=1000)
# ydata = []
# for X in xdata:
#     print(bridge_test.shear_stress(X, bridge_test.L/2))
#     ydata.append(bridge_test.shear_stress(X, bridge_test.L/2))
# plt.plot(xdata, ydata)
# plt.show()