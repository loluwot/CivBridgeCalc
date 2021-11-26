import itertools
import numpy as np
from itertools import accumulate
from bisect import bisect, bisect_left, bisect_right
from matplotlib import pyplot as plt
from numpy.core import numeric
import math
from matplotlib.patches import Rectangle
def get_ybar(bases, lengths):
    endings = list(accumulate(lengths, lambda x, y: x + y))
    centroids = list(lambda i: endings[i] - lengths[i]/2), map(range(len(endings)))
    areas = list(map(lambda x, y: x*y, lengths, bases))
    ybar = list(map(lambda x, y: x*y, areas, centroids))/sum(areas)
    return ybar

class Material:
    def __init__(self, E, max_tension_stress, max_compression_stress, max_shear_stress, u) -> None:
        self.E = E
        self.max_tension_stress = max_tension_stress
        self.max_compression_stress = max_compression_stress
        self.max_shear_stress = max_shear_stress
        self.u = u
        pass

def lerp(L0, X0, L1, X1, L):
    return X0 + (L - L0)*(X1 - X0)/(L1 - L0)


class CrossSection: #assuming can be simplified to horizontally symmetric group of positive and negative area rectangles
    def __init__(self, lengths, bases, xdisp) -> None:
        self.lengths = lengths
        self.bases = bases
        self.xdisp = xdisp
        self.initialize()

    def initialize(self):
        self.endings = list(accumulate(self.lengths, lambda x, y: x + y))
        self.starting = list(accumulate([0] + self.lengths, lambda x, y: x + y))[:-1]
        self.centroids = [self.endings[i] - self.lengths[i]/2 for i in range(len(self.endings))]
        self.areas = [x*y for x, y in zip(self.lengths, self.bases)]
        self.ybar = sum([x*y for x, y in zip(self.areas, self.centroids)])/sum(self.areas)
        self.I = sum([x*y**3/12 for x, y in zip(self.bases, self.lengths)]) + sum([a*(self.ybar - c)**2 for a, c in zip(self.areas, self.centroids)])

    def copy(self):
        return CrossSection(self.lengths, self.bases, self.xdisp)

    def get_Qy(self, y1): #y0 assumed to be bottom 
        idx = bisect_left(self.endings, y1)
        new_section = CrossSection(self.lengths[:idx] + [self.lengths[idx] - (self.endings[idx] - y1)], self.bases[:idx+1])
        return sum([a*(c - self.ybar) for a, c in zip(new_section.areas, new_section.centroids)])
    
    def get_shear_per_force(self, y1):
        idx = bisect_left(self.endings, y1)
        return self.get_Qy(y1)/self.bases[idx]/self.I

    def draw(self):
        ax = plt.gca()
        for b, l, s, x in zip(self.bases, self.lengths, self.starting, self.xdisp):
            if x != 0:
                ax.add_patch(Rectangle((x - b/4,s),b/2,l,linewidth=1,edgecolor='r',facecolor='none'))
                ax.add_patch(Rectangle((-x - b/4,s),b/2,l,linewidth=1,edgecolor='r',facecolor='none'))
            else:
                ax.add_patch(Rectangle((-b/2,s),b,l,linewidth=1,edgecolor='r',facecolor='none'))
        ax.relim()
        ax.autoscale_view()
        plt.show()

    def split(self):
        new_lens = self.lengths[::]
        new_bases = self.bases[::]
        new_xdisp = self.xdisp[::]
        for i in range(len(self.bases)):
            if self.starting[i] < self.ybar < self.endings[i]:
                new_bases.insert(i, new_bases[i])
                new_xdisp.insert(i, new_xdisp[i])
                new_lens.insert(i, self.ybar - self.starting[i])
                new_lens[i + 1] = self.endings[i] - self.ybar
                break
        return CrossSection(new_lens, new_bases, new_xdisp)

class CrossGroup:
    def __init__(self, crosses, peaks, changed_members, starting) -> None:
        self.crosses = crosses
        self.peaks = peaks #locs, values
        self.changed_members = changed_members
        self.starting = starting
        
    def get_cross(self, L):
        idx = bisect_left(self.starting, L) - 1
        changed_member = self.changed_members[idx]
        locs, peaks = zip(*self.peaks)
        idx2 = bisect_left(locs, L) - 1
        height = lerp(locs[idx2], peaks[idx2], locs[idx2 + 1], peaks[idx2+1], L)
        # height = peaks[idx2] + (peaks[idx2 + 1] - peaks[idx2])*(L - locs[idx2])/(locs[idx2 + 1] - locs[idx2])
        cr1 = self.crosses[idx].copy()
        cr1.lengths[changed_member] = height*self.crosses[idx].lengths[changed_member]
        cr1.initialize()
        return cr1

class Bridge: #Constant x thickness, non constant y thickness hollow member
    def __init__(self, cross, L, applied_loads, reaction_locs, material) -> None: #assume 2 reaction locations
        self.L = L
        self.cross = cross
        self.applied_loads = applied_loads #location, forces
        self.reaction_locs = reaction_locs
        moments = sum(list(map(lambda x: x[1]*(x[0] - self.reaction_locs[0]), self.applied_loads)))
        rxn2 = moments/(self.reaction_locs[1] - self.reaction_locs[0])
        self.reaction_forces = [-(sum([x[1] for x in self.applied_loads]) - rxn2), -rxn2]
        self.reaction_loads = list(zip(self.reaction_locs, self.reaction_forces))
        all_loads = self.applied_loads + self.reaction_loads
        self.all_loads = sorted(all_loads, key=lambda x: x[0])
        self.material = material
        
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

    def get_moment(self, L):
        locs, moments = zip(*self.bmd())
        idx = bisect_left(locs, L) - 1
        return lerp(locs[idx], moments[idx], locs[idx+1], moments[idx+1], L)
        # return moments[idx] + (moments[idx + 1] - moments[idx])*(L - self.locs[idx])/(self.locs[idx + 1] - self.locs[idx])

    def get_max_moment(self, L0, L1):
        locs, moments = zip(*self.bmd())
        idx1 = bisect_left(locs, L0) - 1
        idx2 = bisect_left(locs, L1) - 1
        LL = []
        if idx1 != idx2:
            LL = list(zip(moments[idx1 + 1:idx2 + 1], locs[idx1 + 1: idx2 + 1]))
        LL += [(self.get_moment(L0), L0), (self.get_moment(L1), L1), (0, -1)]
        ex_min, ex_max = min(LL), max(LL)
        extras = list(filter(lambda v: (v[0] == ex_min[0]) or (v[0] == ex_max[0]), LL))
        return extras

    def plot_bmd(self):
        xdata, ydata = zip(*self.bmd())
        plt.plot(xdata, ydata)
        plt.show()

    def curvature(self):
        bm = self.bmd()
        return list(map(lambda x: (x[0], x[1]/self.material.E/self.cross.get_cross(x[0]).I*1000), bm))

    def plot_curvature(self):
        xdata, ydata = zip(*self.curvature())
        plt.plot(xdata, ydata)
        plt.show()

    def shear_stress(self, y, L):
        locs, forces = zip(*self.all_loads) 
        location = bisect_left(locs, L)
        nforces = list(accumulate(forces, lambda x, y: x + y))
        V = nforces[location]
        return self.cross.get_cross(L).get_shear_per_force(y)*V

    def plate_buckling(self, L):
        cur_cross = self.cross.get_cross(L).split()
        cur_cross.draw()
        def isbounded(idx):
            return 0 <= idx < len(cur_cross.bases)
        mcrit = math.inf
        for i in range(len(cur_cross.bases)):
            is_compressed = (cur_cross.starting[i] >= cur_cross.ybar)^(self.get_moment(L) > 0)
            b, l = cur_cross.bases[i], cur_cross.lengths[i]
            print(is_compressed)
            if not is_compressed:
                continue
            if l > b:
                print('CASE 3')
                crit = (b/l)**2*6*math.pi**2*self.material.E/12/(1-self.material.u**2)
                mcrit = min(mcrit, crit)
            else:
                neighbours = list(filter(isbounded, [i - 1, i + 1]))
                crit = None
                poss = list(filter(lambda x: cur_cross.xdisp[x] != 0, neighbours))
                if len(poss) > 0:
                    print('CASE 1 AND 2')
                    b2 = 2*cur_cross.xdisp[poss[0]]
                    crit = 0.425*(l/(b - b2))**2*math.pi**2*self.material.E/12/(1-self.material.u**2)
                    crit = min(crit, 4*(l/b2)**2*math.pi**2*self.material.E/12/(1-self.material.u**2))
                else:
                    crit = 0.425*(l/b)**2*math.pi**2*self.material.E/12/(1-self.material.u**2)
                mcrit = min(mcrit, crit)
        return mcrit

    

A = CrossSection([1.27*2, 100, 1.27*2], [80, 2*1.27, 120], [0,10,0])
A.split().draw()
B = CrossSection([1.27*2, 100, 1.27*2], [80, 2*1.27, 80], [0,10,0])

LBridge = (50+150+400)*2
C = CrossGroup([A, B, A], [(0, 1), (LBridge, 0.5)], [1, 1, 1], [0, 50+150, LBridge - 50 - 150])

matboard = Material(4000, 30, 6, 4, 0.2)
print(A.I)
print(A.ybar)
bridge_test = Bridge(C, (50+150+400)*2, [(50, 100/2), (50+150+400, 100), ((50+150+400)*2 - 50, 100/2)], [200, (50+150+400)*2 - 200], matboard)
print(bridge_test.reaction_loads)
bridge_test.plot_sfd()
bridge_test.plot_bmd()
bridge_test.plot_curvature()
print(bridge_test.plate_buckling(LBridge/2))
print(bridge_test.get_max_moment(0, LBridge/2 - 1))
# xdata = np.linspace(0, 120, num=1000)
# ydata = []
# for X in xdata:
#     print(bridge_test.shear_stress(X, bridge_test.L/2))
#     ydata.append(bridge_test.shear_stress(X, bridge_test.L/2))
# plt.plot(xdata, ydata)
# plt.show()