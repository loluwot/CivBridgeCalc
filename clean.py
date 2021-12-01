import itertools 
import matplotlib 
import numpy as np 
from itertools import accumulate
from functools import reduce
from bisect import bisect ,bisect_left ,bisect_right 
from matplotlib import pyplot as plt 
import math 
from matplotlib .patches import Rectangle 
import matplotlib .patches as mpat 
from matplotlib .patches import FancyArrowPatch ,Polygon
import matplotlib.animation as animation


from sympy import *
from sympy .solvers import solve 
from sympy.plotting import plot
from decimal import Decimal 
import inspect
from tqdm import tqdm
from decimal import *

# getcontext().prec = 

def to_decimal(x):
    if type(x) == list:
        if len(x) == 0:
            return []
        if type(x[0]) == tuple:
            return list(map(lambda v: (Decimal(str(v[0])), Decimal(str(v[1]))), x))
        elif type(x[0]) == int or type(x[0]) == float:
            return list(map(lambda v: Decimal(str(v)), x))
        else:
            return x
    elif type(x) == int or type(x) == float or type(x) == np.float64:
        return Decimal(str(x))
    # if type(x) != Decimal:
        # print(x, type(x))
    return x

def fix_legend(names):
    ax = plt.gca()
    for i, p in enumerate(ax.get_lines()):    
        if p.get_label() in names[:i]:    
            idx = names.index(p.get_label())     
            p.set_c(ax.get_lines()[idx].get_c())  
            p.set_label('_' + p.get_label())  

def get_ybar (bases ,lengths ):
    endings =list (accumulate (lengths ,lambda x ,y :x +y ))
    centroids =list (lambda i :endings [i ]-lengths [i ]/Decimal ('2')),map (range (len (endings )))
    areas =list (map (lambda x ,y :x *y ,lengths ,bases ))
    ybar =list (map (lambda x ,y :x *y ,areas ,centroids ))/sum (areas )
    return ybar 

class Material :
    def __init__ (self ,E ,max_tension_stress ,max_compression_stress ,max_shear_stress ,u , max_glue_stress)->None :
        names = inspect.getfullargspec(self.__init__)[0]
        for name in names:
            if name == 'self':
                continue
            setattr(self, name, to_decimal(locals()[name]))
        # self .E = to_decimal(E) 
        # self .max_tension_stress = max_tension_stress 
        # self .max_compression_stress =max_compression_stress 
        # self .max_shear_stress =max_shear_stress 
        # self .u =u 
        
        pass 

def lerp (L0 ,X0 ,L1 ,X1 ,L ):
    if L < L0:
        return X0
    elif L > L1:
        return X1
    return to_decimal(X0) +to_decimal(L -L0 )*to_decimal(X1 -X0 )/to_decimal(L1 -L0 )

def sign (x ):
    return Decimal ('0')if x ==Decimal ('0')else abs (x )/x 

class CrossSection :#assuming can be simplified to horizontally symmetric group of positive and negative area rectangles
    def __init__ (self ,lengths ,bases ,xdisp, glue)->None :
        names = inspect.getfullargspec(self.__init__)[0]
        for name in names:
            if name == 'self':
                continue
            setattr(self, name, to_decimal(locals()[name]))
        self.initialize ()

    def initialize (self ):
        self .endings =list (accumulate (self .lengths ,lambda x ,y :x +y ))
        # self.endings = list(map(lambda x: float("{:.5f}".format(x)), self.endings))
        self .starting =list (accumulate ([0 ]+self .lengths ,lambda x ,y :x +y ))[:-1 ]
        # self.starting = list(map(lambda x: float("{:.5f}".format(x)), self.starting))
        self .centroids =[self .endings [i ]-self .lengths [i ]/Decimal ('2')for i in range (len (self .endings ))]
        self .areas =[x *y for x ,y in zip (self .lengths ,self .bases )]
        self .ybar =sum ([x *y for x ,y in zip (self .areas ,self .centroids )])/sum (self .areas )
        self .I =sum ([x *y **Decimal ('3')/Decimal ('12')for x ,y in zip (self .bases ,self .lengths )])+sum ([a *(self .ybar -c )**Decimal ('2')for a ,c in zip (self .areas ,self .centroids )])

    def copy (self ):
        return CrossSection (self .lengths ,self .bases ,self .xdisp, self.glue)

    def get_Qy (self ,y1 ):#y0 assumed to be bottom 
        idx = bisect_left (self .endings, y1 )
        new_section =CrossSection (self .lengths [:idx ]+[self .lengths [idx ]-(self .endings [idx ]-y1 )],self .bases [:idx +1 ],self .xdisp [:idx +1 ], self.glue)
        return sum ([a *(c -self .ybar )for a ,c in zip (new_section .areas ,new_section .centroids )])

    def get_shear_per_force (self ,y1 ):
        idx = bisect_left (self .endings ,y1 )
        # print(y1, self.endings, idx, self.bases[idx])
        return self .get_Qy (y1 )/self .bases [idx ]/self .I 

    def draw (self ):
        ax =plt .gca ()
        for b ,l ,s ,x in zip (self .bases ,self .lengths ,self .starting ,self .xdisp ):
            if x !=Decimal ('0'):
                ax .add_patch (Rectangle ((x -b /Decimal ('4'),s ),b /Decimal ('2'),l ,linewidth =1 ,edgecolor ='r',facecolor ='none'))
                ax .add_patch (Rectangle ((-x -b /Decimal ('4'),s ),b /Decimal ('2'),l ,linewidth =1 ,edgecolor ='r',facecolor ='none'))
            else :
                ax .add_patch (Rectangle ((-b /Decimal ('2'),s ),b ,l ,linewidth =1 ,edgecolor ='r',facecolor ='none'))
        ax .relim ()
        ax .autoscale_view ()
        plt .show ()

    def split (self ):
        new_lens =self .lengths [::]
        new_bases =self .bases [::]
        new_xdisp =self .xdisp [::]
        for i in range (len (self .bases )):
            if self .starting [i ]<self .ybar <self .endings [i ]:
                new_bases .insert (i ,new_bases [i ])
                new_xdisp .insert (i ,new_xdisp [i ])
                new_lens .insert (i ,self .ybar -self .starting [i ])
                new_lens [i +1 ]=self .endings [i ]-self .ybar 
                break 
        return CrossSection (new_lens ,new_bases ,new_xdisp, self.glue)

class CrossGroup :
    def __init__ (self ,crosses ,peaks ,changed_members ,starting )->None : 
        names = inspect.getfullargspec(self.__init__)[0]
        for name in names:
            if name == 'self':
                continue
            setattr(self, name, to_decimal(locals()[name]))
        self .changed_members =changed_members 

    def get_cross (self ,L ):

        idx = bisect_right (self .starting, L)-1 
        changed_member =self .changed_members [idx ]
        locs ,peaks =zip (*self .peaks )
        idx2 =bisect_right (locs ,L )-1
        # print(idx2, L, locs) 
        if idx2 + 1 >= len(peaks):
            height = peaks[idx2]
        else:
            height =lerp (locs [idx2 ],peaks [idx2 ],locs [idx2 +1 ],peaks [idx2 +1 ],L )
        # height = peaks[idx2] + (peaks[idx2 + 1] - peaks[idx2])*(L - locs[idx2])/(locs[idx2 + 1] - locs[idx2])
        cr1 =self .crosses [idx ].copy ()
        cr1 .lengths [changed_member ]=height *self .crosses [idx ].lengths [changed_member ]
        cr1 .initialize ()
        return cr1 

class Bridge :#Constant x thickness, non constant y thickness hollow member
    def __init__ (self ,cross ,L ,applied_loads ,reaction_locs ,material ,diaphragms )->None :#assume 2 reaction locations
        names = inspect.getfullargspec(self.__init__)[0]
        for name in names:
            if name == 'self':
                continue
            setattr(self, name, to_decimal(locals()[name]))
        self.initialize()

    def initialize(self):
        total_force = 1 if sum(list (map (lambda x :x [1 ], self .applied_loads ))) == 0 else sum(list (map (lambda x :x [1 ], self .applied_loads )))
        average_force = sum(list (map (lambda x :x [1 ]*x[0], self .applied_loads )))/Decimal(total_force)
        average_force = (average_force, sum(list (map (lambda x :x [1 ], self .applied_loads ))))
        # moments =sum (list (map (lambda x :x [1 ]*(x [0 ]-self .reaction_locs [0 ]),self .applied_loads )))
        moments = average_force[1] * (average_force[0] - self.reaction_locs[0])
        rxn2 =moments / (self .reaction_locs [1 ]-self .reaction_locs [0 ])
        self.reaction_forces =[-(sum ([x [1 ]for x in self .applied_loads ])-rxn2 ),-rxn2 ]
        # print(self.reaction_forces)
        self.reaction_loads =list (zip (self .reaction_locs ,self .reaction_forces ))
        all_loads =self .applied_loads +self .reaction_loads 
        self .all_loads =sorted (all_loads ,key =lambda x :x [0 ])

    def sfd (self ):
        return list (accumulate (self .all_loads ,lambda x ,y :(y [0 ],x [1 ]+y [1 ])))

    def plot_sfd (self ):
        net_loads =self .sfd ()
        # print('NET LOADS', net_loads)
        new_plot =[(net_loads [i ][0 ],net_loads [i -1 ][1 ])for i in range (1 ,len (net_loads ))]
        # print('NEW PLOT', new_plot)
        tot =[net_loads ,new_plot ]
        interlaced =[(Decimal ('0'),0 )]+[tot [i % 2][i // 2]for i in range (len (new_plot )+len (net_loads ))]
        # print(interlaced)
        xdata ,ydata =zip (*interlaced )
        plt .plot (xdata ,ydata , 'b', label='SFD')
        # plt .show ()

    def bmd (self ):
        sf =self .sfd ()
        # print('SFD', sf)
        bmd =[(sf [i ][0 ],sf [i -1 ][1 ]*(sf [i ][0 ]-sf [i -1 ][0 ]))for i in range (1 ,len (sf ))]
        bmd =[(Decimal ('0'),0 )]+list (accumulate (bmd ,lambda x ,y :(y [0 ],y [1 ]+x [1 ])))
        if bmd[-1] != (Decimal(1280), Decimal(0)):
            bmd.append((Decimal(1280), Decimal(0)))
        return bmd 

    def get_moment (self ,L ):
        locs, moments =zip (*self .bmd ())
        idx =bisect_right (locs ,L )-1 
        # print(locs, L, idx)
        if idx + 1 >= len(locs):
            return moments[idx]
        ans =lerp (locs [idx ],moments [idx ],locs [idx +1 ],moments [idx +1 ],L )
        # print('ANS', ans)
        return ans 
        # return moments[idx] + (moments[idx + 1] - moments[idx])*(L - self.locs[idx])/(self.locs[idx + 1] - self.locs[idx])

    def get_shear(self, L):
        locs ,shears =zip (*self .sfd ())
        # print('shear ---------------', L)
        
        idx =bisect_right (locs ,L )-1 
        # print(shears[idx])
        return shears[idx]

    def get_max_moment (self ,L0 ,L1 ):
        locs ,moments =zip (*self .bmd ())
        idx1 =bisect_right (locs ,L0 )-1 
        idx2 =bisect_right (locs ,L1 )-1 
        LL =[]
        if idx1 !=idx2 :
            LL =list (zip (moments [idx1 +1 :idx2 +1 ],locs [idx1 +1 :idx2 +1 ]))
        LL +=[(self .get_moment (L0 ),L0 ),(self .get_moment (L1 ),L1 ),(Decimal ('0'),-1 )]
        ex_min ,ex_max =min (LL ),max (LL )
        extras =list (filter (lambda v :(v [0 ]==ex_min [0 ])or (v [0 ]==ex_max [0 ]),LL ))
        return extras 

    def plot_bmd (self ):
        BMD_plot =self .bmd ()
        # print(BMD_plot)
        xdata ,ydata =zip (*BMD_plot )
        artist, = plt .gca().plot (xdata ,ydata, label='BMD')
        # plt .show ()
        return artist
    def curvature (self ):
        bm =self .bmd ()
        return list (map (lambda x :(x [0 ],x [1 ]/self .material .E /self .cross .get_cross (x [0 ]).I *Decimal ('1000')),bm ))

    def plot_curvature (self ):
        xdata ,ydata =zip (*self .curvature ())
        plt .plot (xdata ,ydata )
        plt .show ()

    def flex_stress (self ,y ,L ):
        M =self .get_moment (L )
        return M *(self.cross.get_cross (L ).ybar -y )/self .cross .get_cross (L ).I 

    def shear_stress (self ,y ,L ):
        locs ,forces =zip (*self .all_loads )
        location =bisect_right (locs ,L )-1 
        # print(location, L, locs)
        nforces =list (accumulate (forces ,lambda x ,y :x +y ))
        V =nforces [location ]
        return self .cross .get_cross (L ).get_shear_per_force (y )*V 

    def plate_buckling (self ,L ):
        cur_cross =self .cross .get_cross (L ).split ()
        # cur_cross.draw()
        def isbounded (idx ):
            return Decimal ('0')<=idx <len (cur_cross .bases )
        mcrit =math .inf 
        for i in range (len (cur_cross .bases )):
            # print(i)
            is_compressed =not ((cur_cross .starting [i ]>=cur_cross .ybar )^(self .get_moment (L )>Decimal ('0')))
            b ,l =cur_cross .bases [i ],cur_cross .lengths [i ]
            if not is_compressed :
                continue 
            if l >b :
                crit =(b /Decimal ('2')/l )**Decimal ('2')*Decimal ('6')*Decimal(math .pi) **Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
                mcrit =min (mcrit ,crit )
                # print('POSS CRIT', crit)
            else :
                neighbours =list (filter (isbounded ,[i -1 ,i +1 ]))
                crit =None 
                poss =list (filter (lambda x :cur_cross .xdisp [x ]!=Decimal ('0'),neighbours ))
                if len (poss ) >Decimal ('0'):
                    # print('CASE 1 AND 2')
                    b2 =Decimal ('2')*cur_cross .xdisp [poss [0 ]]-cur_cross .bases [poss [0 ]]/Decimal ('2')#two held end
                    oneheld =(b -b2 -cur_cross .bases [poss [0 ]])/Decimal ('2')
                    twoheld =b2 
                    if oneheld < Decimal ('0') or twoheld < Decimal ('0'):
                        continue
                    # print(twoheld, oneheld, l)
                    crit = math.inf if oneheld == 0 else Decimal ('0.425')*(l /(oneheld ))**Decimal ('2')*Decimal(math .pi)**Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
                    crit2 =math.inf if twoheld == 0 else Decimal ('4')*(l /(twoheld ))**Decimal ('2')*Decimal(math .pi)**Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
                    # print('POSS CRITS', crit, crit2)
                    crit =min (crit ,crit2 )
                else :
                    # print('CASE 2')
                    crit =Decimal ('0.425')*(l /b )**Decimal ('2')*Decimal(math .pi) **Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
                    # print('POSS CRITS', crit)
                mcrit =min (mcrit ,crit )        
        return mcrit 


    def shear_buckling (self ,L ):
        cross =self .cross .get_cross (L )
        min_shear = Decimal(math .inf) 
        # a = max (list (map (lambda i :self .diaphragms [i ]-self .diaphragms [i -1 ],range (1 ,len (self .diaphragms )))))
        aidx = bisect_right(self.diaphragms, L) - 1
        if aidx+1 >= len(self.diaphragms):
            a = self.L - self.diaphragms[aidx]
        elif aidx < 0:
            a = self.diaphragms[aidx+1]
        else:
            a = self.diaphragms[aidx+1] - self.diaphragms[aidx]
        print('A', a, aidx)
        for b ,l ,x in zip (cross .bases ,cross .lengths ,cross .xdisp ):
            if x !=Decimal ('0'):
                # print ('B, L',b /Decimal ('2'),l )
                min_shear =min (Decimal ('5')*Decimal(math .pi) **Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))*((b /Decimal ('2')/l )**Decimal ('2')+(b /Decimal ('2')/a )**Decimal ('2')),min_shear )
            else :
                # print ('B, L',b ,l )
                min_shear =min (Decimal ('5')*Decimal(math .pi) **Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))*((b /l )**Decimal ('2')+(b /a )**Decimal ('2')),min_shear )
        return min_shear 

    def draw (self ):
        ax =plt .gca ()
        # print(self.cross.peaks)
        points =[[0 ,0 ]]+[[x [0 ],-x [1 ]]for x in self .cross .peaks ]+[[self .L ,0 ]]
        points =np .array (points )
        # print(points.shape)
        ax.add_patch (mpat .Polygon (points ))
        for loc ,force in self .all_loads :
            ax.add_patch (FancyArrowPatch ((loc ,-force),(loc ,Decimal ('0')),mutation_scale =10))
        ax .relim ()
        ax .autoscale_view ()
        plt .show ()
        pass 
        
    def max_shear_forces(self, x, use_cross=False, idx=None):
        cross = self.cross.crosses[idx] if use_cross else self.cross.get_cross(x)
        maxes = []
        # matboard shear failure
        # print('Max Qy:', cross.get_Qy(cross.ybar), 'Force: ', cross.get_shear_per_force(cross.ybar))
        maxes.append(abs(self.material.max_shear_stress/cross.get_shear_per_force(cross.ybar)))

        # print('Max', maxes[-1])
        #glue shear failure
        # print('Qy glue:', [cross.get_Qy(g) for g in cross.glue], 'Force: ', [cross.get_shear_per_force(g) for g in cross.glue])
        maxes.append(min([abs(self.material.max_glue_stress/cross.get_shear_per_force(g)) for g in cross.glue]))
        #shear buckling failure
        # print(type(self.shear_buckling(x)), self.shear_buckling(x))
        maxes.append(abs(self.shear_buckling(x)/cross.get_shear_per_force(cross.ybar)))
        return maxes
    
    def plot_annotated_sfd(self):
        self.plot_sfd()
        xaxis, _ = zip(*self.sfd())
        xaxis = list(xaxis)
        xaxis += self.diaphragms + self.cross.starting
        xaxis = sorted(xaxis)
        yaxes = [[] for _ in range(6)]
        labels = ['Matboard Shear Failure', 'Glue Shear Failure', 'Shear Buckling Failure']
        res = [None]*6
        res[::2] = labels
        res[1::2] = labels
        for x in xaxis:
            ys = []
            for V in self.max_shear_forces(x):
                # print(V)
                ys.extend([V, -V])
                
            yaxes = [a + [b] for a, b in zip(yaxes, ys)]
        for yax, labs in zip(yaxes, res):
            net_loads = list(zip(xaxis, yax))
            new_plot =[(net_loads [i ][0 ],net_loads [i -1 ][1 ])for i in range (1 ,len (net_loads ))]
            interlaced = [None]*(len(net_loads) + len(new_plot))
            interlaced[::2] = net_loads
            interlaced[1::2] = new_plot
            xax, yax = zip(*interlaced)
            plt.plot(xax, yax, label=labs)
            # plt.plot(xaxis, yax, label=labs)
        fix_legend(['SFD'] + res)
        plt.xlabel('Location on Bridge (mm)')
        plt.ylabel('Shear force (N)')
        plt.gca().legend()
        

    def plot_annotated_bmd(self):
        artists = [self.plot_bmd()]
        # print(self.bmd())
        xaxis, _ = zip(*self.bmd())
        xaxis = list(xaxis)
        xaxis2, _ = zip(*self.sfd())
        xaxis += xaxis2
        xaxis += self.diaphragms + self.cross.starting
        xaxis += list(map(Decimal, np.linspace(0, float(self.L), num=100)))
        xaxis = sorted(set(xaxis))
        yaxes = [[] for _ in range(4)]
        res = ['Tension Failure', 'Compression Failure', 'Plate Buckling Failure']

        for x in xaxis:
            ys = []
            for V, R in zip(self.max_bending_moments(x), res):
                # print(f'MAX {R}', V)
                ys.append(V)
            yaxes = [a + [b] for a, b in zip(yaxes, ys)]
        # artists = []
        for yax, labs in zip(yaxes, res):
            net_loads = list(zip(xaxis, yax))
            new_plot =[(net_loads [i ][0 ],net_loads [i -1 ][1 ])for i in range (1 ,len (net_loads ))]
            interlaced = [None]*(len(net_loads) + len(new_plot))
            interlaced[::2] = net_loads
            interlaced[1::2] = new_plot
            xax, yax = zip(*interlaced)
            artist, = plt.gca().plot(xax, yax, label=labs)
            artists.append(artist)

        fix_legend(['BMD'] + res)
        plt.xlabel('Location on Bridge (mm)')
        plt.ylabel('Moment (Nmm)')
        plt.gca().legend()
        return artists

    def max_bending_moments(self, x, use_cross=False, idx=None):
        S = int(sign(self.get_moment(x)))
        # print(S)
        if S == 0:
            return [Decimal(math.inf) for _ in range(3)]
        cross = self.cross.crosses[idx] if use_cross else self.cross.get_cross(x)
        y_opts = [cross.ybar - sum(cross.lengths), cross.ybar]
        # print(y_opts)
        maxes = []
        #tension failure
        y = y_opts[(S + 1)//2]
        maxes.append(self.material.max_tension_stress*cross.I/y)
        # print(maxes[0], cross.I, y)
        #compression failure
        y = y_opts[(1 - S)//2]
        # print(y)
        maxes.append(self.material.max_compression_stress*cross.I/y)
        
        #flexural buckling failure
        maxes.append(-self.plate_buckling(x)*cross.I/y)
        if abs(maxes[-1])//10000 == 4:
            print(S)
            print(self.plate_buckling(x), cross.I, y, maxes[-1])
        # print('MOMENTS', maxes)
        return maxes

    def bridge_capacity(self):
        tests1 = set(self.diaphragms + self.cross.starting)
        res1 = [self.max_shear_forces(test) for test in tests1]
        def reduce_mins(res1):
            return reduce(lambda x, y: [a if abs(a) < abs(b) else b for a, b in zip(x, y)], res1)
        lowest_shear = reduce_mins(res1)
        tests2 = map(Decimal, np.linspace(0, float(self.L), num=100))
        res2 = [self.max_bending_moments(test) for test in tests2]
        lowest_moment = reduce(lambda x, y: [a if abs(a) < abs(b) else b for a, b in zip(x, y)], res2)
        lowest_moment_pos = reduce_mins(map(lambda x: [math.inf if v < 0 else v for v in x], res2))
        lowest_moment_neg = reduce_mins(map(lambda x: [math.inf if v >= 0 else v for v in x], res2))
        print('Shear Failures:')
        print('\n'.join(list(map(lambda x, y: f'{x}: {y} N', ['Matboard Shear Failure', 'Glue Shear Failure', 'Shear Buckling Failure'], lowest_shear))))
        print('Moment Failures:')
        print('\n'.join(list(map(lambda x, y: f'{x} (Positive): {y} Nmm', ['Tension Failure', 'Compression Failure', 'Plate Buckling Failure'], lowest_moment_pos))))
        print('\n'.join(list(map(lambda x, y: f'{x} (Negative): {y} Nmm', ['Tension Failure', 'Compression Failure', 'Plate Buckling Failure'], lowest_moment_neg))))
        return lowest_shear, lowest_moment_pos, lowest_moment_neg
    
    def failure_force_at_x (self, x):
        # print('X', x)
        SFs = self.max_shear_forces(x)
        BMs = self.max_bending_moments(x)
        shearP = math.inf if self.get_shear(x) == 0 else min(SFs)/abs(self.get_shear(x))
        M1 = self.get_moment(x)
        # possible_BM = list(map(abs, filter(lambda V: sign(V) == sign(M1), self.max_bending_moments(x))))
        # print(M1)
        momentP = math.inf if M1 == 0 else sign(M1)*min(map(abs, filter(lambda V: sign(V) == sign(M1), BMs)))/M1
        # print(momentP)
        # print('SHEAR ----', shearP, 'MOMENT', momentP)
        
        return min(abs(shearP), abs(momentP)), shearP, momentP, SFs.index(min(SFs)), BMs.index(min(BMs))

    def failure_force (self):
        test_points = list(set([v[0] for v in self.sfd()] + [v[0] for v in self.bmd()] + self.cross.starting + self.diaphragms))
        test_points += list(map(Decimal, np.linspace(0, float(self.L), num=50)))
        test_points = sorted(test_points)
        all_fails = list(map(self.failure_force_at_x, test_points))
        # plt.plot(test_points, all_fails)
        # plt.show()
        # print(all_fails)
        # self.plot_bmd()
        # print(min(all_fails))
        return min(all_fails)[0]
        # print('l', all_fails_l)
        #         return min(all_fails)[0]
        # return 

    def displacement (self):
        SS = self.sfd() + [(self.L, 0)]
        print(SS)
        z = Symbol('z')
        defs = [(SS[i][1], (z >= SS[i][0]) & (z < SS[i+1][0])) for i in range(len(SS) - 1)] + [(0, True)]
        p = Piecewise(*defs)
        print(self.cross.starting)
        cross_sections = [(self.cross.crosses[i].I*self.material.E, (z >= self.cross.starting[i]) & (z < (self.cross.starting[i + 1] if i + 1 < len(self.cross.starting) else self.L))) for i in range(len(self.cross.starting))] + [(1, True)]
        print(cross_sections)
        cross_sections = Piecewise(*cross_sections)
        # print(p)
        # print(cross_sections)
        plot(p, (z, 0, self.L), show=True)
        BMD = integrate(p, (z, 0, z))
        BMD = BMD - BMD.subs(z, self.L)
        plot(BMD, (z, 0, self.L), show=True)
        # self.plot_bmd()
        # plt.show()
        curvature = BMD/cross_sections
        plot(curvature, (z, 0, self.L), show=True)
        slopes = integrate(curvature, (z, 0, z))
        print(slopes)
        disp = integrate(slopes, (z, 0, z))
        print(disp)
        y1 = disp.subs(z, self.reaction_locs[0])
        y2 = disp.subs(z, self.reaction_locs[1])
        m = (y2 - y1)/(self.reaction_locs[1] - self.reaction_locs[0])
        disp = disp - (m*(z-self.reaction_locs[0]) + y1)
        print(disp.subs(z, 1250/2))
        plot(disp, (z, 0, self.L), show=True, xlabel='Location of Bridge (mm)', ylabel='Displacement (mm)')
        return disp
        


class BridgeTester:
    def __init__(self, PTrain) -> None:
        self.PTrain = PTrain
        
    
    def get_failure_load(self, bridge):
        bridge.applied_loads = [(550, -1), (1250, -1)]
        bridge.initialize()
        return bridge.failure_force()
    
    def trainFOS(self, bridge):
        LBridge = bridge.L
        LTrain = Decimal((52 + 176 + 164 + 176/2)*2)
        train_loadings = list(map(Decimal, np.linspace(-float(LTrain), float(LBridge), num=50)))
        fail_loads = []
    
        artists = []
        for x in (train_loadings):
            # print(x)
            percentage_of_train = min(lerp(x+52, 1, x + LTrain - 52, 0, 0), lerp(x+52, 0, x+LTrain - 52, 1, LBridge))
            a = 52/LTrain
            if percentage_of_train > 1 - a:
                percentage_of_train = 1
            elif percentage_of_train < a:
                percentage_of_train = 0
            # print(percentage_of_train)
            
            P = Decimal('1')/6
            loading = [(52+x, -P), (52 + 176+x, -P), (52 + 176 + 164+x, -P), (52 + 176 + 164 + 176+x, -P), (52 + 176 + 164 + 176 + 164+x, -P), (52 + 176 + 164 + 176 + 164 + 176+x, -P)]
            loading = list(filter(lambda x: 0 <= x[0] < LBridge, loading))
            bridge.applied_loads = loading
            bridge.initialize()
            # print(bridge.bmd())
            P1 = bridge.failure_force()
            
            if percentage_of_train != 0:
                A = bridge.plot_annotated_bmd()
                artists.append(A)
                # plt.show()
            fail_loads.append(P1)
            # return P1/self.PTrain
            # plt.gcf().canvas.draw_idle()
        ani = animation.ArtistAnimation(plt.gcf(), artists, interval=50, blit=True,
                                repeat_delay=1000)
        plt.show()
        return min(fail_loads)/self.PTrain

P = 400/6
LBridge = 550 + 510 + 190+30
loadings = [[(550, -P), (1250, -P)]]
LTrain = (52 + 176 + 164 + 176/2)*2
# print(LTrain < LBridge)
train_loadings = np.linspace(0, LBridge - LTrain)
train_loadings = list(map(lambda x: [(52+x, -P), (52 + 176+x, -P), (52 + 176 + 164+x, -P), (52 + 176 + 164 + 176+x, -P), (52 + 176 + 164 + 176 + 164+x, -P), (52 + 176 + 164 + 176 + 164 + 176+x, -P)], train_loadings))
total_loadings = loadings + train_loadings
LBot = 50
N = 3
HEIGHT = 70
# A = CrossSection([1.27, 75 - 1.27*3, 1.27, 1.27], [80, 2*1.27, 20 + 2*1.27, 100], [0, 40-1.27/2, 35 - 1.27/2, 0], [75 - 1.27])
# A = CrossSection([1.27*2, 1.27, 75 - 1.27*6, 1.27, 1.27*2], [80, 2*1.27 + 20, 2*1.27, 2*1.27 + 20, 100], [0, 35-1.27/2, 40 - 1.27/2, 35-1.27/2, 0], [75 - 1.27*2, 1.27*2])

A = CrossSection([1.27, 1.27, HEIGHT - (3 + 2)*1.27, 1.27, 2*1.27], [LBot, 4*1.27 + 20, 4*1.27, 20 + 4*1.27, 100], [0, LBot/2 - 1.27 - 5 ,LBot/2-1.27,LBot/2-1.27 - 5, 0], [1.27, HEIGHT - 2*1.27])
B = CrossSection([1.27*2, 1.27, HEIGHT - (3 + 2)*1.27, 1.27, 1.27], [100, 20 + 4*1.27, 4*1.27, 20 + 4*1.27, 100], [0, 50 - 1.27 - 5, 50 - 1.27, 50 - 1.27 - 5, 0], [2*1.27, HEIGHT - 1.27])
# A.draw()
# print(A.get_Qy(Decimal(75 - 1.27)))
# B.draw()
# C = CrossGroup([A], [(0, 1), (LBridge, 1)], [1], [0])
C = CrossGroup([A, B], [(0, 1), (LBridge, 1)], [1, 1], [0, 780])

matboard = Material(4000, 30, -6, 4, 0.2, 2)
# fail_loads = []
# train_loadings = np.linspace(-LTrain, LBridge, num=50)
# for x in (train_loadings):
# print(x)
# percentage_of_train = min(lerp(x+52, 1, x + LTrain - 52, 0, 0), lerp(x+52, 0, x+LTrain - 52, 1, LBridge))
# a = 52/LTrain
# if percentage_of_train > 1 - a:
#     percentage_of_train = 1
# elif percentage_of_train < a:
#     percentage_of_train = 0

# print(percentage_of_train)
# P = 400*percentage_of_train/6
# loading = [(52+x, -P), (52 + 176+x, -P), (52 + 176 + 164+x, -P), (52 + 176 + 164 + 176+x, -P), (52 + 176 + 164 + 176 + 164+x, -P), (52 + 176 + 164 + 176 + 164 + 176+x, -P)]
# loading = list(filter(lambda x: 0 <= x[0] < LBridge, loading))
P = 100
# bridge_test = Bridge(C, LBridge, [(550, -P), (1250, -P)], [0, 550 + 510], matboard, [0, 550, 550 + 510, 550 + 510 + 190])
diaphragm_lengths = [15] + [115]*4 + [75] + [170]*3 + [190] + [15]
diaphragms1 = list(accumulate(diaphragm_lengths, lambda x, y: x + y))
# print(diaphragms1)
bridge_test = Bridge(C, LBridge, [(550, -P), (1250, -P)], [0, 550 + 510], matboard, diaphragms1)
# bridge_test.draw()
# bridge_test.displacement()

bridge_test.bridge_capacity()
# bridge_test.plot_annotated_sfd()
bridge_test.plot_annotated_bmd()
plt.show()
bridge_test.plot_annotated_sfd()
# bridge_test.plot_bmd()
plt.show()
bridge_test.displacement()
tester = BridgeTester(400)
print('TRAIN FOS', tester.trainFOS(bridge_test))
print('Failure Load (Case 2)', tester.get_failure_load(bridge_test))
# bridge_test.plot_sfd()
# plt.show()
# print(BridgeTester().get_failure_load(bridge_test))
# print('Bending Moment Peaks: ', bridge_test.bmd()[1:-1])
# bridge_test.draw()
# bridge_test.plot_sfd()
# print(bridge_test.sfd())
# BM = [x[0] for x in bridge_test.bmd()[1:-1]]
# print(bridge_test.failure_force())
# fail_loads.append(bridge_test.failure_force())
# plt.plot(train_loadings, fail_loads)
plt.show()