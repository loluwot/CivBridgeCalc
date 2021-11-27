import itertools 
import matplotlib 
import numpy as np 
from itertools import accumulate 
from bisect import bisect ,bisect_left ,bisect_right 
from matplotlib import pyplot as plt 
import math 
from matplotlib .patches import Rectangle 
import matplotlib .patches as mpat 
from matplotlib .patches import FancyArrowPatch ,Polygon 
from sympy import *
from sympy .solvers import solve 
from decimal import Decimal 
import inspect

def to_decimal(x):
    if type(x) == list:
        if type(x[0]) == tuple:
            return list(map(lambda v: (Decimal(v[0]), Decimal(v[1])), x))
        elif type(x[0]) == int or type(x[0]) == float:
            return list(map(Decimal, x))
        else:
            return x
    elif type(x) == int or type(x) == float:
        return Decimal(x)
    return x

def get_ybar (bases ,lengths ):
    endings =list (accumulate (lengths ,lambda x ,y :x +y ))
    centroids =list (lambda i :endings [i ]-lengths [i ]/Decimal ('2')),map (range (len (endings )))
    areas =list (map (lambda x ,y :x *y ,lengths ,bases ))
    ybar =list (map (lambda x ,y :x *y ,areas ,centroids ))/sum (areas )
    return ybar 

class Material :
    def __init__ (self ,E ,max_tension_stress ,max_compression_stress ,max_shear_stress ,u )->None :
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
    return X0 +(L -L0 )*(X1 -X0 )/(L1 -L0 )

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
        idx =bisect_left (self .endings ,y1 )
        new_section =CrossSection (self .lengths [:idx ]+[self .lengths [idx ]-(self .endings [idx ]-y1 )],self .bases [:idx +1 ],self .xdisp [:idx +1 ], self.glue)
        return sum ([a *(c -self .ybar )for a ,c in zip (new_section .areas ,new_section .centroids )])

    def get_shear_per_force (self ,y1 ):
        idx =bisect_left (self .endings ,y1 )
        # print (y1 ,self .endings ,idx )
        # print (self .bases [idx ])
        # print (self .get_Qy (y1 )/self .bases [idx ]/self .I )
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
        idx =bisect_left (self .starting ,L )-1 
        changed_member =self .changed_members [idx ]
        locs ,peaks =zip (*self .peaks )
        idx2 =bisect_left (locs ,L )-1 
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
        # self.L = L 
        # self.cross =cross 
        # self.applied_loads = applied_loads #location, forces
        # self.reaction_locs = reaction_locs 
        # self.diaphragms = diaphragms 
        moments =sum (list (map (lambda x :x [1 ]*(x [0 ]-self .reaction_locs [0 ]),self .applied_loads )))
        rxn2 =moments /(self .reaction_locs [1 ]-self .reaction_locs [0 ])
        self .reaction_forces =[-(sum ([x [1 ]for x in self .applied_loads ])-rxn2 ),-rxn2 ]
        self .reaction_loads =list (zip (self .reaction_locs ,self .reaction_forces ))
        all_loads =self .applied_loads +self .reaction_loads 

        self .all_loads =sorted (all_loads ,key =lambda x :x [0 ])
        # print (all_loads )
        self .material =material 

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
        plt .plot (xdata ,ydata )
        plt .show ()

    def bmd (self ):
        sf =self .sfd ()
        # print('SFD', sf)
        bmd =[(sf [i ][0 ],sf [i -1 ][1 ]*(sf [i ][0 ]-sf [i -1 ][0 ]))for i in range (1 ,len (sf ))]
        bmd =[(Decimal ('0'),0 )]+list (accumulate (bmd ,lambda x ,y :(y [0 ],y [1 ]+x [1 ])))
        if bmd[-1] != (Decimal(1280), Decimal(0)):
            bmd.append((Decimal(1280), Decimal(0)))
        return bmd 

    def get_moment (self ,L ):
        locs ,moments =zip (*self .bmd ())
        idx =bisect_left (locs ,L )-1 
        # print(locs, L, idx)
        ans =lerp (locs [idx ],moments [idx ],locs [idx +1 ],moments [idx +1 ],L )
        # print('ANS', ans)
        return ans 
        # return moments[idx] + (moments[idx + 1] - moments[idx])*(L - self.locs[idx])/(self.locs[idx + 1] - self.locs[idx])

    def get_max_moment (self ,L0 ,L1 ):
        locs ,moments =zip (*self .bmd ())
        idx1 =bisect_left (locs ,L0 )-1 
        idx2 =bisect_left (locs ,L1 )-1 
        LL =[]
        if idx1 !=idx2 :
            LL =list (zip (moments [idx1 +1 :idx2 +1 ],locs [idx1 +1 :idx2 +1 ]))
        LL +=[(self .get_moment (L0 ),L0 ),(self .get_moment (L1 ),L1 ),(Decimal ('0'),-1 )]
        ex_min ,ex_max =min (LL ),max (LL )
        extras =list (filter (lambda v :(v [0 ]==ex_min [0 ])or (v [0 ]==ex_max [0 ]),LL ))
        return extras 

    def plot_bmd (self ):
        BMD_plot =self .bmd ()
        xdata ,ydata =zip (*BMD_plot )
        plt .plot (xdata ,ydata )
        plt .show ()

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
        location =bisect_left (locs ,L )-1 
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
            is_compressed =not ((cur_cross .starting [i ]>=cur_cross .ybar )^(self .get_moment (L )>Decimal ('0')))
            b ,l =cur_cross .bases [i ],cur_cross .lengths [i ]
            if not is_compressed :
                continue 
            if l >b :
            # print('CASE 3')
            # print('L', l)
            # print('B', b)
                crit =(b /Decimal ('2')/l )**Decimal ('2')*Decimal ('6')*Decimal(math .pi) **Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
                mcrit =min (mcrit ,crit )
                # print('POSS CRIT', crit)
            else :
                neighbours =list (filter (isbounded ,[i -1 ,i +1 ]))
                crit =None 
                poss =list (filter (lambda x :cur_cross .xdisp [x ]!=Decimal ('0'),neighbours ))
                if len (poss )>Decimal ('0'):
                # print('CASE 1 AND 2')
                    b2 =Decimal ('2')*cur_cross .xdisp [poss [0 ]]-cur_cross .bases [poss [0 ]]/Decimal ('2')#two held end
                    oneheld =(b -b2 -cur_cross .bases [poss [0 ]])/Decimal ('2')
                    twoheld =b2 
                    if oneheld <Decimal ('0')or twoheld <Decimal ('0'):
                        continue 
                    crit =Decimal ('0.425')*(l /(oneheld ))**Decimal ('2')*Decimal(math .pi)**Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
                    crit2 =Decimal ('4')*(l /(twoheld ))**Decimal ('2')*Decimal(math .pi)**Decimal ('2')*self .material .E /Decimal ('12')/(1 -self .material .u **Decimal ('2'))
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
        min_shear =math .inf 
        a =max (list (map (lambda i :self .diaphragms [i ]-self .diaphragms [i -1 ],range (1 ,len (self .diaphragms )))))
        # print('a', a)
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
            ax.add_patch (FancyArrowPatch ((loc ,-sign (force )*Decimal ('4')),(loc ,Decimal ('0')),mutation_scale =Decimal ('10')))
        ax .relim ()
        ax .autoscale_view ()
        plt .show ()
        pass 

    def get_max_stresses(self, x):
        maximumstresses = [self.material.max_compression_stress, self.material.max_tension_stress]
        v1 = self.plate_buckling(x)
        idx =int((sign(v1) + 1)//2)
        maximumstresses[idx] = sign(v1)*min(abs(v1), abs(maximumstresses[idx]))
        # print('MAX STRESSES', maximumstresses)
        return maximumstresses

    def failure_force_at_x (self, x):
        bottom = self.flex_stress(0, x)
        top = self.flex_stress(sum(self.cross.get_cross(x).lengths), x)
        stresses = [bottom, top]
        max_stresses = self.get_max_stresses(x)
        # print('MAX STRESSES----------', max_stresses)
        stresses = sorted(stresses)
        # print('STRESSES', stresses)
        max_glue_stress = list(map(lambda g: self.shear_stress(g, x), self.cross.get_cross(x).glue))
        max_shear_stress = [self.shear_stress(self.cross.get_cross(x).ybar, x)]
        # stresslimit = max(list(map(lambda x, y: math.inf if x == 0 else y/x, stresses, max_stresses)))
        # shearstresslimit = max(list(map(lambda x: math.inf if x == 0 else 2/x, max_glue_stress)) + list(map(lambda x: math.inf if x == 0 else min(self.material.max_shear_stress, self.shear_buckling(x)), max_shear_stress)))
        megalist = list(map(lambda x, y: math.inf if x == 0 else y/x, stresses, max_stresses)) + list(map(lambda x: math.inf if x == 0 else 2/x, max_glue_stress)) + list(map(lambda x: math.inf if x== 0 else min(self.material.max_shear_stress, self.shear_buckling(x))/x, max_shear_stress))
        megalist = list(map(abs, megalist))
        # print('MEGA LIST OF FORCES', megalist)
        return min(megalist)
        # return min(stresslimit, shearstresslimit)

    def failure_force (self):
        test_points = list(map(Decimal, np.linspace(0, float(self.L))))
        return min(list(map(self.failure_force_at_x, test_points)))
        # return 
P = 1
LBridge = 550 + 510 + 190
loadings = [[(550, -P), (LBridge, -P)]]
LTrain = (52 + 176 + 164 + 176/2)*2
# print(LTrain < LBridge)
train_loadings = np.linspace(0, LBridge - LTrain)
train_loadings = list(map(lambda x: [(52+x, -P), (52 + 176+x, -P), (52 + 176 + 164+x, -P), (52 + 176 + 164 + 176+x, -P), (52 + 176 + 164 + 176 + 164+x, -P), (52 + 176 + 164 + 176 + 164 + 176+x, -P)], train_loadings))
total_loadings = loadings + train_loadings

for loading in total_loadings:
    LBot = 50
    N = 3
    HEIGHT = 65
    # A = CrossSection([1.27, 1.27, HEIGHT - (3 + 2)*1.27, 1.27, 2*1.27], [LBot, 4*1.27 + 20, 4*1.27, 20 + 4*1.27, 100], [0, LBot/2 - 1.27 - 5 ,LBot/2-1.27,LBot/2-1.27 - 5, 0], [1.27, HEIGHT - 2*1.27])
    A = CrossSection([1.27, 75 - 1.27*3,1.27, 1.27], [80, 2*1.27, 20 + 2*1.27, 100], [0, 40-1.27/2, 35 - 1.27/2, 0], [75 - 1.27])
    # B = CrossSection([1.27*3, 1.27, HEIGHT - (3 + 3)*1.27, 1.27, 1.27], [100, 20 + 4*1.27, 4*1.27, 20 + 4*1.27, 100], [0, 50 - 1.27 - 5, 50 - 1.27, 50 - 1.27 - 5, 0], [2*1.27, HEIGHT - 1.27])
    # A.draw()
    # B.draw()
    # C = CrossGroup([A, B], [(0, 1), (LBridge, 1)], [1, 1], [0, 780])
    C = CrossGroup([A], [(0, 1), (LBridge, 1)], [1], [0])
    matboard = Material(4000, 30, -6, 4, 0.2)
    P = 1
    bridge_test = Bridge(C, LBridge, loading, [0, LBridge - 190], matboard,  [0, 550, 550 + 510, LBridge])
    # print('Bending Moment Peaks: ', bridge_test.bmd()[1:-1])
    # bridge_test.plot_sfd()
    # BM = [x[0] for x in bridge_test.bmd()[1:-1]]
    print(bridge_test.failure_force())