#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')

import mesa as ms
import pylab as pyl
import numpy as np
from astropy import units as u
from astropy import constants as const
import math
import copy


model_paths = ["/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/3.0/LOGS/", 
               "/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/4.0/LOGS/",
               "/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/5.0/LOGS/",
               "/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/7.0/LOGS/",]


h1=ms.history_data(model_paths[3])
h1.kippenhahn(1,"model")
pyl.show()

p1=ms.mesa_profile(model_paths[0],100,num_type='model')
mass = p1.get('mass') * u.Msun # solar mass
NN = mass.size
dm = copy.deepcopy(mass)
dm[:-1] = mass[:-1]-mass[1:]
print dm
print mass

r = p1.get('radius') * u.Rsun# solar radius
dr = copy.deepcopy(r)
dr[:-1] = r[:-1]-r[1:]
csound = p1.get('csound') * u.cm/u.s #cm/s
t_sc = dr/(csound)
mdot_sc =  dm/t_sc

t_ff = math.pi/2. * r**(3./2.) / np.sqrt(2.*const.G*(mass))
mdot_ff = dm/t_ff

fig = pyl.figure(2)

pyl.semilogy(mass.to('Msun').value,mdot_ff.to('Msun/yr').value,label="$local dm/dt_{free\ fall}$",linewidth=2)
pyl.semilogy(mass.to('Msun').value,mdot_sc.to('Msun/yr').value, label="$local dm/dt_{sound\ crossing}$",linewidth=2)
pyl.xlabel("Mass coordinate m ($M_{\odot}$)",fontsize=14)
pyl.ylabel("Mass-transfer rate ($M_{\odot}/yr$)",fontsize=14)

pyl.legend()

#Put figure window on top of all other windows
fig.canvas.manager.window.attributes('-topmost', 1)
# #After placing figure window on top, allow other windows to be on top of it later
fig.canvas.manager.window.attributes('-topmost', 0)

pyl.show()