#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')

import mesa as ms
import pylab as pyl
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy import units as u
from astropy import constants as const
import math
import copy
import matplotlib.gridspec as gridspec



model_paths = ["/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/3.0/LOGS/", 
               "/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/4.0/LOGS/",
               "/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/5.0/LOGS/",
               "/Users/tassos/Dropbox/Projects/CE/SINGLEgrid/7.0/LOGS/",]


profileN = [100, 1000, 1500, 2000, 2500]



fig = plt.figure(figsize=(8.267, 11.692))

# gridspec inside gridspec
gs = gridspec.GridSpec(5, 4, wspace=0.1, hspace=0.1)

for i in range(4):
	h=ms.history_data(model_paths[i])
	for j in range(5):
		ax = plt.subplot(gs[j,i])
		print i,j
		p=ms.mesa_profile(model_paths[i],profileN[j],num_type='nearest_model')
		mass = p.get('mass') * u.Msun # solar mass
		dm = copy.deepcopy(mass)
		dm[:-1] = mass[:-1]-mass[1:]

		r = p.get('radius') * u.Rsun# solar radius
		dr = copy.deepcopy(r)
		dr[:-1] = r[:-1]-r[1:]
		csound = p.get('csound') * u.cm/u.s #cm/s
		t_sc = dr/(csound)
		mdot_sc =  dm/t_sc
		t_ff = math.pi/2. * r**(3./2.) / np.sqrt(2.*const.G*(mass))
		mdot_ff = dm/t_ff

		plt.semilogy(mass.to('Msun').value,mdot_ff.to('Msun/yr').value,label="$local\ dm/dt_{free\ fall}$",linewidth=2)
		plt.semilogy(mass.to('Msun').value,mdot_sc.to('Msun/yr').value, label="$local\ dm/dt_{sound\ crossing}$",linewidth=2)
		if j == 0:
			plt.xlabel("Mass coord. m ($M_{\odot}$)",fontsize=12)
		else:
			plt.xlabel("")
			matplotlib.ticker.NullFormatter
		if i== 0:
			plt.ylabel("dM/dt ($M_{\odot}/yr$)",fontsize=12)
		else:
			plt.xlabel("")
			matplotlib.ticker.NullFormatter


		#plt.legend()

#Put figure window on top of all other windows
fig.canvas.manager.window.attributes('-topmost', 1)
# #After placing figure window on top, allow other windows to be on top of it later
fig.canvas.manager.window.attributes('-topmost', 0)

plt.show()






# h1=ms.history_data(model_paths[0])
# h1.kippenhahn(1,"model")
# pyl.show()

# h2=ms.history_data(model_paths[1])
# h2.kippenhahn(2,"model")
# pyl.show()

# h3=ms.history_data(model_paths[2])
# h3.kippenhahn(3,"model")
# pyl.show()

# h4=ms.history_data(model_paths[3])
# h4.kippenhahn(4,"model")
# pyl.show()

# fig = plt.figure(figsize=(20, 15))

# # gridspec inside gridspec
# gs = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.2)



# for i in range(4):
# 	ax = plt.subplot(gs[i])
# 	h = ms.history_data(model_paths[i])
# 	h.kippenhahn(-1,"model",symbol_size=4,print_legend=True)
# 	fig.add_subplot(ax)

# plt.show()


# p1=ms.mesa_profile(model_paths[0],100,num_type='model')
# mass = p1.get('mass') * u.Msun # solar mass
# NN = mass.size
# dm = copy.deepcopy(mass)
# dm[:-1] = mass[:-1]-mass[1:]
# print dm
# print mass

# r = p1.get('radius') * u.Rsun# solar radius
# dr = copy.deepcopy(r)
# dr[:-1] = r[:-1]-r[1:]
# csound = p1.get('csound') * u.cm/u.s #cm/s
# t_sc = dr/(csound)
# mdot_sc =  dm/t_sc

# t_ff = math.pi/2. * r**(3./2.) / np.sqrt(2.*const.G*(mass))
# mdot_ff = dm/t_ff

# fig = pyl.figure(2)

# pyl.semilogy(mass.to('Msun').value,mdot_ff.to('Msun/yr').value,label="$local dm/dt_{free\ fall}$",linewidth=2)
# pyl.semilogy(mass.to('Msun').value,mdot_sc.to('Msun/yr').value, label="$local dm/dt_{sound\ crossing}$",linewidth=2)
# pyl.xlabel("Mass coordinate m ($M_{\odot}$)",fontsize=14)
# pyl.ylabel("Mass-transfer rate ($M_{\odot}/yr$)",fontsize=14)

# pyl.legend()

# #Put figure window on top of all other windows
# fig.canvas.manager.window.attributes('-topmost', 1)
# # #After placing figure window on top, allow other windows to be on top of it later
# fig.canvas.manager.window.attributes('-topmost', 0)

# pyl.show()


