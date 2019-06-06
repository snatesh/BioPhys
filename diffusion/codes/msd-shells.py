import os
import numpy as np
from numpy import *
from scipy import loadtxt, optimize
import matplotlib.pyplot as plt
from matplotlib.ticker import *
import pylab

# 10ps 100ps 500ps 1000ps 5000ps 10000ps 20000ps 100000ps
wind_sizes = [3,4,5]
water_wind_sizes = [2*wind_size for wind_size in wind_sizes]
ps_per_frame = 2

#axes = plt.axes()
#first = 0
#fig,axes = plt.subplots(1,len(wind_sizes))
thickness = 2
rad = 10;
nshells = int(rad/thickness)
r = thickness*np.arange(1,nshells+1)
for w in wind_sizes:
  t = np.linspace(0,w*ps_per_frame,w)
  tw = np.linspace(0,w*ps_per_frame,2*w)
  x_diffs = np.zeros(nshells)
  y_diffs = np.zeros(nshells)
  z_diffs = np.zeros(nshells)
  for nshell in range(1,nshells+1):
    folder = "shell"+str(nshell)+"/tau"+str(w)+"f/"
    #t = np.linspace(first*ps_per_frame,w*ps_per_frame,w-first)
    if os.path.exists(folder):
      x_diff = loadtxt(folder+'unwrap-diff-x_out.txt', unpack=True,skiprows=0,usecols=[0])
      y_diff = loadtxt(folder+'unwrap-diff-y_out.txt', unpack=True,skiprows=0,usecols=[0])
      z_diff = loadtxt(folder+'unwrap-diff-z_out.txt', unpack=True,skiprows=0,usecols=[0])
      x_diffs[nshell-1] = x_diff[-1]
      y_diffs[nshell-1] = y_diff[-1]
      z_diffs[nshell-1] = z_diff[-1]
  plt.rcParams["figure.figsize"] = [16,9]
  axes = plt.axes()
  axes.plot(r,x_diffs,'s-',linewidth=1,label='x')	
  axes.plot(r,y_diffs,'s-',linewidth=1,label='y')	
  axes.plot(r,z_diffs,'s-',linewidth=1,label='z')
  axes.legend(loc='lower center',prop={'size':17}, fancybox=False, shadow=False,frameon=False,ncol=3)
  axes.set_title('2M4J pore water MSD as distance to pore surface varies',fontsize=25)
  axes.set_xlabel('Distance to pore surface  ($\\AA$)',fontsize=25)
  axes.set_ylabel('$\langle MSD(\\tau=$' + str(w*ps_per_frame) + ' ps$)$' + '$\\rangle$ ($\\AA^2$)',fontsize=25)
  plt.show()
  #plt.savefig("2m4j_msd_tau"+str(w*ps_per_frame)+"ps_compGap.png")
  axes.clear()
  #axes.plot(t,x_diff[first:w],'o-',linewidth=1.7,label='x, ' + 'tau = ' + str(w*ps_per_frame))
  #axes.plot(t,y_diff[first:w],'^-',linewidth=1.7,label='y, ' + 'tau = ' + str(w*ps_per_frame))
  #axes.plot(t,z_diff[first:w],'o-',linewidth=1.7,label='z, ' + 'tau = ' + str(w*ps_per_frame))
  #first = w

	


#x2 = np.arange(0,40,.001)
