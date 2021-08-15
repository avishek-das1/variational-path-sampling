import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle,Wedge
import matplotlib.colors as mcol
import matplotlib.cm as cm
import math as m
from numba import jit
from multiprocessing import Pool
import multiprocessing
import matplotlib.pylab as pl
import random
import datetime
import scipy.special as sc
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy import signal
from scipy.interpolate import CubicSpline

#mpl.style.use("classic")
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams.update({'font.size': 28})
plt.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(figsize = (18, 4))
gs_top = fig.add_gridspec(1,2,top=0.95,hspace=0.2,wspace=0.2)
gs_bott = fig.add_gridspec(1,2,hspace=0.,wspace=0.55)

meta=np.loadtxt('meta.txt',skiprows=1)
omega=np.loadtxt('omega.txt',skiprows=1)
kvar=np.loadtxt('kvar.txt',skiprows=1)

ax0 = plt.subplot(gs_top[0,0])
temp=0
ax0.plot(np.arange(101),meta[:,4],lw=2,label=r'$\mathrm{Metadynamics}$')
ax0.plot(np.arange(1000)+100,omega[1:1001,4],lw=2,label=r'$\mathrm{MCVB-T}$')
ax0.plot(np.arange(1000)+1100,omega[1002:2002,4],lw=2,label=r'$\mathrm{MCVB}$')
ax0.hlines(-4.8035+np.log(0.2),0,2100,color='k',linestyle='--',lw=3)
ax0.set_xlim(0,2100)
ax0.legend(frameon='False',framealpha=0.0,loc='center right',bbox_to_anchor=(1.0, 0.3),handletextpad=0.5)
ax0.set_ylabel(r'$\ln\overline{h_{B|A}}-\overline{\Delta U_{\boldsymbol{\lambda}}}\bigr| _{B|A}$')
ax0.set_xlabel(r'$\mathrm{training~steps}$')

ax1 = plt.subplot(gs_top[0,1])
temp=0
for i in range(2,12):
    ax1.plot(np.arange(1000)+temp,omega[i*1001+1:(i+1)*1001,4],lw=2)
    temp+=1000
    ax1.vlines(temp,-19,-7,color='silver',alpha=0.3)
    ax1.hlines(kvar[i-2,0]+np.log(0.2),temp-1000,temp,color='k',linestyle=':',lw=3,zorder=10)

ax1.set_xlim(0,10000)
ax1.set_ylim(-18.5,-7.5)
ax1.set_xlabel(r'$\mathrm{training~steps}$')

ax0.annotate(r'$\mathrm{(a)}$',xy=(0.05, 0.94), xycoords='figure fraction')
ax0.annotate(r'$\mathrm{(b)}$',xy=(0.52, 0.94), xycoords='figure fraction')

plt.savefig('sm1.png',bbox_inches='tight',pad_inches=0.1)
#plt.show()

