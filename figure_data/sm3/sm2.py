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
gs_top = fig.add_gridspec(1,3,top=0.95,hspace=0.2,wspace=0.4)

M=200
xlow=-40
xhigh=70
ax0 = plt.subplot(gs_top[0,0])

td=np.loadtxt('tilted_distribution.txt',skiprows=1)
ax0.plot(td[:,0],td[:,1],color='r',marker='o',markersize=15,mew=2,fillstyle='none',linestyle='none',label=r'$e^{-\Delta U_{\boldsymbol{\lambda}}}P_{B|A,\boldsymbol{\lambda}}$')
ud=np.loadtxt('unbiased_distribution.txt',skiprows=1)
ax0.plot(ud[:,0],ud[:,1],color='b',marker='s',markersize=15,mew=2,fillstyle='none',linestyle='none',label=r'$k_{\mathrm{rwt}}t_{f}P_{B|A}$')
ax0.set_yscale('log')
ax0.set_ylim(1e-8,1e-2)
ax0.set_xlim(-20,20)
ax0.set_xlabel(r'$\Delta U_{\boldsymbol{\lambda}}$')
ax0.set_ylabel(r'$P(\Delta U_{\boldsymbol{\lambda}})$')
ax0.legend(frameon='False',framealpha=0.0,loc='lower left',bbox_to_anchor=(-0.14, -0.1),handletextpad=0.0)

ax1= plt.subplot(gs_top[0,1])
k=-4.493094264579371
kexp=np.loadtxt('kexp_Nw.txt',skiprows=1)
krwt=np.loadtxt('krwt_Nw.txt',skiprows=1)
ax1.errorbar(kexp[:,0],kexp[:,1]-k,yerr=kexp[:,2],markersize=15,mew=2,lw=3,fillstyle='none',color='g',marker='s',capsize=5,label=r'$k_{\mathrm{exp}}$')
ax1.errorbar(krwt[:,0],krwt[:,1]-k,yerr=krwt[:,2],markersize=15,mew=2,lw=3,fillstyle='none',color='b',marker='o',capsize=5,label=r'$k_{\mathrm{rwt}}$')
ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.set_xticks([1e4,1e5,1e6])
ax1.set_xlabel(r'$N_{w}$')
ax1.set_ylabel(r'$\mathrm{Error~in~}\ln kt_{f}$',labelpad=-0.2)
ax1.legend(frameon='False',framealpha=0.0,loc='center right',bbox_to_anchor=(1.05, 0.83),handletextpad=0.3)
ax1.set_ylim(-1,1)
ax1.set_yticks([-1,0,1])

ax2= plt.subplot(gs_top[0,2])
bs=np.loadtxt('biased_2000traj_tp_trxn.txt',skiprows=1)
ubs=np.loadtxt('unbiased_2000traj_tp_trxn.txt',skiprows=1)
M=20
ts=np.linspace(0,0.2,M+1)
tts=0.5*(ts[:-1]+ts[1:])
bs1h,temp=np.histogram(bs[:,0],bins=ts,density=True)
bs2h,temp=np.histogram(bs[:,1],bins=ts,density=True)
ubs1h,temp=np.histogram(ubs[:,0],bins=ts,density=True)
ubs2h,temp=np.histogram(ubs[:,1],bins=ts,density=True)
ax2.plot(tts/0.2,bs1h,lw=3,color='r',marker='None',label=r'$P_{B|A,\boldsymbol{\lambda}}(\tau^{\ddag})$')
ax2.plot(tts/0.2,ubs1h,lw=3,color='b',marker='None',label=r'$P_{B|A}(\tau^{\ddag})$')
ax2.plot(tts/0.2,bs2h,lw=3,color='r',linestyle='--',marker='None',label=r'$P_{B|A,\boldsymbol{\lambda}}(t_{\mathrm{rxn}})$')
ax2.plot(tts/0.2,ubs2h,lw=3,color='b',linestyle='--',marker='None',label=r'$P_{B|A}(t_{\mathrm{rxn}})$')
ax2.set_xlabel(r'$t/t_{f}$')
ax2.set_ylabel(r'$\mathrm{Transition~path~PDF}$')
ax2.legend(frameon='False',framealpha=0.0,loc='center right',bbox_to_anchor=(1.87, 0.5),handletextpad=0.3)

ax0.annotate(r'$\mathrm{(a)}$',xy=(0.025, 0.935), xycoords='figure fraction')
ax0.annotate(r'$\mathrm{(b)}$',xy=(0.31, 0.935), xycoords='figure fraction')
ax0.annotate(r'$\mathrm{(c)}$',xy=(0.575, 0.935), xycoords='figure fraction')

plt.savefig('sm2.png',bbox_inches='tight',pad_inches=0.1)
#plt.show()

