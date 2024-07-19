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

fig = plt.figure(figsize = (10, 12))
gs_top = fig.add_gridspec(4,4,top=0.95,bottom=-0.05,hspace=0.7,wspace=0.55)
#gs_bott = fig.add_gridspec(4,4,hspace=0.,wspace=0.7)

x0=np.loadtxt('snapshot_start.txt',skiprows=1)
xt=np.loadtxt('snapshot_end.txt',skiprows=1)
L=(80.0/0.6)**0.5
wca=0.5#2**(1/6)/2
f=0.7*wca #half of length of arrow

ax0 = plt.subplot(gs_top[2,0:2])
ax0.set_xlim(-3.7,5.0)
ax0.set_ylim(-2.3,2.3)
for p in range(1,2):
    for q in range(1,2):
        plt.gca().add_patch(plt.Circle((x0[0,0]+(p-1)*L, x0[0,1]+(q-1)*L), wca, color='r', clip_on=True))
        plt.gca().add_patch(plt.Circle((x0[1,0]+(p-1)*L, x0[1,1]+(q-1)*L), wca, color='r', clip_on=True))
        for i in range(2,80):
            plt.gca().add_patch(plt.Circle((x0[i,0]+(p-1)*L, x0[i,1]+(q-1)*L), wca, color='b', fill=False,clip_on=True))
            plt.arrow(x0[i,0]+(p-1)*L-f*np.cos(x0[i,2]),x0[i,1]+(q-1)*L-f*np.sin(x0[i,2]),2*f*np.cos(x0[i,2]),2*f*np.sin(x0[i,2]),color='b',head_width=0.3,length_includes_head=True)
ax0.set_xticks([])
ax0.set_xticks([], minor=True)
ax0.set_yticks([])
ax0.set_yticks([], minor=True)

ax1 = plt.subplot(gs_top[2,2:])
ax1.set_xlim(-3.7,5.0)
ax1.set_ylim(-2.3,2.3)
for p in range(1,2):
    for q in range(1,2):
        plt.gca().add_patch(plt.Circle((xt[0,0]+(p-1)*L, xt[0,1]+(q-1)*L), wca, color='r', clip_on=True))
        plt.gca().add_patch(plt.Circle((xt[1,0]+(p-1)*L, xt[1,1]+(q-1)*L), wca, color='r', clip_on=True))
        for i in range(2,80):
            plt.gca().add_patch(plt.Circle((xt[i,0]+(p-1)*L, xt[i,1]+(q-1)*L), wca, color='b', fill=False,clip_on=True))
            plt.arrow(xt[i,0]+(p-1)*L-f*np.cos(xt[i,2]),xt[i,1]+(q-1)*L-f*np.sin(xt[i,2]),2*f*np.cos(xt[i,2]),2*f*np.sin(xt[i,2]),color='b',head_width=0.3,length_includes_head=True)
ax1.set_xticks([])
ax1.set_xticks([], minor=True)
ax1.set_yticks([])
ax1.set_yticks([], minor=True)

bf=np.loadtxt('bruteforceestimate.txt',skiprows=1)
es=np.loadtxt('expestimate.txt',skiprows=1)
mb=np.loadtxt('reweighting.txt',skiprows=1)
q=np.loadtxt('heat_corrected.txt',skiprows=1)
#qn=np.loadtxt('list_unconditioned.txt')
#qn=np.loadtxt('list_stratonovich.txt')
#qn=np.loadtxt('list.txt')
#qn[:,1:]/=20

ax2 = plt.subplot(gs_top[0:2,:])
ax2.errorbar(bf[:,0],bf[:,1],yerr=bf[:,2],color='k',marker='x',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$k$',capsize=5)
ax2.errorbar(es[:,0],es[:,1],yerr=es[:,2],color='g',marker='s',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$k_{\mathrm{exp}}$',capsize=5)
ax2.errorbar(mb[:,0],mb[:,1],yerr=mb[:,2],color='b',marker='o',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$k_{\mathrm{rwt}}$',capsize=5)
#ax2.errorbar(q[:,0],q[:,1],yerr=q[:,2],color='r',marker='^',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$\frac{Q_{\mathrm{sol}}}{2k_{\mathrm{B}}T^{\prime}}$',capsize=5)

#ax2.fill_between(q[:,0],q[:,1],np.zeros(10)+150,color='r',alpha=0.2)

ax2.hlines(-4.36+np.log(0.2),0,1.0,color='k',linestyle='-',lw=7,zorder=0,alpha=0.4,label=r'$k_{\mathrm{iso}}$')
ax2.set_xlabel(r'$v_{0}\sigma/k_{\mathrm{B}}T$')
ax2.set_ylabel(r'$\ln kt_{f}$')
ax2.set_ylim(-7.5,-1.1)
ax2.set_xlim(0,18)
ax2.set_xticks([0,6,12,18])

axins = inset_axes(ax2, width='50%',height='50%',bbox_to_anchor=(-110,650,750,500))
#axins = inset_axes(ax2, width='50%',height='40%')
 
axins.errorbar(bf[:,0],bf[:,1]-bf[0,1],yerr=bf[:,2]+bf[0,2],color='k',marker='x',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$k$',capsize=5)
axins.errorbar(es[:,0],es[:,1]-bf[0,1],yerr=es[:,2]+bf[0,2],color='g',marker='s',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$k_{\mathrm{exp}}$',capsize=5)
axins.errorbar(mb[:,0],mb[:,1]-bf[0,1],yerr=mb[:,2]+bf[0,2],color='b',marker='o',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$k_{\mathrm{rwt}}$',capsize=5)
axins.errorbar(q[:,0],(q[:,1]-q[0,1]),yerr=(q[:,2]+q[0,2]),color='r',marker='^',lw=3,fillstyle='none',markersize=15,mew=2,linestyle='none',label=r'$\frac{Q_{\mathrm{sol}}}{2k_{\mathrm{B}}T^{\prime}}$',capsize=5)

axins.set_xlabel(r'$v_{0}\sigma/k_{\mathrm{B}}T$',labelpad=-20)
axins.set_ylabel(r'$\ln k/k_{0}$',labelpad=-10)
axins.set_ylim(1e-2,q[-1,1]-q[0,1]+100)
axins.set_xlim(1,18)
axins.set_xticks([2,18])
axins.set_yscale('log')
axins.fill_between(q[:,0],q[:,1]-q[0,1],np.zeros(10)+210,color='r',alpha=0.2)

ax0.annotate(r'$\mathrm{(a)}$',xy=(0.03, 0.94), xycoords='figure fraction')
ax0.annotate(r'$\mathrm{(b)}$',xy=(0.03, 0.22), xycoords='figure fraction')
ax0.annotate(r'$\mathrm{(c)}$',xy=(0.53, 0.22), xycoords='figure fraction')

plt.savefig('fig3_corrected.png',bbox_inches='tight',pad_inches=0.1)
#plt.show()

