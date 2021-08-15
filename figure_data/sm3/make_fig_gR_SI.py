import numpy as np
import matplotlib.pyplot as plt
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
plt.rc('text', usetex=True)


def compare_gR_v_ABTS():
    c = 4
    width, height = c*2, c*2
    fig = plt.figure(figsize = (width, height))
    gs = fig.add_gridspec(3,2)

    ymax   = 2.25

    vo_tab = [0,3,6,9]
    h_tab = [3.5]
    h = 3.5

    # A,B
    ax0 = fig.add_subplot(gs[0,0])
    plt.tight_layout()
    vo =  vo_tab[3]
    fname = 'gR_A_cut_2h' + str(int(2*h)) + '_' + str(int(vo))+'.npy'
    data = np.load(fname, allow_pickle = True)
    gR1, gRe1, rtab, thetatab, rmax = data
    vo =  vo_tab[0]
    fname = 'gR_A_cut_2h' + str(int(2*h)) + '_' + str(int(vo))+'.npy'
    data = np.load(fname, allow_pickle = True)
    gR2, gRe2, rtab, thetatab, rmax = data


    tind = np.argmin(np.absolute(thetatab-m.pi/2.))
    zmin, zmax = -3.,22.
    zmin2, zmax2 = -.8,2.2

    print('zmin, zmax = ',zmin, zmax)
    levels = MaxNLocator(nbins=70).tick_values(zmin, zmax)
    levels2 = MaxNLocator(nbins=70).tick_values(zmin2, zmax2)
    cmap = plt.get_cmap('hsv')
    norm = BoundaryNorm(levels2, ncolors=cmap.N, clip=True)
    x = thetatab
    y = rtab

    z = gR1-gR2
    cf0 = ax0.contourf(x,y, z, levels=levels,cmap=cmap)
    ax0.set_ylim([0,ymax])
    ax0.set_xticks([0,np.pi/2,np.pi])
    ax0.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$'])
    ax0.set_yticks([0,1,2])
    ax0.set_yticklabels([r'$0$',r'$1$',r'$2$'])
    plt.xlabel(r'  $ \phi_1 $', {'color': 'k', 'fontsize': 20})

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.legend( bbox_to_anchor=(.45,.15,.5, .502),
              mode="expand", borderaxespad=0.,fontsize=20,frameon = False)

    plt.ylabel(r'  $ r/\sigma $', {'color': 'k', 'fontsize': 20})

    cbaxes = ax0.inset_axes( [0.3, 0.85, 0.3, 0.03])
    cbar = plt.colorbar(cf0,cax=cbaxes, ticks = [zmin, zmax],orientation='horizontal')
    cbaxes.tick_params(labelsize=16)


    ax1 = fig.add_subplot(gs[0,1])


    z = gRe1-gRe2
    cf1 = ax1.contourf(x,y, z, levels=levels2,cmap=cmap)
    ax1.set_ylim([0,ymax])
    ax1.set_xticks([0,np.pi/2,np.pi])
    ax1.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$'])
    ax1.set_yticks([0,1,2])
    ax1.set_yticklabels([r'$0$',r'$1$',r'$2$'])
    plt.xlabel(r'  $ \phi_2 $', {'color': 'k', 'fontsize': 20})

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.legend( bbox_to_anchor=(.45,.15,.5, .502),
              mode="expand", borderaxespad=0.,fontsize=20,frameon = False)

    plt.ylabel(r'  $ r/\sigma $', {'color': 'k', 'fontsize': 20})

    cbaxes = ax1.inset_axes( [0.3, 0.85, 0.3, 0.03])
    cbar = plt.colorbar(cf1,cax=cbaxes, ticks = [zmin2, zmax2],orientation='horizontal')
    cbaxes.tick_params(labelsize=16)




    # C,D
    ax2 = fig.add_subplot(gs[1,0])
    # plt.tight_layout()
    vo =  vo_tab[3]
    fname = 'gR_tpe_cut_2h' + str(int(2*h)) + '_' + str(int(vo))+'.npy'
    data = np.load(fname, allow_pickle = True)
    gR1, gRe1, rtab, thetatab, rmax = data
    vo =  vo_tab[0]
    fname = 'gR_tpe_cut_2h' + str(int(2*h)) + '_' + str(int(vo))+'.npy'
    data = np.load(fname, allow_pickle = True)
    gR2, gRe2, rtab, thetatab, rmax = data


    tind = np.argmin(np.absolute(thetatab-m.pi/2.))
    zmin, zmax = -3.,22.
    zmin2, zmax2 = -.8,2.2

    print('zmin, zmax = ',zmin, zmax)
    levels = MaxNLocator(nbins=70).tick_values(zmin, zmax)
    levels2 = MaxNLocator(nbins=70).tick_values(zmin2, zmax2)
    cmap = plt.get_cmap('hsv')
    norm = BoundaryNorm(levels2, ncolors=cmap.N, clip=True)
    x = thetatab
    y = rtab

    z = gR1-gR2
    cf2 = ax2.contourf(x,y, z, levels=levels,cmap=cmap)
    ax2.set_ylim([0,ymax])
    ax2.set_xticks([0,np.pi/2,np.pi])
    ax2.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$'])
    ax2.set_yticks([0,1,2])
    ax2.set_yticklabels([r'$0$',r'$1$',r'$2$'])
    plt.xlabel(r'  $ \phi_1 $', {'color': 'k', 'fontsize': 20})

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.legend( bbox_to_anchor=(.45,.15,.5, .502),
              mode="expand", borderaxespad=0.,fontsize=20,frameon = False)

    plt.ylabel(r'  $ r/\sigma $', {'color': 'k', 'fontsize': 20})

    cbaxes = ax2.inset_axes( [0.3, 0.85, 0.3, 0.03])
    cbar = plt.colorbar(cf2,cax=cbaxes, ticks = [zmin, zmax],orientation='horizontal')
    cbaxes.tick_params(labelsize=16)


    ax3 = fig.add_subplot(gs[1,1])


    z = gRe1-gRe2
    cf3 = ax3.contourf(x,y, z, levels=levels2,cmap=cmap)
    ax3.set_ylim([0,ymax])
    ax3.set_xticks([0,np.pi/2,np.pi])
    ax3.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$'])
    ax3.set_yticks([0,1,2])
    ax3.set_yticklabels([r'$0$',r'$1$',r'$2$'])
    plt.xlabel(r'  $ \phi_2 $', {'color': 'k', 'fontsize': 20})


    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.legend( bbox_to_anchor=(.45,.15,.5, .502),
              mode="expand", borderaxespad=0.,fontsize=20,frameon = False)

    plt.ylabel(r'  $ r/\sigma $', {'color': 'k', 'fontsize': 20})

    cbaxes = ax3.inset_axes( [0.3, 0.85, 0.3, 0.03])
    cbar = plt.colorbar(cf3,cax=cbaxes, ticks = [zmin2, zmax2],orientation='horizontal')
    cbaxes.tick_params(labelsize=16)




    # E,F
    ax4 = fig.add_subplot(gs[2,0])
    vo =  vo_tab[3]
    fname = 'gR_B_cut_2h' + str(int(2*h)) + '_' + str(int(vo))+'.npy'
    data = np.load(fname, allow_pickle = True)
    gR1, gRe1, rtab, thetatab, rmax = data
    vo =  vo_tab[0]
    fname = 'gR_B_cut_2h' + str(int(2*h)) + '_' + str(int(vo))+'.npy'
    data = np.load(fname, allow_pickle = True)
    gR2, gRe2, rtab, thetatab, rmax = data


    tind = np.argmin(np.absolute(thetatab-m.pi/2.))
    zmin, zmax = -3.,22.
    zmin2, zmax2 = -.8,2.2

    print('zmin, zmax = ',zmin, zmax)
    levels = MaxNLocator(nbins=70).tick_values(zmin, zmax)
    levels2 = MaxNLocator(nbins=70).tick_values(zmin2, zmax2)
    cmap = plt.get_cmap('hsv')
    norm = BoundaryNorm(levels2, ncolors=cmap.N, clip=True)
    x = thetatab
    y = rtab

    z = gR1-gR2
    cf4 = ax4.contourf(x,y, z, levels=levels,cmap=cmap)
    ax4.set_ylim([0,ymax])
    ax4.set_xticks([0,np.pi/2,np.pi])
    ax4.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$'])
    ax4.set_yticks([0,1,2])
    ax4.set_yticklabels([r'$0$',r'$1$',r'$2$'])
    plt.xlabel(r'  $ \phi_1 $', {'color': 'k', 'fontsize': 20})

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.legend( bbox_to_anchor=(.45,.15,.5, .502),
              mode="expand", borderaxespad=0.,fontsize=20,frameon = False)

    plt.ylabel(r'  $ r/\sigma $', {'color': 'k', 'fontsize': 20})

    cbaxes = ax4.inset_axes( [0.3, 0.85, 0.3, 0.03])
    cbar = plt.colorbar(cf4,cax=cbaxes, ticks = [zmin, zmax],orientation='horizontal')
    cbaxes.tick_params(labelsize=16)


    ax5 = fig.add_subplot(gs[2,1])


    z = gRe1-gRe2
    cf5 = ax5.contourf(x,
                      y, z, levels=levels2,
                      cmap=cmap)
    ax5.set_ylim([0,ymax])
    ax5.set_xticks([0,np.pi/2,np.pi])
    ax5.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$'])
    ax5.set_yticks([0,1,2])
    ax5.set_yticklabels([r'$0$',r'$1$',r'$2$'])
    plt.xlabel(r'  $ \phi_2 $', {'color': 'k', 'fontsize': 20})

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.legend( bbox_to_anchor=(.45,.15,.5, .502),
              mode="expand", borderaxespad=0.,fontsize=20,frameon = False)

    plt.ylabel(r'  $ r/\sigma $', {'color': 'k', 'fontsize': 20})

    cbaxes = ax5.inset_axes( [0.3, 0.85, 0.3, 0.03])
    cbar = plt.colorbar(cf5,cax=cbaxes, ticks = [zmin2, zmax2],orientation='horizontal')
    cbaxes.tick_params(labelsize=16)

    plt.tight_layout()

    plt.savefig('SI_fig_gR.pdf',bbox_inches='tight',pad_inches=0.1)






    plt.show()

compare_gR_v_ABTS()
