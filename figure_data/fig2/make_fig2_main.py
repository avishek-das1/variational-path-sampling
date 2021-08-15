import numpy as np
import math as m
from numba import jit
from multiprocessing import Pool
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams.update({'font.size': 40})
plt.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']= r"\usepackage{amsmath}"

ls0 = 30
ms0 = 30
mew0=4.
c = 1.
height, width = 18,15
fig = plt.figure(figsize = (width, height))

c_tab = np.array(['firebrick', 'seagreen', 'steelblue'])

gs = fig.add_gridspec(8,8)
ax0 = plt.subplot(gs[0:4,0:4])

data = np.load('data_fig2a.npy', allow_pickle = True)
tsteps,dL,dLe  = data

data = np.load('data_fig2b.npy', allow_pickle = True)
imax_tabn,imax_tab, dU_cgf, dU1, dU2, dU3, dU4, dU1e, dU2e, dU3e, dU4e, lograte,lograte_err  = data

ax0.errorbar(tsteps,dL,\
                yerr=dLe,\
                fmt='o--', color = c_tab[0], alpha = 1.,lw=5,fillstyle='none',ms=ms0,mew=mew0,capsize=5)

h = 10.
FE = -h
lplus, lminus = np.ones(len(tsteps))*(lograte+lograte_err),np.ones(len(tsteps))*(lograte-lograte_err)
p6,=ax0.plot(tsteps,np.ones(len(tsteps))*lograte,'k--',lw=5,fillstyle='none',ms=ms0,mew=mew0 )
ax0.fill_between(np.arange(len(tsteps)), np.ones(len(tsteps))*(lograte-lograte_err),\
            np.ones(len(tsteps))*(lograte+lograte_err),facecolor="blue",color='blue', alpha=0.35)

ax0.set_ylabel(r'  $ \ln k t_\mathrm{f}$')
ax0.set_xlabel(r'  $ \mathrm{training~steps}$')
ax0.tick_params(axis='y', color  = 'k',labelsize=ls0)
ax0.tick_params(axis='x', color  = 'k',labelsize=ls0)

ax0.set_xlim([2,30])
ax0.set_ylim([-60,-3])

ax1 = plt.subplot(gs[0:4, 4:])

data = np.load('data_fig2b.npy', allow_pickle = True)
imax_tabn,imax_tab, dU_cgf, dU1, dU2, dU3, dU4, dU1e, dU2e, dU3e, dU4e, lograte,lograte_err  = data

p1=ax1.errorbar(np.concatenate((imax_tab[0:-1:4],imax_tab[-2:])),np.concatenate((dU1[0:-1:4],dU1[-2:])),yerr = np.concatenate((dU1e[0:-1:4],dU1e[-2:])),fmt='o--',color=c_tab[0],lw=5,fillstyle='none',ms=ms0,mew=mew0,alpha = 1.,capsize=5)[0]
p2=ax1.errorbar(np.concatenate((imax_tab[0:-1:4],imax_tab[-2:])),np.concatenate((dU2[0:-1:4],dU2[-2:])),yerr = np.concatenate((dU2e[0:-1:4],dU2e[-2:])),fmt='^--',color=c_tab[2],lw=5,fillstyle='none',ms=ms0,mew=mew0,alpha = 1.,capsize=5)[0]
# p3=ax1.errorbar(np.concatenate((imax_tab[0:-1:4],imax_tab[-2:])),np.concatenate((dU3[0:-1:4],dU3[-2:])),yerr = np.concatenate((dU3e[0:-1:4],dU3e[-2:])),fmt='X--',color=c_tab[1],label=r'  $\kappa_3$',lw=3,fillstyle='none',ms=15,mew=2,alpha = 1.,capsize=5)[0]
p4=ax1.errorbar(np.concatenate((imax_tab[0:-1:4],imax_tab[-2:])),np.concatenate((dU4[0:-1:4],dU4[-2:])),yerr = np.concatenate((dU4e[0:-1:4],dU4e[-2:])),fmt='s--',color='darkorange',lw=5,fillstyle='none',ms=ms0,mew=mew0,alpha = 1.,capsize=5)[0]
#
# p5,=ax1.plot(np.concatenate((imax_tabn[0:-1:3],imax_tab[-2:])), np.concatenate((dU_cgf[0:-1:3],dU_cgf[-2:])),"*--",color = 'k', alpha = 1.\
#                 ,label=r'  $ k_{\mathrm{exp}}$',lw=3,fillstyle='none',ms=15,mew=2,zorder=10)
lr_tab, lr_err_tab = lograte*np.ones(len(imax_tab)), lograte_err*np.ones(len(imax_tab))
imax_tab1 = imax_tab
imax_tab1[0]  = imax_tab1[0]-1
ax1.plot(imax_tab,lr_tab,'k--',lw=5,fillstyle='none',ms=ms0,mew=mew0 )

ax1.fill_between(imax_tab, np.ones(len(imax_tab))*(lograte-lograte_err),\
            np.ones(len(imax_tab))*(lograte+lograte_err),facecolor="blue",color='blue', alpha=0.15)
ax1.set_ylabel(r'  $ \ln k t_\mathrm{f}$')
ax1.set_xlabel(r'  $M_{R},M_{t}$')
ax1.tick_params(axis='y', color  = 'k',labelsize=ls0)
ax1.tick_params(axis='x', color  = 'k',labelsize=ls0)
ax1.set_xlim([2,40])
ax1.set_ylim([-15,-5])
ax1.set_yticks([-5,-10,-15])

ax3 = fig.add_subplot(gs[4:, :])

data = np.load('data_fig2c.npy', allow_pickle = True)
h_tab, lograte_tab, dS2,dS2e, dS4,dS4e,dS1,dS1e, dQ1, dQ1e = data


ax3.plot(h_tab,lograte_tab,'--',color = 'k',lw=5,fillstyle='none',ms=ms0,mew=mew0)
#p6,=ax3.plot(h_tab,-h_tab,'p--',color = 'goldenrod',label=r'$h$',lw=3,fillstyle='none',ms=15,mew=2,zorder=10)
ax3.errorbar(h_tab, dS1,\
                yerr=dS1e,\
                fmt='o', color = c_tab[0], alpha = 1.,capsize=5,lw=5,fillstyle='none',ms=ms0,mew=mew0)


ax3.errorbar(h_tab, dS2,\
                yerr=dS2e,\
                fmt='^', color = c_tab[2],alpha = 1, capsize=5,lw=5,fillstyle='none',ms=ms0,mew=mew0)
print('test')
# plt.legend([p1,p2,p4,p6],[r'$\kappa_{1}$',r'$\kappa_{2}$',r'$\kappa_{4}$',r'$\ln k t_\mathrm{f}$',],frameon='False',framealpha=0.0,loc='upper center',bbox_to_anchor=(0.5, 2.8),ncol=4,handletextpad=1,columnspacing=3.5)

ax3.set_xlabel(r'  $\Delta V/k_{\mathrm{B}}T$')
ax3.set_ylabel(r'  $\ln k t_\mathrm{f}$')
ax3.tick_params(axis='y', color  = 'k',labelsize=ls0)
ax3.tick_params(axis='x', color  = 'k',labelsize=ls0)


plt.subplots_adjust(hspace=18.0,wspace=15.0)

ax0.annotate(r'$\mathrm{(a)}$',xy=(0.02, 0.92), xycoords='figure fraction')
ax1.annotate(r'$\mathrm{(b)}$',xy=(0.51, 0.92), xycoords='figure fraction')
ax3.annotate(r'$\mathrm{(c)}$',xy=(.02, .45), xycoords='figure fraction')
plt.tight_layout()

plt.savefig('fig2.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()
