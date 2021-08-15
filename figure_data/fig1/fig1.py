import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.colors as mcol
import matplotlib.cm as cm

#mpl.style.use("classic")
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams.update({'font.size': 28})
plt.rc('text', usetex=True)

fig = plt.figure(figsize=(10,4))
plt.subplots_adjust(hspace=0.0,wspace=0.4)
gs = gridspec.GridSpec(1, 2)
axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])]

cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["b","r"])
cpick = cm.ScalarMappable(cmap=cm1)
cpick.set_array([])
#axs[0].imshow([[0,1], [0,1]],cmap = cm1,interpolation = 'bicubic',extent=[1.42,1.72,20,23],aspect='auto')
axs[0].imshow([[0,1], [0,1]],cmap = cm1,interpolation = 'bicubic',extent=[1.37,1.67,19,22],aspect='auto')

#axs[0].annotate(r'$0$',xy=(1.4,16), xycoords='data')
#axs[0].annotate(r'$t_{f}$',xy=(1.7,16), xycoords='data')
axs[0].annotate(r'$0$',xy=(1.27,19.5), xycoords='data')
axs[0].annotate(r'$t_{f}$',xy=(1.7,19.5), xycoords='data')

ts=6
for i in range(ts):
    xV=np.loadtxt('fig1_a_t'+str(i*0.02)+'.txt',skiprows=1)
    axs[0].plot(xV[:,0],xV[:,1],color=(i/(ts-1),0,1-i/(ts-1)),lw=3)

traj=np.reshape(np.loadtxt('fig1_b_traj_100_200.txt',skiprows=1),(100,200))
for i in range(100):
    axs[1].plot((np.arange(200)+1)*1e-3/0.2,traj[i,:],'k',alpha=0.5)

axs[0].set_xlabel(r'$R/\sigma$',labelpad=-0.3)
axs[0].set_ylabel(r'$V_{t}(R)/k_{\mathrm{B}}T$')
axs[1].set_xlabel(r'$t/t_{f}$',labelpad=-0.3)
axs[1].set_ylabel(r'$R/\sigma$')
axs[1].set_xlim(0,1)
axs[1].set_ylim(0.9,2.3)
axs[0].set_xlim(0.9,2.3)
axs[0].set_ylim(-2,25)

axs[1].add_patch(Rectangle((0, 1.25), 1.0, -6,alpha=0.4,facecolor='lightpink'))
axs[1].add_patch(Rectangle((0, 1.85), 1.0, 6,alpha=0.4,facecolor='lightskyblue'))
axs[0].add_patch(Rectangle((1.25, -10), -6, 40,alpha=0.4,facecolor='lightpink'))
axs[0].add_patch(Rectangle((1.85, -10), 6, 40,alpha=0.4,facecolor='lightskyblue'))

axs[0].annotate(r'$\mathrm{(a)}$',xy=(0.02, 0.9), xycoords='figure fraction')
axs[1].annotate(r'$\mathrm{(b)}$',xy=(0.5, 0.9), xycoords='figure fraction')

plt.savefig('fig1.png',bbox_inches='tight',pad_inches=0.1)

#plt.show()
