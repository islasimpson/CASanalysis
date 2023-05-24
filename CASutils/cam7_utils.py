import matplotlib.pyplot as plt
import numpy as np
from dycoreutils import colormap_utils as mycolors
import sys

def plot_pre_mon(fig, data, pre, ci, cmin, cmax, expname, x1=None, x2=None, y1=None, y2=None, oplot=False, ax=None, cmap='precip', taxis='time', paxis='lev'):
    """
    Plot seasonal cycle, pressure versus time.
    """

    # move the time axis to the first
    if (data.dims[1] != taxis):
        data = data.transpose(..., taxis)

    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    if (cmap == "precip_nowhite"):
        mymap = mycolors.precip_cmap_nowhite(nlevs)


    # if overplotting, check for axis input
    if (oplot and (not ax)):
        print("This isn't going to work.  If overplotting, specify axis")
        sys.exit()

    plt.rcParams['font.size'] = '14'

    monticks_temp = np.arange(0,12,1)
    monticks2_temp = np.arange(0,12,1)+0.5

    monticks = monticks_temp
    monticks2 = np.zeros([len(monticks2_temp)+2])
    monticks2[0] = -0.5 ; monticks2[len(monticks2)-1] = 12.5
    monticks2[1:len(monticks2)-1] = monticks2_temp

    dataplot = np.zeros([data[paxis].size,len(monticks2)])
    dataplot[:,0] = data[:,11]
    dataplot[:,len(monticks2)-1] = data[:,0]
    dataplot[:,1:len(monticks2)-1] = data[:,:]


    if not oplot:
        if (x1):
            ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
        else:
            ax = fig.add_axes()

    ax.contourf(monticks2, -np.log10(data[paxis]), dataplot, levels=clevs, cmap=mymap, extend='max')
    ax.set_ylim(-np.log10(100),-np.log10(3))
    ax.set_yticks([-np.log10(100),-np.log10(30),-np.log10(10),-np.log10(3)])
    ax.set_yticklabels(['100','30','10','3'])
    ax.set_xlim([0,12])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2[1:13], minor=True)
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], minor=True, fontsize=14)
    ax.set_title(expname, fontsize=16)

    return ax

