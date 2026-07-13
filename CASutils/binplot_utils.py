import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Arial')
from math import nan

from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm

from CASutils import colormap_utils as mycolors

def plot_binned_ai_prptile(fig, vp, aridity, ci, cmin, cmax, titlestr, x1, x2, y1, y2, signifdat=None, ylabel=False, xlabel=False, hatching=None, fsize=11, stipplesize=2):
    ax = fig.add_axes([x1, y1, (x2-x1), (y2-y1)])
    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.precip_cmap(nlevs)
    norm = BoundaryNorm(clevs, ncolors=mymap.N, clip=True)
    x = vp.pr + (vp.pr[1]-vp.pr[0])/2.
    y = vp.ai + (vp.ai[1]-vp.ai[0])/2.

    aridity0p05 = np.argmin(np.abs(np.array(aridity)-0.05))
    aridity0p5 = np.argmin(np.abs(np.array(aridity)-0.5))

#    ax.pcolormesh( (np.arange(0,vp.pr.size,1)/vp.pr.size)*100., (np.arange(0,vp.ai.size,1)/vp.ai.size)*100., vp, cmap=mymap, norm=norm)
    ax.pcolormesh( (x/len(x))*100., (y/len(y))*100., vp, cmap=mymap, norm=norm, edgecolor='none')
    if (hatching is not None):
        ax.fill_between([0,100],[0,0],[100,100],hatch='/////', color='none', edgecolor='white')
        ax.pcolormesh( (x/len(x))*100., (y/len(y))*100., hatching, cmap=mymap, norm=norm, edgecolor='none')
    ax.set_xlim(0,100)
    ax.set_ylim(0,100)

    ax.set_yticks([0,20,40,60,80,100])

    if (ylabel):
        ax.set_yticklabels(['0','20','40','60','80','100'], fontsize=fsize-2)
        ax.set_ylabel('Aridity Index percentile', fontsize=fsize, labelpad=0)
    else:
        ax.set_yticklabels([' ',' ',' ',' ',' ',' '])

    ax.set_xticks([0,50,100])
    ax.set_xticklabels(['0','50','100'], fontsize=fsize-2)
    if (xlabel):
        ax.set_xlabel('$\overrightarrow{pr}$ percentile', fontsize=fsize, labelpad=0)
    ax.set_title(titlestr, fontsize=fsize)

    ax.plot( (x/len(x))*100., x*0 + (y[aridity0p05]/len(y))*100., color='black', linestyle='dashed')
    ax.plot( (x/len(x))*100., x*0 + (y[aridity0p5]/len(y))*100., color='black', linestyle='dashed')

    if (signifdat is not None):
        signifdat = signifdat.stack(z=("ai","pr"))
        signifdat = signifdat.dropna("z")
        x = signifdat.pr + (vp.pr[1]-vp.pr[0])/2.
        y = signifdat.ai + (vp.ai[1]-vp.ai[0])/2.
        ax.scatter( (x/vp.pr.size)*100., (y/vp.ai.size)*100., s=stipplesize, color='black')
        #ax.scatter((signifdat.pr/vp.pr.size)*100., (signifdat.ai/vp.ai.size)*100., s=9, color='black')

    return ax


def plot_binned_ai_monthlyai(fig, vp, aridity, ci, cmin, cmax, titlestr, x1, x2, y1, y2, signifdat=None, ylabel=False, xlabel=False, hatching=None, cmap='precip', fsize=11, stipplesize=2):
    ax = fig.add_axes([x1, y1, (x2-x1), (y2-y1)])
    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    norm = BoundaryNorm(clevs, ncolors=mymap.N, clip=True)
    x = vp.month + (vp.month[1]-vp.month[0])/2.
    y = vp.ai + (vp.ai[1]-vp.ai[0])/2.

    aridity0p05 = np.argmin(np.abs(np.array(aridity)-0.05))
    aridity0p5 = np.argmin(np.abs(np.array(aridity)-0.5))

    ax.pcolormesh( x, (y/len(y))*100., vp, cmap=mymap, norm=norm, edgecolor='none')
    if (hatching is not None):
        ax.fill_between([0,12],[0,0],[100,100],hatch='/////', color='none', edgecolor='white')
        ax.pcolormesh( x, (y/len(y))*100., hatching, cmap=mymap, norm=norm, edgecolor='none')

    ax.set_yticks([0,20,40,60,80,100])

    if (ylabel):
        ax.set_yticklabels(['0','20','40','60','80','100'], fontsize=fsize-2)
        ax.set_ylabel('Aridity Index percentile', fontsize=fsize-1, labelpad=0)
    else:
        ax.set_yticklabels([' ',' ',' ',' ',' ',' '])

    ax.set_title(titlestr, fontsize=fsize)

    ax.plot( x, x*0 + (y[aridity0p05]/len(y))*100., color='black', linestyle='dashed')
    ax.plot( x, x*0 + (y[aridity0p5]/len(y))*100., color='black', linestyle='dashed')

    if (signifdat is not None):
        signifdat = signifdat.stack(z=("ai","month"))
        signifdat = signifdat.dropna("z")
        x = signifdat.month + (vp.month[1]-vp.month[0])/2.
        y = signifdat.ai + (vp.ai[1]-vp.ai[0])/2.
        ax.scatter( x, (y/vp.ai.size)*100., s=stipplesize, color='black')

    ax.set_xlim(0,12)
    ax.set_xticks([0,12])
    ax.set_xticklabels(['Most\n Arid','Most\n Humid'], fontsize=fsize-5)

    if (xlabel):
        ax.set_xlabel('Month', fontsize=fsize-1, labelpad=0)

    return ax





