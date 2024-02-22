import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from CASutils import colormap_utils as mycolors
import matplotlib.colors as colors
import warnings
warnings.filterwarnings('ignore')

def plotbox(fig, x1,x2,y1,y2, color='red', alpha=1):
    """plot a colored box"""
    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    box = Rectangle((0.,0.),1,1,facecolor=color, alpha=alpha)
    ax.axis('off')
    ax.add_patch(box)
    return ax


def plot2dhisto(fig, dat1, dat2, xvals, yvals, ci, cmin, cmax, x1, x2, y1, y2, titlestr, xticks=None, yticks=None, xticknames=None, yticknames=None, xlabel=None, ylabel=None):
    """plot a 2D histogram"""
    
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.blue2red_cmap(nlevs)
    
    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    
    ax.hist2d(dat1, dat2, bins=[xvals, yvals], cmap = mymap, vmin=cmin, vmax=cmax)
    ax.set_title(titlestr, fontsize=16)    

    if (xticks):
        ax.set_xticks(xticks)
    if (xticknames):
        ax.set_xticklabels(xticknames, fontsize=14)
    if (yticks):
        ax.set_yticks(yticks)
    if (yticknames):
        ax.set_yticklabels(yticknames, fontsize=14)
    if (xlabel):
        ax.set_xlabel(xlabel, fontsize=16)
    if (ylabel):
        ax.set_ylabel(ylabel, fontsize=16)


    return ax


def plotWK(fig, dat, x, y, ci, cmin, cmax, titlestr, x1, x2, y1, y2, cmap="wk", xlim=[-15,15], ylim=[0,0.5], xlabel=True, ylabel=True, xticks=[-10,0,10], xticklabels=['-10','0','10'], yticks=[0,0.1,0.2,0.3,0.4,0.5],yticklabels=['0','0.1','0.2','0.3','0.4','0.5'], contourlines=True,contourlinescale=1, speclevs=None, posonly=False):
    """ ???? """

    if (speclevs is None):
        nlevs = (cmax - cmin)/ci + 1
        clevs = np.arange(cmin, cmax+ci, ci)
    else:
        clevs = speclevs
        nlevs = np.size(clevs)

    import importlib
    importlib.reload(mycolors)


    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs, posonly=posonly)
    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)
    if (cmap == "wk"):
        mymap = mycolors.wkcmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if (ylabel):
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=14)
        ax.set_ylabel('Frequency (day$^{-1}$)', fontsize=14)
    else:
        ax.set_yticks(yticks)
        ax.set_yticklabels([])

    if (xlabel):
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=14)
        ax.set_xlabel('Wavenumber', fontsize=14)
    else:
        ax.set_xticks(xticks)
        ax.set_xticklabels([])

    ax.set_title(titlestr, fontsize=16)

    if (speclevs is not None):
        norm = colors.BoundaryNorm(boundaries=clevs, ncolors=256)
        ax.contourf(x,y,dat,levels=clevs,cmap=mymap,extend='both', norm=norm)
    else:
        ax.contourf(x, y, dat, levels=clevs, cmap = mymap, extend='both')

    if (contourlines):
        if (speclevs is None):
            clevlines = clevs*contourlinescale
            clevlines = clevlines[np.abs(clevlines) > (ci/2.)]
        else:
            clevlines = clevs
        if (posonly):
            ax.contour(x, y, dat, levels=clevlines[1:len(clevlines)+1], colors='black', linestyles='solid')
        else:
            ax.contour(x, y, dat, levels=clevlines, colors='black')


    kmidpos = np.min(dat.k.values[dat.k > 0])
    kmidneg = np.max(dat.k.values[dat.k < 0])
    ax.fill_between([kmidneg,kmidpos],[ylim[0],ylim[0]],[ylim[1],ylim[1]], color='lightgray')

    return ax

