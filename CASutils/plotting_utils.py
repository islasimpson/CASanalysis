import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from CASutils import colormap_utils as mycolors


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



