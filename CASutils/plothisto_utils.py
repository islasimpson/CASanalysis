import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def plothisto(fig, dat, bins, x1=0, x2=1, y1=0, y2=1, percent=False, xlim=None, ylim=None,
     yticklabels=None, ytitle=None, xticklabels=None, xtitle=None, color='lightgray', label=None, addlines=False, title=None):
    """ Plot a histogram 
    Input: fig = the figure
           dat = data to plot
           bins = bins
           percent=False for number of values, True for percentage of values
    """
    histo, binedges = np.histogram(np.array(dat), bins)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])

    if (xlim):
        ax.set_xlim(xlim)
    if (ylim):
        ax.set_ylim(ylim)
    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=14)
    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)
    if (xticklabels):
        ax.set_xticklabels(xticklabels, fontsize=14)
    if (xtitle):
        ax.set_xlabel(xtitle, fontsize=14)
    if (title):
        ax.set_title(title, fontsize=16)


    if (percent):
        histo = (histo/dat.size)*100.

    ax.bar(binedges[0:np.size(binedges)-1], histo, width=binedges[1]-binedges[0],bottom=0, 
                     edgecolor='black',color=color,label=label)


    if (addlines):
        print("adding lines")
        binwidth=binedges[1]-binedges[0]
        for ibin in np.arange(0,np.size(binedges)-1,1):
            if (histo[ibin] > 1):
                for iy in np.arange(1,histo[ibin],1):
                    ax.plot([binedges[ibin]-binwidth/2., binedges[ibin]+binwidth/2.], [iy, iy], color='black', linewidth=1)

    return ax

def oplothisto(ax, dat, bins, percent=False, color='lightgray',label=None, alpha=1):
    """ Overplot a histogram
    input: ax = the figure axis
           dat = data to plot
           bins = bins
           percent = false for number of data points, true for percentage of data points
    """
    histo, binedges = np.histogram(np.array(dat), bins)
    
    if (percent):
        histo = (histo/dat.size)*100.
    
    ax.bar(binedges[0:np.size(binedges)-1], histo, width=binedges[1]-binedges[0], bottom=0,
                    edgecolor='black', color=color, label=label, alpha=alpha)

    return ax


def calchisto(dat, bins, percent=False):
    histo, binedges = np.histogram(np.array(dat), bins)
    if (percent):
        histo = (histo/dat.size)*100.
    return histo, binedges


