import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def plothisto(fig, dat, bins, x1=0, x2=1, y1=0, y2=1, percent=False, xlim=None, ylim=None,
     yticks=None, yticklabels=None, ytitle=None, xticks=None, xticklabels=None, xtitle=None, 
     color='lightgray', label=None, addlines=False, title=None,orient='horizontal',alpha=1, 
     pointsonly=False, markercolor='black'):
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
    if (yticks):
        ax.set_yticks(yticks)
    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=14)
    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)
    if (xticks):
        ax.set_xticks(xticks)
    if (xticklabels):
        ax.set_xticklabels(xticklabels, fontsize=14)
    if (xtitle):
        ax.set_xlabel(xtitle, fontsize=14)
    if (title):
        ax.set_title(title, fontsize=16)


    if (percent):
        histo = (histo/dat.size)*100.

    if ( orient == 'horizontal' ):
        if (pointsonly == False):
            ax.bar(binedges[0:np.size(binedges)-1], histo, width=binedges[1]-binedges[0],bottom=0, 
                     edgecolor='black',color=color,label=label,alpha=alpha)
        else:
            ax.plot(binedges[0:np.size(binedges)-1]+(binedges[1]-binedges[0])/2., histo, "o", 
                markersize=8, color=markercolor)
    elif (orient == 'vertical'):
        ax.bar(0, binedges[1:np.size(binedges)] - binedges[0:np.size(binedges)-1], width=histo, 
         bottom = binedges[0:np.size(binedges)-1], edgecolor='black', color=color, label=label, align='edge',alpha=alpha)

    #!!! Note - doesn't work for the vertical orientation
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


def oplotbar(ax, x, y, bottom=0, color='lightgray', label='', hatch=False):

    if (hatch):
        ax.bar(x,y,width=1,bottom=bottom,edgecolor=color,facecolor='None',
           color=color, label=label, hatch='///')
    else:
        ax.bar(x,y,width=1,bottom=bottom,edgecolor='black',color=color,label=label)
    return ax


def calchisto(dat, bins, percent=False):
    histo, binedges = np.histogram(np.array(dat), bins)
    if (percent):
        histo = (histo/dat.size)*100.
    return histo, binedges

def compute_kde(dat, bsize, bmin, bmax):
    bins = np.arange(bmin, bmax, bsize)
    kernel = stats.gaussian_kde(dat)
    datpdf = kernel(bins)*100.*bsize
    datpdf = xr.DataArray(datpdf, dims=['bins'], coords=[bins])
    return datpdf
