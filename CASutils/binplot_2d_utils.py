import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan

from CASutils import colormap_utils as mycolors
from CASutils import linfit_utils as linfit
from matplotlib.colors import BoundaryNorm
import sys


def meanw(dat,axis):
    meanw = (dat*np.cos(np.deg2rad(dat.lat))).sum(axis) / np.sum(np.cos(np.deg2rad(dat.lat)))
    return meanw

def bin_data_x_y(dat, xdat = None, xbins = None, ydat = None, ybins = None, nxsplit=None, nysplit=None, compute_trend_aft_xbin=False):
    """ Binning a two dimensional data array.  
        dat = the data to be binned
        xdat = the data to be used to bin the x axis
        xbins (optional) bins for the x-axis binning
        ydat = the data to be used to bin the y axis
        ybins (optional) bins for the y-axis binning
        nxsplit (optional) the number of equal area bins to bin in the x direction
        nysplit (optional) the number of equal area bins to bin in the y direction
        compute_trend_aft_xbin = set to True if you want to compute the trend over years of the xbin-ed data 
    """
 
    if xdat is None:
        xdat = dat

    dims = dat.dims
    xdim = dims[1]
    ydim = dims[0]

    if xbins is not None:
        print("you're binning into bins")
        digitized = np.digitize(xdat, xbins)
        dat_xbind = [ dat.where( digitized == i).mean(xdim) for i in range(1,len(xbins)+1) ]
        dat_xbind = xr.concat(dat_xbind, dim='x')
   
    if xdat is not None:
        print("you're binning the x axis based on a variable")

    if (compute_trend_aft_xbin == True):
        dat_xbind = xr.apply_ufunc( linfit.compute_slope, dat_xbind, vectorize=True,
                                    input_core_dims=[['year']])*dat_xbind.year.size

    if ydat is not None:
        print("you're binning the y axis based on a variable")
        # splitting into nysplit equal area bins
        indices = np.argsort(ydat)
        dat_xbind = dat_xbind.transpose(ydim,'x')
        dat_xbind_sort = dat_xbind[indices,:]
        area = np.ones(len(indices))*np.cos(np.deg2rad(dat_xbind_sort.lat))
        cum_area = area.cumsum()/area.sum()
        idx = np.searchsorted(cum_area, np.linspace(0,1,nysplit, endpoint=False)[1:])
        area_chunks = np.split(area,idx)
        dat_xbind_chunks = np.split(dat_xbind_sort, idx)
        dat_xybind=[]
        for i in np.arange(0,nysplit,1):
            dat_xybind.append(meanw(dat_xbind_chunks[i], axis=ydim))
        dat_xybind = xr.concat(dat_xybind, dim='y')

    return dat_xybind
       


def plot_bin_data_x_y(fig, bindat, xbins, ybins, ci, cmin, cmax, x1, x2, y1, y2,titlestr=None,xlabel=None,ylabel=None):
    ax = fig.add_axes([x1, y1, (x2-x1), (y2-y1)])
    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.precip_cmap(nlevs)
    norm = BoundaryNorm(clevs, ncolors=mymap.N, clip=True)

    ax.pcolormesh(xbins, ybins, bindat, cmap=mymap, norm=norm, edgecolor='none')
    if titlestr is not None:
        ax.set_title(titlestr)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    return ax
    


