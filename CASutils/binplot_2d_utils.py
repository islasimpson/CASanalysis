import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan

from xhistogram.xarray import histogram
from CASutils import colormap_utils as mycolors
from CASutils import linfit_utils as linfit
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Rectangle

from CASutils import colormap_utils as mycolors
import sys


def meanw(dat,axis):
    meanw = (dat*np.cos(np.deg2rad(dat.lat))).sum(axis) / np.sum(np.cos(np.deg2rad(dat.lat)))
    return meanw

def bin_2d_area_weighted(dat, xdat, ydat, xbins, ybins, dim="z", weight=True): 
    """
    Bin `dat` in two dimensions according to `xdat` and `ydat`,
    using cosine-latitude area weighting.

    Inputs:
    - dat = xr.DataArray must contain dimension 'dim'
    - xdat, ydat = xr.DataArrays defining x and y values for binning
    - xbins, ybins = arrays containing bin edges
    - dim = str, the spatial dimension over which to bin
 
    Returns:
    """

    # Cosine-latitude area weights
    if (weight):
        weights = np.cos(np.deg2rad(dat.lat))
        weights = weights.broadcast_like(dat)
    else:
        weights = xr.ones_like(dat)
        #weights = np.zeros([dat.size]) + 1
    #weights = xr.DataArray(weights, dims=[dim], name='wgt')

    valid = (
       dat.notnull()
       & xdat.notnull()
       & ydat.notnull()
       & weights.notnull()
    )

    dat = dat.where(valid)
    xdat = xdat.where(valid)
    ydat = ydat.where(valid)
    weights = weights.where(valid)

    # Sum of weighted dat in each bin
    weighted_sum = histogram(
        xdat, 
        ydat, 
        bins=[xbins, ybins],
        weights = dat*weights,
        dim=[dim]) 

    # Sum of area weights in each bin
    area_sum = histogram(
        xdat,
        ydat,
        bins=[xbins, ybins],
        weights=weights,
        dim=[dim],
    )

    # Number of valid grid points in each bin
    count = histogram(
        xdat,
        ydat,
        bins = [xbins,ybins],
        dim=[dim],
    )

    # Area mean
    area_mean = weighted_sum / area_sum

    # Fraction of the total valid area in each bin
    area_fraction = area_sum / weights.sum(dim)

    # Compute the standard error 
    if (weight == False):
        dat1sum = weighted_sum
    else:
        dat1sum = histogram(
            xdat,
            ydat, 
            bins=[xbins,ybins],
            weights = dat,
            dim = [dim])
    dat2sum = histogram(
            xdat,
            ydat,
            bins=[xbins,ybins],
            weights = dat**2,
            dim = [dim])
    mean_unweighted = dat1sum / count
    std = np.sqrt( (dat2sum / count) - mean_unweighted**2)
    stderr = std / np.sqrt(count)


    datout = xr.Dataset(
        {
           "dat_mean": area_mean,
           "area": area_sum,
           "area_fraction": area_fraction,
           "number": count,
           "standard_deviation": std,
           "standard_error": stderr
        }
    )


    return datout


def plot_2d_binned(fig, dat, xbins, ybins, ci, cmin, cmax, titlestr, x1, x2, y1, y2,
     xlim=None, ylim=None, xticks=None, xticklabels=None, yticks=None, yticklabels=None, 
     xtitle=None, ytitle=None, cmap='blue2red'):
    """ Plot a 2D binned plot, already computed with bin_2d_area_weighted """

    ax = fig.add_axes([x1, y1, (x2-x1), (y2-y1)])
    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)
    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs) 

    norm = BoundaryNorm(clevs, ncolors = mymap.N, clip=True)

    ax.pcolormesh( xbins, ybins, dat, cmap=mymap, norm=norm, edgecolor='none')
  
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xticks is not None:
        ax.set_xticks(xticks)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels, fontsize=14)
    if yticks is not None:
        ax.set_yticks(yticks)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels, fontsize=14)
    if xtitle is not None:
        ax.set_xlabel(xtitle, fontsize=14)
    if ytitle is not None:
        ax.set_ylabel(ytitle, fontsize=14)


    return ax



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
    


