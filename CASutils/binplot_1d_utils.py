import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan

from CASutils import colormap_utils as mycolors
from CASutils import linfit_utils as linfit
from matplotlib.colors import BoundaryNorm
import sys
from matplotlib import colors as c

def meanw(dat):
    nanmask = xr.where( ~np.isnan(dat), 1, nan)
    meanw  = ((dat*np.cos(np.deg2rad(dat.lat))).sum()) / (nanmask*np.cos(np.deg2rad(dat.lat))).sum()
    return meanw

def meanw_specifydim(dat,axis):
    nanmask = xr.where( ~np.isnan(dat), 1, nan)
    meanw = (dat*np.cos(np.deg2rad(dat.lat))).sum(axis) / (nanmask*np.cos(np.deg2rad(dat.lat))).sum(axis)
    return meanw

def sort_by_dat_over_space(dat,sortdat,nbins):
    """Sort dataset dat over the spatial ('z') dimension (assuming 1D.
       Sorting based on the values in sortdat which should be a 1D array with size z
       Splitting into nbins over space
    """
    indices = np.argsort(sortdat)
    dat_sort = dat[np.array(indices)]
    area = np.ones(len(indices))*np.cos(np.deg2rad(dat_sort.lat))
    cum_area = area.cumsum() / area.sum()
    idx = np.searchsorted(cum_area, np.linspace(0,1,nbins, endpoint=False)[1:])

    area_chunks = np.split(area, idx)
    dat_chunks = np.split(dat_sort, idx)
    
    allchunks=[]
    for i in np.arange(0,len(dat_chunks),1):
        allchunks.append(meanw_specifydim(dat_chunks[i],'z'))
    allchunks = xr.concat(allchunks, dim='y')
    allchunks['y'] = np.linspace(0,100,nbins) - ((100/nbins)/2.)
    return allchunks


def bin_data_x_y_1d(dat, xdat, ydat, xbins = None, ybins = None, nxsplit=None, nysplit=None):
    """ Binning a 1D array
        dat = the data to be binned
        xdat = the data to be used to bin the x axis
        xbins (optional) bins for the x-axis binning
        ydat = the data to be used to bin the y axis
        ybins (optional) bins for the y-axis binning
        nxsplit (optional) the number of equal area bins to bin in the x direction
        nysplit (optional) the number of equal area bins to bin in the y direction
    """

    #---order the data according to y dat
    indices = np.argsort(ydat)
    datsort = dat[np.array(indices)]
    xdatsort = xdat[np.array(indices)]
    ydatsort = ydat[np.array(indices)]

    if ybins is None: # binning into equal area bins in the y direction
        #---divide the y axis up into equal area bins
        area = np.ones(len(indices))*np.cos(np.deg2rad(datsort.lat))
        cum_area = area.cumsum() / area.sum()
        idx = np.searchsorted(cum_area, np.linspace(0,1,nysplit, endpoint=False)[1:])

        #---divide into equal area chunks
        area_chunks = np.split(area, idx)
        datsort_chunks = np.split(datsort, idx)
        ydatsort_chunks = np.split(ydatsort, idx)
        xdatsort_chunks = np.split(xdatsort, idx)

    
    if xbins is None: # binning into equal area bins in the x direction
        dat_bind=[]
        std_bind=[]
        number=[]
        for i in np.arange(0,len(datsort_chunks),1):
            indices = np.argsort(xdatsort_chunks[i])
            xdatsort2 = xdatsort_chunks[i][np.array(indices)]
            ydatsort2 = ydatsort_chunks[i][np.array(indices)]
            datsort2 = datsort_chunks[i][np.array(indices)]
            area = np.ones(len(indices))*np.cos(np.deg2rad(datsort2.lat))
            cum_area = area.cumsum() / area.sum()
            idx = np.searchsorted(cum_area, np.linspace(0,1,nxsplit, endpoint=False)[1:])

            #---divide into equal area chunks
            area_chunks2 = np.split(area, idx)
            datsort_chunks2 = np.split(datsort2, idx)
            ydatsort_chunks2 = np.split(ydatsort2, idx)
            xdatsort_chunks2 = np.split(xdatsort2, idx)
           
            dat_xbind=[]
            std_xbind=[] 
            number_xbind=[]
            for j in np.arange(0,nxsplit,1):
                datuse = datsort_chunks2[j]
                dat_xbind.append(meanw(datuse))
                #! note - not area weighted
                std_xbind.append(np.std(datuse))
                #---fractional area of the planet in each bin
                number_xbind.append(area_chunks2[j].sum() / np.sum((dat*0 + 1)*np.cos(np.deg2rad(dat.lat))))
            dat_xbind = xr.concat(dat_xbind, dim='x')
            std_xbind = xr.concat(std_xbind, dim='x')
            number_xbind = xr.concat(number_xbind, dim='x')
            dat_bind.append(dat_xbind)
            std_bind.append(std_xbind)
            number.append(number_xbind)
        dat_bind = xr.concat(dat_bind, dim='y')
        dat_bind = dat_bind.rename('dat_bind')
        std_bind = xr.concat(std_bind, dim='y')
        std_bind = std_bind.rename('std_bind')
        number = xr.concat(number, dim='y')
        number = number.rename('area')


    if xbins is not None:
        dat_bind=[]
        std_bind=[]
        number=[]
        for i in np.arange(0,len(xdatsort_chunks),1):
            digitized = np.digitize(xdatsort_chunks[i],xbins) # bin according to xdat

            if (sum == True):
                #!!! Note, not a weighted sum
                dat_bindt = [ datsort_chunks[i].where( digitized == j ).sum() for j in range(1,len(xbins)+1) ]
                dat_bindt = xr.concat(dat_bindt, dim='x')
                numbert = [ xr.DataArray(xr.where( digitized == j, 1, 0).sum()) for j in range(1,len(xbins)+1) ]
                numbert = xr.concat(numbert, dim='x')
            else:
#                dat_bindt = [ datsort_chunks[i].where( digitized == j ).mean() for j in range(1,len(xbins)+1) ]
                dat_bindt = [ meanw(datsort_chunks[i].where(digitized == j)) for j in range(1,len(xbins)+1) ]
                dat_bindt = xr.concat(dat_bindt, dim='x')
                #!!! Note, number is jus tthe number of grid points, not area weighted
                numbert = [ xr.DataArray(xr.where( digitized == j, 1, 0).sum()) for j in range(1,len(xbins)+1) ]
                numbert = xr.concat(numbert, dim='x')

            std_bindt = [ xr.DataArray(np.std(datsort_chunks[i].where( digitized == j))) for j in range(1,len(xbins)+1) ]
            std_bindt = xr.concat(std_bindt, dim='x')


            dat_bind.append(dat_bindt)
            std_bind.append(std_bindt)
            number.append(numbert)

        dat_bind = xr.concat(dat_bind, dim='y')
        dat_bind = dat_bind.rename('dat_bind')
        dat_bind['x'] = xbins
        std_bind = xr.concat(std_bind, dim='y')
        std_bind = std_bind.rename('std_bind')
        std_bind['x'] = xbins
        number = xr.concat(number, dim='y')
        number = number.rename('number')
        number['x'] = xbins

    if (ybins is None):
        dat_bind['y'] = np.arange(0,nysplit,1)
        std_bind['y'] = np.arange(0,nysplit,1)
        number['y'] = np.arange(0,nysplit,1)
    else:
        dat_bind['y'] = ybins
        std_bind['y'] = ybins
        number['y'] = ybins


    if (xbins is None):
        dat_bind['x'] = np.arange(0,nxsplit,1)
        std_bind['x'] = np.arange(0,nxsplit,1)
        number['x'] = np.arange(0,nxsplit,1)
    else:
        dat_bind['x'] = xbins
        std_bind['x'] = xbins
        number['x'] = xbins

    datout = xr.merge([dat_bind, std_bind, number])

    return datout


def plotbindata(fig,dat,ci,cmin,cmax,title,xtitle,ytitle,x1,x2,y1,y2,cmap='blue2red', yticks=None, xticks=None,
               nanvals = None, signifdat=None, stipplesize=2, xlim=None):
    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])
    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    if (cmap == 'blue2red'):
        mymap = mycolors.blue2red_cmap(nlevs)
    if (cmap == 'precip'):
        mymap = mycolors.precip_cmap(nlevs)
    graymap = c.ListedColormap(['lightgray'])   
     
    norm = BoundaryNorm(clevs, ncolors=mymap.N, clip=True)

    dims = dat.dims
    ydim = dims[0] ; xdim = dims[1]
    x = dat[xdim] + (dat[xdim][1] - dat[xdim][0])/2.
    y = dat[ydim] + (dat[ydim][1] - dat[ydim][0])/2.

    ax.pcolormesh(x, y, dat, cmap=mymap, norm=norm, edgecolor='none')
   
    if (signifdat is not None):
        dims = signifdat.dims
        ydim = dims[0] ; xdim = dims[1]
        ydimvals = signifdat[ydim]
        xdimvals = signifdat[xdim]
        signifdat = signifdat.stack(z=(ydim,xdim))
        signifdat = signifdat.dropna('z')
        xscatter = signifdat[xdim] + (xdimvals[1]-xdimvals[0])/2.
        yscatter = signifdat[ydim] + (ydimvals[1]-ydimvals[0])/2.
        ax.scatter(xscatter,yscatter,s=stipplesize,color='black')
 
    if (xticks is not None):
        ax.set_xticks([xticks])
        ax.set_xticklabels([xticks])
    if (yticks is not None):
        ax.set_yticks([yticks])
        ax.set_yticklabels([yticks])

    if (nanvals is not None):
        ax.pcolormesh(x,y,nanvals,cmap=graymap)

    if (xlim is not None):
        ax.set_xlim(xlim)

    ax.set_title(title, fontsize=16)
    ax.set_xlabel(xtitle, fontsize=14)
    ax.set_ylabel(ytitle, fontsize=14)

    return ax

