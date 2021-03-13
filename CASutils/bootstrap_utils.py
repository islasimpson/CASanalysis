import xarray as xr
import numpy as np
import sys
from math import nan

def bootgen(darray, nsamples=None, seed=None, nboots=1000):
    """Generate nboots bootstrap samples from darray with nsamples within each bootstrap.
    Sampling is done from the left most dimension
    If nsamples = None then nsamples = the length of the left most dimension.

    """

    ### exit if darray is a dataset.
    if (str(type(darray)) == "xarray.core.dataset.Dataset"):
        print("this function doesn't accept datasets, convert to data array")
        sys.exit()

    ###if it's an xarray dataset, set up the dimensions
    try:
        dims = darray.dims
        if nsamples is None:
            nsamples = darray[dims[0]].size
 
        nmemin = darray[dims[0]].size

        dimboot = [nsamples*nboots]
        dimboot2d = [nboots, nsamples]
        bootcoords = [('iboot', np.arange(0,nboots,1)), ('isample', np.arange(0,nsamples,1))]
        for icoord in range(1,len(dims)):
            dimboot.append(darray[dims[icoord]].size)
            dimboot2d.append(darray[dims[icoord]].size)
            bootcoords.append( (dims[icoord], darray[dims[icoord]] ))

        print("you are using an xarray dataarray")
   
    except:
        if nsamples is None:
            nsamples = darray.shape[0]
        nmemin = darray.shape[0]

        dimboot = [nsamples*nboots]
        dimboot2d = [nboots, nsamples]
        for icoord in range(1,len(darray.shape)):
            dimboot.append(darray.shape[icoord])
            dimboot2d.append(darray.shape[icoord])

        print("you are using a numpy array")

    ### generate random number for bootstrapping
    if (seed):
        np.random.seed(seed)

    ### do the resampling
    ranu = np.random.uniform(0,nmemin,nboots*nsamples)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.array(darray[ranu])
    bootdat = bootdat.reshape(dimboot2d)

    try:
        bootdat = xr.DataArray(bootdat, coords=bootcoords)
    except:
        pass

    return bootdat



