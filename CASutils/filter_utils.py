import numpy as np
from scipy.fft import fft, ifft
import xarray as xr
import sys

def calc_season_nharm(darray, nharms, dimtime=0):
    """ calculate the seasonal cycle defined as the first n-harmonics of the annual 
        time series.  Assumes the first dimension is time unless specified

    Input: darray = a data array 
    output: seascycle = the seasonal cycle
    !!!! Not totally confident this works for arrays with >2 dimensions at this point!!!

    """
    # get the dimensions of the input array
    dims = darray.dims


    # convert to a numpy array
    darray_np = np.array(darray)

    # reorder the axes if dimtime != 0
    if (dimtime != 0):
        darray_np = np.moveaxis(darray_np, dimtime, 0)

    ntime = darray_np.shape[0]
    shapein = darray_np.shape
    # collapse dimensions 1 to n
    if (len(dims) > 1):
        darray_np = darray_np.reshape( darray_np.shape[0], np.prod(darray_np.shape[1::]))

    # check for NaNs
    if (len(dims) > 1):
        numnans = [np.count_nonzero(np.isnan(darray_np[:,i])) for i in range(0,darray_np.shape[1],1)]
        nanels = np.where(np.array(numnans) > 0)

        print("numnans=",numnans)

        # if nan's exist, del with them by linear interpolation
        if np.any(nanels):
            print("you got nans")
            i = np.arange(ntime)
            for j in np.array(nanels[0]):
                 if (np.isfinite(darray_np[:,j]).any()):
                    mask = np.isfinite(darray_np[:,j])
                    darray_np[:,j] = np.interp(i,i[mask],darray_np[mask,j])

    else: # doing this differently for a 1D array
        i = np.arange(ntime) 
        nanels = np.count_nonzero(np.isnan(darray_np))
        if (nanels > 0):
            mask = np.isfinite(darray_np[:])
            darray_np[:] = np.interp(i,i[mask],darray_np[mask])


    tempft = fft(darray_np, axis=0)
    tempft2 = np.zeros_like(tempft)
#    tempft2[0:nharms,:] = tempft[0:nharms,:]
#    tempft2[ntime-nharms+1:ntime+1,:] = tempft[ntime-nharms+1:ntime+1,:]
    tempft2[0:nharms] = tempft[0:nharms]
    tempft2[ntime-nharms+1:ntime+1] = tempft[ntime-nharms+1:ntime+1]



    darray_filtered = np.real(ifft(tempft2, axis=0))

    # reshape array to expand dimensions out again
    darray_filtered = darray_filtered.reshape(shapein)

    
    if (dimtime != 0):
        darray_filtered = np.moveaxis(darray_filtered, 0, dimtime)

    darray_filtered_xr = xr.DataArray(darray_filtered, coords = darray.coords)

    return darray_filtered_xr
