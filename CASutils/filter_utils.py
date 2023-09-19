import numpy as np
from scipy.fft import fft, ifft
import xarray as xr
import sys

def filterk(darray, kmin, kmax, dimlon=0):
    """filter a field based on zonal wavenumber.  Include wavenumbers kmin to kmax"""
    dims = darray.dims
    darray_np = np.array(darray)
    nk = kmax - kmin + 1

    # reorder the axes if dimlon != 0
    if (dimlon != 0):
        darray_np  = np.moveaxis(darray_np, dimlon, 0)

    nlon = darray_np.shape[0]
    shapein = darray_np.shape
    # collapse dimensions 1 to n
    if (len(dims) > 1):
        darray_np = darray_np.reshape(darray_np.shape[0], np.prod(darray_np.shape[1::]))

    # check for NaNs
    if (len(dims) > 1):
        numnans = [np.count_nonzero(np.isnan(darray_np[:,i])) for i in range(0,darray_np.shape[1],1)]
        nanels = np.where(np.array(numnans) > 0)

        # if nan's exist, del with them by linear interpolation
        if np.any(nanels):
            print("you got nans")
            i = np.arange(nlon)
            for j in np.array(nanels[0]):
                 if (np.isfinite(darray_np[:,j]).any()):
                    mask = np.isfinite(darray_np[:,j])
                    darray_np[:,j] = np.interp(i,i[mask],darray_np[mask,j])

    else: # doing this differently for a 1D array
        i = np.arange(nlon)
        nanels = np.count_nonzero(np.isnan(darray_np))
        if (nanels > 0):
            mask = np.isfinite(darray_np[:])
            darray_np[:] = np.interp(i,i[mask],darray_np[mask])

    tempft = fft(darray_np, axis=0)
    tempft2 = np.zeros_like(tempft)
    tempft2[kmin:kmax+1] = tempft[kmin:kmax+1]
    tempft2[nlon - kmin - nk + 1: nlon - kmin + 1] = tempft[ nlon - kmin - nk + 1 : nlon - kmin + 1 ] 

    darray_filtered = np.real(ifft(tempft2, axis=0))

    # reshape array to expand dimensions out again
    darray_filtered = darray_filtered.reshape(shapein)

    
    if (dimlon != 0):
        darray_filtered = np.moveaxis(darray_filtered, 0, dimlon)

    darray_filtered_xr = xr.DataArray(darray_filtered, coords = darray.coords)

    return darray_filtered_xr

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

def runningmean(dat, nysm, timeaxis='time', dropna=False):
    """dat = your data with a time axis with name equal to whatever you set "timeaxis" to
       nysm = the number of time values in your running mean
       dropna = False if you don't want to drop the NaN's at the edges
    """

    window_kwargs = {timeaxis:nysm}
    if (dropna):
        datm = dat.rolling(center=True, min_periods=nysm, **window_kwargs).mean(timeaxis).dropna(timeaxis)
    else:
        datm = dat.rolling(center=True, min_periods=nysm, **window_kwargs).mean(timeaxis)
    return datm

def runningmean_cyclic(dat, nysm, timeaxis='time', dropna=False):
    """dat = your data with a time axis with name equal to whatever you set "timeaxis" to
       nysm = the number of time values in your running mean
       dropna = False if you don't want to drop the NaN's at the edges
    """
    window_kwargs = {timeaxis:nysm}
    npad = np.int((nysm - 1)/2.)
    pad_width_kwargs={timeaxis:npad}
    datm = dat.pad(mode='wrap', **pad_width_kwargs).rolling(center=True, min_periods=1, **window_kwargs).mean(timeaxis).dropna(timeaxis)
    return datm

def wkfilter(dat, ftaper, kmin, kmax, pmin, pmax, spd=1):
    """ 
    dat = an xarray DataArray with dimensions time and lon and watever other dimensions
          you like.
    ftaper = the fraction of the data that's used for tapering.  ftaper/2 will be tapered
             at each end.
    kmin = the minimum wavenumber included in the output
    kmax = the maximum wavenumber included in the output
    pmin = the minimum period (in days) included in the output
    pmax = the maximum period (in days) included in the output
    spd = the number of time stamps per day

    Negative wavenumbers = westward propagating
    Positive wavenumbers = eastward propagating

    Only use positive frequencies as the FFT coefficients are multiplied by 2 before inverting

    """

    #---Convert min/max period to min/max frequency
    fmin = 1./pmax ; fmax = 1./pmin

    #---Find the dimensions other than lon and time 
    dims = np.array(dat.dims)
    coords = dat.coords

    #---Check the dataset has the required dimensions and that only positive frequencies are used
    assert (dims == 'time')[ ~((dims == 'time') == False) ][0] == True, f"Error, you need to have a time dimension"
    assert (dims == 'lon')[ ~((dims == 'lon') == False) ][0] == True, f"Error, you need to have a lon dimension"
    assert ( (pmin >= 0) & (pmin >= 0)), f"Error, you can't choose negative periods/frequencies"



    #---Taper the time series
    npts = int(np.rint(ftaper*dat.time.size))
    endtaper = np.hanning(npts)
    taper = np.ones(dat.time.size)
    taper[0:npts//2+1] = endtaper[0:npts//2+1]
    taper[-npts//2+1:] = endtaper[npts//2+1:]
    taper = xr.DataArray(taper, coords=[dat.time], dims=['time'], name='taper')
    dat_tp = taper*dat

    #---Do the FFTs and organize the dimensions
    londim = dat_tp.dims.index('lon')
    timedim = dat_tp.dims.index('time')

    z = np.fft.fft(dat_tp, axis=londim) / dat_tp.lon.size
    z = np.fft.fft(z, axis=timedim) / dat_tp.time.size

    dimsfft = []
    coordsfft = []
    for i in np.arange(0,len(dims),1):
        if (dims[i] == 'time'):
            dimsfft.append('frequency')
            frequency = xr.DataArray(np.fft.fftfreq(dat.time.size, 1./spd), dims=['frequency'])
            coordsfft.append(frequency)
        elif (dims[i] == 'lon'):
            dimsfft.append('wavenumber')
            wavenumber = xr.DataArray(-1.*np.fft.fftfreq(dat.lon.size, 1./dat.lon.size), dims=['wavenumber'])
            coordsfft.append(wavenumber)
        else:
            dimsfft.append(dims[i])
            coordsfft.append(coords[dims[i]])

    z = xr.DataArray(z, dims=dimsfft, coords = coordsfft)

    #---Select the desired wavenumbers and frequencies
    zselect = z.where( ((z.wavenumber >= kmin) & (z.wavenumber <= kmax)) & 
                       ((z.frequency >= fmin) & (z.frequency <= fmax)), 0)

    #---Invert
    izselect = np.fft.ifft(zselect, axis=timedim)*zselect.frequency.size
    izselect = np.fft.ifft(izselect, axis=londim)*zselect.wavenumber.size

    izselect = 2.*np.real(izselect)

    izselect = xr.DataArray( izselect, dims=dims, coords = coords) 

    return izselect


def wkfilter_flux(x1, x2, ftaper, spd=1):
    """ Function to compute the flux (x1'x2') associated with
        eastward and westward propagating waves separately
        Inputs: x1 (first field)
                x2 (second field)
                ftaper (the fraction of the time series used for tapering f/2 is used at each end)
                spd (the number of time stamps per day)
        Outputs: eastward (the flux due to eastward propagating waves)
                 westward (the flux due to westward propagating waves)
    """
    x1e = wkfilter(x1, ftaper, 1, x1.lon.size, spd, x1.time.size, spd=spd)
    x1w = wkfilter(x1, ftaper, -1*x1.lon.size, -1, spd, x1.time.size, spd=spd)

    x2e = wkfilter(x2, ftaper, 1, x2.lon.size, spd, x2.time.size, spd=spd)
    x2w = wkfilter(x2, ftaper, -1*x2.lon.size, -1, spd, x2.time.size, spd=spd)

    fluxe = x1e*x2e
    fluxw = x1w*x2w

    fluxe = fluxe.rename('eastward')
    fluxw = fluxw.rename('westward')

    flux = xr.merge([fluxe, fluxw])
    return flux
