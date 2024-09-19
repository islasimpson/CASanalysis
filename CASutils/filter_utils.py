import numpy as np
from scipy.fft import fft, ifft
import xarray as xr
import sys
import logging
logging.basicConfig(level=logging.INFO)

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

def runningmean(dat, nysm, timeaxis='time', dropna=False, center=True):
    """dat = your data with a time axis with name equal to whatever you set "timeaxis" to
       nysm = the number of time values in your running mean
       dropna = False if you don't want to drop the NaN's at the edges
    """

    window_kwargs = {timeaxis:nysm}
    if (dropna):
        datm = dat.rolling(center=center, min_periods=nysm, **window_kwargs).mean(timeaxis).dropna(timeaxis)
    else:
        datm = dat.rolling(center=center, min_periods=nysm, **window_kwargs).mean(timeaxis)
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

def wkfilter_flux(dat1, dat2, ftaper, kmin = None, kmax = None, pmin = None, pmax = None, spd=1):
    """
    dat1 = an xarray DataArray with dimensions time and lon and whatever other dimensions you like
    dat2 = as dat1
    the output will be the filtered flux dat1*dat2
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
    if ((pmax is not None) & (pmin is not None)):
        fmin = 1./pmax ; fmax = 1./pmin
        assert ( (pmin >= 0) & (pmin >= 0)), f"Error, you can't choose negative periods/frequencies"

    #---Find the dimensions other than lon and time 
    dims = np.array(dat1.dims)
    coords = dat1.coords

    #---Check the dataset has the required dimensions and that only positive frequencies are used
    assert (dims == 'time')[ ~((dims == 'time') == False) ][0] == True, f"Error, you need to have a time dimension"
    assert (dims == 'lon')[ ~((dims == 'lon') == False) ][0] == True, f"Error, you need to have a lon dimension"

    #---Taper the time series
    npts = int(np.rint(ftaper*dat1.time.size))
    endtaper = np.hanning(npts)
    taper = np.ones(dat1.time.size)
    taper[0:npts//2+1] = endtaper[0:npts//2+1]
    taper[-npts//2+1:] = endtaper[npts//2+1:]
    taper = xr.DataArray(taper, coords=[dat1.time], dims=['time'], name='taper')
    dat1_tp = taper*dat1
    dat2_tp = taper*dat2

    #---Dot he FFTs and organize the dimensions
    londim = dat1_tp.dims.index('lon')
    timedim = dat1_tp.dims.index('time')

    dat1fft = np.fft.fft(dat1_tp, axis=londim) / dat1_tp.lon.size
    dat1fft = np.fft.fft(dat1fft, axis=timedim) / dat1_tp.time.size

    dat2fft = np.fft.fft(dat2_tp, axis=londim) / dat2_tp.lon.size
    dat2fft = np.fft.fft(dat2fft, axis=timedim) / dat2_tp.time.size

    dimsfft = []
    coordsfft = []
    for i in np.arange(0,len(dims),1):
        if (dims[i] == 'time'):
            dimsfft.append('frequency')
            frequency = xr.DataArray(np.fft.fftfreq(dat1.time.size, 1./spd), dims=['frequency'])
            coordsfft.append(frequency)
        elif (dims[i] == 'lon'):
            dimsfft.append('wavenumber')
            wavenumber = xr.DataArray(-1.*np.fft.fftfreq(dat1.lon.size, 1./dat1.lon.size), dims=['wavenumber'])
            coordsfft.append(wavenumber)
        else:
            dimsfft.append(dims[i])
            coordsfft.append(coords[dims[i]])

    dat1fft = xr.DataArray(dat1fft, dims=dimsfft, coords=coordsfft)
    dat2fft = xr.DataArray(dat2fft, dims=dimsfft, coords=coordsfft)

    if ( ( kmin is None ) | (kmax is None) | (pmin is None) | (pmax is None) ):
        print('Decomposing into easterly and westerly because not all of kmin, kmax, pmin, pmax were specified')
        dat1_e = dat1fft.where( (( dat1fft.wavenumber >= 1) & (dat1fft.wavenumber <= np.max(dat1fft.wavenumber))), 0)
        dat1_e = dat1_e.where( (dat1_e.frequency >= 0) & (dat1_e.frequency <= np.max(dat1_e.frequency)),0)
        
        dat2_e = dat2fft.where( (( dat2fft.wavenumber >= 1) & (dat2fft.wavenumber <= np.max(dat2fft.wavenumber))), 0)
        dat2_e = dat2_e.where( (dat2_e.frequency >= 0) & (dat2_e.frequency <= np.max(dat2_e.frequency)),0)
        
        dat1_w = dat1fft.where( (( dat1fft.wavenumber >= np.min(dat1fft.wavenumber)) & (dat1fft.wavenumber <= -1)), 0)
        dat1_w = dat1_w.where( (dat1_w.frequency >= 0) & (dat1_w.frequency <= np.max(dat1_w.frequency)),0)
        
        dat2_w = dat2fft.where( (( dat2fft.wavenumber >= np.min(dat2fft.wavenumber)) & (dat2fft.wavenumber <= -1)), 0)
        dat2_w = dat2_w.where( (dat2_w.frequency >= 0) & (dat2_w.frequency <= np.max(dat2_w.frequency)),0)
        
        idat1_e = np.fft.ifft(dat1_e, axis=timedim)*dat1_e.frequency.size
        idat1_e = np.fft.ifft(idat1_e, axis=londim)*dat1_e.wavenumber.size
        idat1_e = 2*np.real(idat1_e)
        
        idat2_e = np.fft.ifft(dat2_e, axis=timedim)*dat2_e.frequency.size
        idat2_e = np.fft.ifft(idat2_e, axis=londim)*dat2_e.wavenumber.size
        idat2_e = 2*np.real(idat2_e)
        
        idat1_w = np.fft.ifft(dat1_w, axis=timedim)*dat1_w.frequency.size
        idat1_w = np.fft.ifft(idat1_w, axis=londim)*dat1_w.wavenumber.size
        idat1_w = 2*np.real(idat1_w)
        
        idat2_w = np.fft.ifft(dat2_w, axis=timedim)*dat2_w.frequency.size
        idat2_w = np.fft.ifft(idat2_w, axis=londim)*dat2_w.wavenumber.size
        idat2_w = 2*np.real(idat2_w)
        
        flux_e = idat1_e*idat2_e
        flux_w = idat1_w*idat2_w

        flux_e = xr.DataArray(flux_e, dims=dims, coords=coords, name='eastward')
        flux_w = xr.DataArray(flux_w, dims=dims, coords=coords, name='westward')

        datout = xr.merge([flux_e, flux_w])
    else:
        dat1_select = dat1fft.where( ((dat1fft.wavenumber >= kmin) & (dat1fft.wavenumber <= kmax)) &
                       ((dat1fft.frequency >= fmin) & (dat1fft.frequency <= fmax)), 0)
        dat2_select = dat2fft.where( ((dat2fft.wavenumber >= kmin) & (dat2fft.wavenumber <= kmax)) & 
                       ((dat2fft.frequency >= fmin) & (dat2fft.frequency <= fmax)), 0)

        idat1 = np.fft.ifft(dat1_select, axis=timedim)*dat1_select.frequency.size
        idat1 = np.fft.ifft(idat1, axis=londim)*dat1_select.wavenumber.size
        idat1 = 2*np.real(idat1)

        idat2 = np.fft.ifft(dat2_select, axis=timedim)*dat2_select.frequency.size
        idat2 = np.fft.ifft(idat2, axis=londim)*dat2_select.wavenumber.size
        idat2 = 2*np.real(idat2)

        flux = idat1*idat2

        datout = xr.DataArray(flux, dims=dims, coords=coords)
    

    return datout























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

def flanczos(dat, delt, nw, pcoff, frequency=None):
    """
    Filter data using a lanczos filter
    Inputs:
        dat=the data you want to filter (must have axis "time")
        delt=the timestep in the same units as period
        nw=the number of weights to be used (should be an odd number)
        pcoff=the cut of period 
        frequency (optional)=output frequency for response function

    Output: datout which contains
        datlow=low pass filtered data
        dathigh = high pass filtered data
        response = response function

    Warning: hasn't been extensively tested for delt!=1, but I think it works
    """
    
    fcoff=(1/pcoff)*delt # cut of frequency in units of delt

    # check that the cutoff frequency is less than the nyquist frequency
    if ( fcoff > 0.5 ):
        logging.info("Cut of frequency {fcoff/delt} is greater than nyquist")
        sys.exit()

    # compute the lanczos weights
    n=int((nw-1)/2)
    k=np.arange(1,n)
    w=np.zeros([nw])
    w[n]=2*fcoff
    whalf=(np.sin(2. * np.pi * fcoff * k) / (np.pi * k))*(np.sin(np.pi * k / n) * n / (np.pi * k))
    w[n-1:0:-1]=whalf
    w[n+1:-1]=whalf
    w=xr.DataArray(w,dims=['window'])

    # compute the response function
    if frequency is None:
        frequency=np.arange(0,0.5,0.01)
    fdelt=frequency*delt # frequency in units of delt
    r=[ (w[n] + 2*np.sum(whalf[0:n]*np.cos(k[0:n]*i*2*np.pi))) for i in fdelt ]
    r = xr.DataArray(r, dims=['frequency'], coords=[frequency], name='response')

    # low pass filtered data
    datlow = dat.rolling(time=len(w), center=True).construct('window').dot(w)

    # high pass filtered data
    dathigh = dat - datlow

    datlow = datlow.rename('lowf')
    dathigh = dathigh.rename('highf')

    datout = xr.merge([datlow, dathigh, r])
    return datout,w


 








