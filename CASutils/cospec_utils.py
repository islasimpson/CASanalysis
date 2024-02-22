import xarray as xr
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
import sys
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit

#****************************************************************************
#* Calculation of wavenumber-frequency cospectra for two time series        *
#* Inputs must have dimentions time and lon                                 *
#*                                                                          *
#* Isla Simpson (26th Nov 2023)                                             *
#****************************************************************************

def cospeccalc(x1, x2, lat_bnds=None, dosymmetries=False, dosegs=False, segsize=None, noverlap=0, ftaper=0.1, spd=1, deseas=True, nharm=4, detrend=True, dw=None, dk=None):
    """ Calculating eddy cos-spectra
        Inputs: x1 (first data series (nt,nlat,nlon))
                x2 (second data series (nt,nlat,nlon))
                lat_bnds (latitude bounds over which to calculate cospectra)
                --- If (lat_bnds = None) cospectra are calculated over all lats
                dosymmetries (if True cospectra are calculated for 
                              symmetric and antisymmetric components separately)
                dosegs (if True then divide into segments)
                segsize (if not None then divide up the data into segments of length segsize)
                noverlap = the number of days to overlap from one segment to the next
                ftaper = the fraction of the segment that is tapered. Set to None if you don't
                         want to taper (not recommended unless testing with idealized waves)
                spd = number of samples per day
                deaseas = if True then data is deseasonalized using nharm harmonics 
                          of the seasonal cycle.
                nharm = the number of harmonics to use for the seasonal cycle when deseasonalizing
                detrend = if True then data is detrended.
                dw = width of frequency bins for output.  If None then the default is used
                dk = width of wavenumber bins for output.  If None then the default is used

        Output:
                cospec = wavenumber-frequency cospectra retaining the latitude dimension
                cospec_avg = cospec averaged over the latitude dimension with area weighting

        The units of the output are such that if you sum over wavenumber and frequency, 
        you should approximately recover the eddy flux, minus the effects of the tapering
        and the stationary waves.
    """

    #--- Check the time series are the same length
    if (x1.time.size != x2.time.size):
        logging.info(
         f"time series don't have the same length.  ntime1={x1.time.size}, ntime2={x2.time.size}")
        logging.info('exiting')
        sys.exit()
    if (x1.lon.size != x2.lon.size):
        logging.info(
         f"time series don't have the same lons.  nlon1={x1.lon.size}, nlon2={x2.lon.size}")
        logging.info('exiting')
        sys.exit()
    if ( ((dosegs) & (segsize is None) | (~(dosegs) & (segsize is not None) ) )):
        logging.info(
         f"dosegs = {dosegs} but segsize = {segsize}.")
        logging.info('exiting')

    #--- Select latitude bounds for the calculation
    if lat_bnds is not None:
        assert isinstance(lat_bnds, tuple)
        x1 = x1.sel(lat=slice(*lat_bnds))
        x2 = x2.sel(lat=slice(*lat_bnds))

    #--- Deseasonalize if necessary
    if (deseas):
        x1 = deseasonalize(x1,nharm)
        x2 = deseasonalize(x2,nharm)

    if (detrend):
        logging.info(
         f"linearly detrending")
        x1 = linfit.lineardetrend(x1, 'time')
        x2 = linfit.lineardetrend(x2, 'time')

    #--- Divide data into symmetric and asymmetric components if necessary
    if (dosymmetries):
        x1 = decompose2SymAsym(x1)
        x2 = decompose2SymAsym(x2)

    #--- Divide up into segments if necessary
    if (dosegs):
        x1_win = getwindows(x1, segsize, noverlap)
        x2_win = getwindows(x2, segsize, noverlap)
    else:
        x1_win = x1
        x2_win = x2
        segsize=x1.time.size

    #--- Apply the taper if necessary
    if (ftaper is not None):
        x1_wintap = apply_taper(x1_win, ftaper)
        x2_wintap = apply_taper(x2_win, ftaper)
    else:
        x1_wintap = x1_win
        x2_wintap = x2_win


    #--- Get longitude and time axis info
    lonaxis = x1_wintap.dims.index('lon')
    timeaxis = x1_wintap.dims.index('time')

    lon_size = x1_wintap.lon.size
    time_size = x1_wintap.time.size

    logging.info(f'lon axis is {lonaxis}')
    logging.info(f'time axis is {timeaxis}')

    #--- Get the wavenumbers and frequencies
    k = np.fft.fftfreq(x1_wintap.lon.size, 1./x1_wintap.lon.size)
    w = np.fft.fftfreq(x1_wintap.time.size, 1./spd)

    #--- Do the FFT in longitude
    x1lonft = np.fft.fft(x1_wintap, axis=lonaxis) / lon_size
    x2lonft = np.fft.fft(x2_wintap, axis=lonaxis) / lon_size

    #--- Get the sine and cosine parts
    x1cos = np.real(x1lonft)
    x1sin = -1.*np.imag(x1lonft)

    x2cos = np.real(x2lonft)
    x2sin = -1.*np.imag(x2lonft)

    #--- perform FFT of the sine and cosine part in time
    x1cosft = np.fft.fft(x1cos, axis=timeaxis) / time_size
    x1ac = np.real(x1cosft)
    x1bc = -1.*np.imag(x1cosft)

    x1sinft = np.fft.fft(x1sin, axis=timeaxis) / time_size
    x1as = np.real(x1sinft)
    x1bs = -1.*np.imag(x1sinft)

    x2cosft = np.fft.fft(x2cos, axis=timeaxis) / time_size
    x2ac = np.real(x2cosft)
    x2bc = -1.*np.imag(x2cosft)

    x2sinft = np.fft.fft(x2sin, axis=timeaxis) / time_size
    x2as = np.real(x2sinft)
    x2bs = -1.*np.imag(x2sinft)

    #--- Calculate the cospectrum
    west=(((x1ac-x1bs)*(x2ac-x2bs))+((-x1bc-x1as)*(-x2bc-x2as)))*(2.)
    east=(((x1ac+x1bs)*(x2ac+x2bs))+((x1bc-x1as)*(x2bc-x2as)))*(2.)

    #--- Set up the coordinates for the westward and eastward xarray data arrays
    outcoords=[]
    dimnames=x1_wintap.dims
    for icoord in range(0,len(dimnames),1):
        if (dimnames[icoord] == 'lon'):
            outcoords.append( ('k', k))
        elif (dimnames[icoord] == 'time'):
            outcoords.append( ('w', w))
        else:
            outcoords.append( (dimnames[icoord], x1_wintap[dimnames[icoord]].values))

    #--- Set up the xr versions of west and east
    west_xr = xr.DataArray(west, coords=outcoords, name='west')
    east_xr = xr.DataArray(east, coords=outcoords, name='east')

    #---Only take the frequencies greater than zero
    #---Delineating westward and eastward by +ve and -ve k
    west_xr = west_xr.where( west_xr.w >= 0, drop=True)
    east_xr = east_xr.where( east_xr.w >= 0, drop=True)


    #---Currently k is ordered 0,1,2,.....,-2,-1
    #---Currently w is ordered 0,1,2,....nlon/2-1,-nlon/2....segsize/2.-1
    #---Instead ordering k from -nlon /2, nlon/2
    kro = np.arange(-lon_size/2, (lon_size/2)+1, 1, dtype=int) # reordered k's


    #--- set up the dimensions and coordinates of the output
    #    Replace k in the current arrays with kro
    dimnames = west_xr.dims
    dimshape = list(west_xr.shape)
    coords=west_xr.coords
    kindex = west_xr.dims.index('k')
    dimshape[kindex]+=1 # increasing length of k dimension by 1

    outcoords=[]
    for icoord in range(0,len(dimnames),1):
        if (dimnames[icoord] == 'k'):
            outcoords.append( (dimnames[icoord], kro) )
        else:
            outcoords.append( (dimnames[icoord], west_xr[dimnames[icoord]].values))

    #---setting up the data array for the output
    cospec = xr.DataArray(np.zeros(dimshape),
                          coords=outcoords, name='cospec')


    #---transpose cospec, west_xr, east_xr to have the last two dimensions as (...,w,k)
    cospec = cospec.transpose(...,'w','k')
    west_xr = west_xr.transpose(...,'w','k')
    east_xr = east_xr.transpose(...,'w','k')

    #---assign westward and eastward propraging waves to cospec
    #---taking only the +(w,k) parts of the eastward and westward arrays.
    cospec[...,0:int(cospec.k.size/2)] =  np.array(west_xr[...,int(west_xr.k.size/2.):0:-1])
    cospec[...,int((cospec.k.size/2)+1):cospec.k.size] = \
         np.array(east_xr[...,1:int(east_xr.k.size/2)+1])


    #--- Get the default width of frequency and wavenumber bins
    delw = cospec.w.isel(w=2).values - cospec.w.isel(w=1).values # default w width
    delk = cospec.k.isel(k=2).values - cospec.k.isel(k=1).values # default k width

    #--- Bin and sum onto alternate grid if needed
    if dw is not None:
        if (dw < delw):
            logging.info(f"!! Warning, you're trying to upscale the w axis"\
                          +f"Default delw={delw}, requested output dw={dw}")
        cospec = binw(cospec, dw, delw)

    if dk is not None:
        if (dk < delk):
            logging.info(f"!! Warning, you're trying to upscale the k axis"\
                         +f"Default delk={delk}, dequested output dk={dk}")
        cospec = bink(cospec, dk, delk)

    #--- Set up the output units
    unitsout = getunits(x1,x2,dw,dk,delw,delk)
    logging.info(f"Output will be in units {unitsout}")

    #---Assign units
    cospec = cospec.assign_attrs({'units':unitsout})

    print(cospec)

    #---Take the latitudinal average
    cospec_avg = latavgcospec(cospec,dosymmetries,dosegs)
    cospec_avg = cospec_avg.assign_attrs({'units':unitsout})

    #---Determing the output (doesn't make sense) to output cospec if dosymmetries=True
    cospec = cospec.rename('cospec')
    cospec_avg = cospec_avg.rename('cospec_latavg')
    if (dosymmetries == False):
        out = xr.merge([cospec,cospec_avg])
    else:
        out = cospec_avg

    return out

#--------------- Functions called by cospeccalc ------------
def split_hann_taper(time4taper, ftaper):
    """ Implements `split cosing bell` taper of length ntime where only the
        fraction ftaper of poitns are tapered (combined on both ends).
        This returns a function that tapers to zero on the ends.  
        Taken from Brian Medeiros' wavenumber-frequency code (21st Nov 2023)
    """
    ntime = len(time4taper)
    npts = int(np.rint(ftaper*ntime)) # total size of taper
    taper = np.hanning(npts)
    tapervals = np.ones(ntime)
    tapervals[0:npts//2+1] = taper[0:npts//2 + 1]
    tapervals[-npts//2+1:] = taper[npts//2 + 1:]
    tapervals = xr.DataArray(tapervals, coords=[time4taper.values], 
                   dims=time4taper.dims, name=time4taper.coords)
    return tapervals 


def decompose2SymAsym(dat):
    """ Decompose into symmstric and asymmetric parts
        Mimic Brian Medeiros' code which mimics an NCL function

        Input: dat[...,lat,...]
        Return: array of dimensions [...,lat,...] with the 
                symmetric components in the SH latitudes
                and the asymmetric components in the NH lattudes

        Brian noted that this function produces indistinguishable results from NCL
    """

    lat_dim = dat.dims.index('lat')
    dat_sym = 0.5*(dat.values + np.flip(dat.values, axis=lat_dim))
    dat_asy = 0.5*(dat.values - np.flip(dat.values, axis=lat_dim))
    dat_sym = xr.DataArray(dat_sym, dims=dat.dims, coords=dat.coords)
    dat_asy = xr.DataArray(dat_asy, dims=dat.dims, coords=dat.coords)
    out = dat.copy(deep=True)
    out.loc[{'lat':dat['lat'][dat['lat']<0]}] = dat_sym.isel(lat=dat_sym.lat<0)
    out.loc[{'lat':dat['lat'][dat['lat']>0]}] = dat_asy.isel(lat=dat_asy.lat>0) 
    return out

def deseasonalize(dat, nharm):
    logging.info(
     f"deseasonalizing with {nharm} harmonics of the seasonal cycle")
    # remove leap year Feb 29ths if necessary
    dat = dat.sel(time=~((dat.time.dt.month == 2) & (dat.time.dt.day == 29)))
    daystr = xr.DataArray(dat.indexes['time'].strftime('%m-%d'), 
                           coords=dat.time.coords, name="daystr")
    datseas = dat.groupby(daystr).mean('time')
    dimtime = dat.dims.index('time')
    datseascyc = filt.calc_season_nharm(datseas, nharm, dimtime=dimtime)
    dat = dat.groupby(daystr) - datseascyc
    return dat

def getwindows(dat, segsize, noverlap):
    """Sort out the windows of size segsize with noverlap timestamps between them""" 
    datseg = dat.rolling(time=segsize, min_periods=segsize)
    dat_win = datseg.construct('segtime')
    assert segsize-noverlap > 0, \
      f'Error, segsize={segsize}, noverlap={noverlap}, which gives stride> -' 
    dat_win = dat_win.isel(time=slice(segsize-1,None,segsize-noverlap))
    logging.info(f"dat_win shape is {dat_win.shape}")
    dat_win = dat_win.rename({'time':'iseg'})
    dat_win = dat_win.rename({'segtime':'time'}) 
    return dat_win

def apply_taper(dat, ftaper):
    """Applying the hanning taper.  First subtract the mean, 
       then do the taper then add the mean back"""
    time4taper = dat.time 
    taper = split_hann_taper(time4taper,ftaper)

    #----Remove the time mean before tapering
    dat = dat - dat.mean('time')

    #----Apply the taper
    dattap = dat*taper

    #----Add back in the time average
    dattap = dattap + dat.mean('time')
    return dattap

def binw(dat,dw,delw):
    wvals = dat.w
    wout = np.arange(0-dw/2,np.max(wvals) + delw+ dw/2,dw) # output frequency axis
    dat_binned = dat.groupby_bins('w',wout).sum('w')
    waxis = [ v.mid for v in dat_binned['w_bins'].values]
    dat_binned.coords['w_bins'] = waxis
    dat_binned = dat_binned.rename({'w_bins':'w'})
    return dat_binned

def bink(dat,dk,delk):
    kvals = dat.k
    kout = np.arange(np.min(kvals) - dk/2,np.max(kvals) + delk + dk/2,dk) # output wavenumber axis
    dat_binned = dat.groupby_bins('k',kout).sum('k')
    kaxis = [ v.mid for v in dat_binned['k_bins'].values]
    dat_binned.coords['k_bins'] = kaxis
    dat_binned = dat_binned.rename({'k_bins':'k'})
    return dat_binned


def getunits(x1,x2,dw,dk,delw,delk):
    """ Set up the units of the output 
        Search for units in the input if available
        otherwise set to "x1units" and "x2units"
    """
    try:
        units1 = x1.units
        units2 = x2.units
    except:
        units1 = 'x1units'
        units2 = 'x2units'


    if dw is not None:
        if dk is not None:
            unitsout = units1+'*'+units2+' per '+str(dw)+' w by '+str(dk)+' k box'
        else:
            unitsout = units1+'*'+units2+' per '+str(dw)+' w by '+str(delk)+' k box'
    else:
        if dk is not None:
            unitsout = units1+'*'+units2+' per '+str(delw)+' w by '+str(dk)+' k box'
        else:
            unitsout = units1+'*'+units2+' per '+str(delw)+' w by '+str(delk)+' k box'

    return unitsout

def latavgcospec(cospec,dosymmetries,dosegs):
    """ Take the area weighted latitudinal average """
    if (dosymmetries):
        cospec_sym = cospec.isel(lat=cospec.lat < 0)
        cospec_asy = cospec.isel(lat=cospec.lat > 0)
        cospec_symw = cospec_sym.weighted(np.cos(np.deg2rad(cospec_sym.lat)))
        cospec_asyw = cospec_asy.weighted(np.cos(np.deg2rad(cospec_asy.lat)))
        if (dosegs):
            cospec_sym = cospec_symw.mean(['iseg','lat'])
            cospec_asy = cospec_asyw.mean(['iseg','lat'])
        else:
            cospec_sym = cospec_symw.mean(['lat'])
            cospec_asy = cospec_asyw.mean(['lat'])
        cospec_avg = xr.concat([cospec_sym, cospec_asy], "component")
        cospec_avg = cospec_avg.assign_coords({"component":["symmetric","antisymmetric"]})
    else:
        cospecw = cospec.weighted(np.cos(np.deg2rad(cospec.lat)))
        if (dosegs):
            cospec_avg = cospecw.mean(['iseg','lat'])
        else:
            cospec_avg = cospecw.mean(['lat'])

    return cospec_avg
