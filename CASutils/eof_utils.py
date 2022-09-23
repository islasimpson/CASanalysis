import xarray as xr
import numpy as np
from eofs.xarray import Eof

def eofcalc_pcnorm(dat, w='sqrtcoslat', neofs=1, timeaxis='time', lonneg=None, latneg=None):
    """ Perform EOF analysis across time for an array with lat and lon dimensions
    Output is normalized such that the PC time series has unit variance 
    and the EOF pattern is in the units of dat.
    Input: dat = array of the form (time,lat,lon) where ... represent the spatial dimensions
           w = If set to sqrtcoslat then sqrt(cos(lat)) weighting is used.
                     Otherwise, set this to an array of the form (lat,lon) that contains the 
                     weights with same spatial dimensions as dat 
           neofs = the number of EOFs to be calculated.

    Output: 
           pcs = PC time series with unit standard deviation
           eofs = EOF pattern with the same units as dat.  Note that this still includes
                  the sqrt(cos(lat))) weighting (I think - determined this by what is needed to 
                  project dat back onto the eof and get the same pc time series back).
    Isla (islas@ucar.edu) Sept 30th, 2021. 
    """

    if (timeaxis != 'time'):
        dat = dat.rename({timeaxis:'time'})

    #move the time axis to the first 
    if (dat.dims[0] != 'time'):
        dat = dat.transpose("time",...)


    if (w == 'sqrtcoslat'):
        w = np.sqrt(np.cos( (dat.lat/180.)*np.pi))
        w = w.expand_dims(dim={'lon':dat.lon.size})
        w = w.transpose()
        w['lon'] = dat.lon

    # Do the EOF calculation
    solver = Eof(dat, weights=w, center=True)
    pcs = solver.pcs(npcs=neofs, pcscaling=1)
    eofs = solver.eofs(neofs=neofs, eofscaling=2)

    # Flipping the EOF based on a reference point if provided.
    # Note this would need to be altered if different reference points were needed for different EOFS
    # Right now all EOFs are flipped.
    if (lonneg and latneg):
        if (eofs.sel(lon=lonneg, lat=latneg, method='nearest') > 0):
            eofs = -1.*eofs
            pcs = -1.*pcs

    return pcs, eofs

def eofcalc_eofnorm(dat, w='sqrtcoslat', neofs=1, timeaxis='time', lonneg=None, latneg=None):
    """ Perform EOF analysis across time for an array with lat and lon dimensions
    Output is normalized such that the PC time series has unit variance 
    and the EOF pattern is in the units of dat.
    Input: dat = array of the form (time,lat,lon) where ... represent the spatial dimensions
           w = If set to sqrtcoslat then sqrt(cos(lat)) weighting is used.
                     Otherwise, set this to an array of the form (lat,lon) that contains the 
                     weights with same spatial dimensions as dat 
           neofs = the number of EOFs to be calculated.

    Output: 
           pcs = PC time series with unit standard deviation
           eofs = EOF pattern with the same units as dat.  Note that this still includes
                  the sqrt(cos(lat))) weighting (I think - determined this by what is needed to 
                  project dat back onto the eof and get the same pc time series back).
    Isla (islas@ucar.edu) Sept 30th, 2021. 
    """

    if (timeaxis != 'time'):
        dat = dat.rename({timeaxis:'time'})

    #move the time axis to the first 
    if (dat.dims[0] != 'time'):
        dat = dat.transpose("time",...)


    if (w == 'sqrtcoslat'):
        w = np.sqrt(np.cos( (dat.lat/180.)*np.pi))
        w = w.expand_dims(dim={'lon':dat.lon.size})
        w = w.transpose()
        w['lon'] = dat.lon

    # Do the EOF calculation
    solver = Eof(dat, weights=w, center=True)
    pcs = solver.pcs(npcs=neofs, pcscaling=0)
    eofs = solver.eofs(neofs=neofs, eofscaling=2)

    # Flipping the EOF based on a reference point if provided.
    # Note this would need to be altered if different reference points were needed for different EOFS
    # Right now all EOFs are flipped.
    if (lonneg and latneg):
        if (eofs.sel(lon=lonneg, lat=latneg, method='nearest') > 0):
            eofs = -1.*eofs
            pcs = -1.*pcs

    return pcs, eofs










def proj_onto_eof(dat, eof, w='sqrtcoslat'):
    """ project array dat onto the eof that has been obtained using "eofcalc_pcnorm".
        Input: dat = the array to be projected onto the EOF of the form (time, lat, lon)
               eof = the EOF calculated using "eofcalc_pcnorm) of the form (lat,lon)
               w = If set to sqrtcoslat then sqrt(cos(lat)) weighting is used.
                   Otherwise, set this to an array of the form (lat, lon) that contins 
                   the weights with the same spatial dimensions as dat
    """

    if (w == 'sqrtcoslat'):
        weights = np.sqrt(np.cos( dat.lat/180.)*np.pi)
        weights = weights.expand_dims(dim={'lon':dat.lon.size})
        weights = weights.transpose()
        weights['lon'] = dat.lon

    num = (dat*weights).dot(eof, dims=['lat','lon'])
    denom = (eofs).dot( (eofs), dims=['lat','lon'])
    proj = num/denom
    
    return proj



