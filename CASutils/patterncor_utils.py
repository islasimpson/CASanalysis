import xarray as xr
import numpy as np
from math import nan

def weightedspatialcor(dat1, dat2, w, londim='lon', latdim='lat'):
    """ 
    Calculate the weighted spatial pattern correlation between datasets dat1 and dat2
    Inputs:
        dat1 = the first data array with last two dimensions (latdim, londim)
        dat2 = the second data array with last two dimensions (latdim, londim)
        w = the weights with dimensions (latdim, londim)
        londim = option to specify alternate longitude dimension
        latdim = option to specify alternate latitude dimension

    Outputs:
        r = the weighted spatial pattern correlation between the two datasets.

    Note that locations that have NaNs are given zero weight.
    """

    # give zero weight to the locations that have NaNs
    w = w.where( ~(np.isnan(dat1) | np.isnan(dat2) | np.isnan(w)), 0)
    # set data to zero where you have zero weight to avoid having NaNs
    dat1 = dat1.where( w != 0, 0)
    dat2 = dat2.where( w != 0, 0)

    # weight the arrays
    dat1w = dat1.weighted(w)
    dat2w = dat2.weighted(w)
    
    # calculate the weighed mean
    dat1m = dat1w.mean(('lon','lat'))
    dat2m = dat2w.mean(('lon','lat'))
    
    # calculate the anomalies from the weighted means
    dat1anoms = dat1 - dat1m
    dat2anoms = dat2 - dat2m
    
    # weight the anomalies
    dat1anomsw = dat1anoms.weighted(w)
    dat2anomsw = dat2anoms.weighted(w)
    
    # calculate the variance of the weighted anomalies
    dat1anomsw_var = dat1anomsw.var(('lon','lat'))
    dat2anomsw_var = dat2anomsw.var(('lon','lat'))
    
    # calculate the weighted covariance of the anomalies
    cov = (dat1anoms*dat2anoms*w).sum(('lon','lat')) / w.sum(('lon','lat'))
    r = cov / np.sqrt(dat1anomsw_var*dat2anomsw_var)
    
    return r
