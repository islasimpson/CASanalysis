import xarray as xr
import numpy as np
from math import nan

def barotropicvortgradient(u):
    """Compute the barotropic vorticity gradient (beta - uyy)"""
    omega = 2*np.pi / (24.*60.*60.)
    a = 6.371e6
    beta = 2*omega*np.cos(np.deg2rad(u.lat))

    latdim = u.dims.index('lat')

    ucosphi = u*np.cos(np.deg2rad(u.lat))
    temp = np.gradient(ucosphi, np.deg2rad(u.lat), axis=latdim)
    temp = xr.DataArray(temp, coords=u.coords, dims=u.dims)
    uy = (1./(a*np.cos(np.deg2rad(u.lat))))*temp
    uy = uy.transpose(*u.dims)
    uyy = np.gradient(uy, np.deg2rad(u.lat), axis=latdim)
    uyy = xr.DataArray(uyy, coords = u.coords, dims=u.dims)

    output = beta - uyy
    return output


