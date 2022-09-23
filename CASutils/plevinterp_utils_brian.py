import numpy as np
import xarray as xr

from numba import njit, prange

@njit(parallel=True)  # NOTE: if this is breaking, get rid of the decorator, and change 'prange' to range
def vert_remap(x_mdl, p_mdl, plev):
    """Apply simple 1-d interpolation to a field, x
       given the pressure p and the new pressures plev.
       x_mdl, p_mdl are numpy arrays of shape (nlevel, spacetime).

       NOTE: we try to use Numba here to parallelize the loop.
       If Numba is not available, I tried to 'spoof' the decoratory,
       but I have not tested much. If it breaks, remove the njit decorator,
       and change `prange` to `range' in the `for` statement.
    """
    out_shape = (plev.shape[0], x_mdl.shape[1])
    output = np.full(out_shape, np.nan)
    for i in prange(out_shape[1]):
        output[:,i] = np.interp(plev, p_mdl[:,i], x_mdl[:,i])
    return output


def pres_from_hybrid(psfc, hya, hyb, p0=100000.):
    """Derive the pressure field on hybrid-sigma levels."""
    # p = a(k)*p0 + b(k)*ps.
    return hya*p0 + hyb*psfc

def sigma_from_hybrid(psfc, hya, hyb, p0=100000.):
    """Derive the sigma levels from hybrid-sigma levels."""
    return (hya*p0 + hyb*psfc)/psfc


def lev_to_plev(data, ps, hyam, hybm, P0=100000., new_levels=None):
    """Interpolate data from hybrid-sigma levels to isobaric levels.
    
    data : DataArray with a 'lev' coordinate
    ps   : DataArray of surface pressure (Pa), same time/space shape as data
    hyam, hybm : hybrid coefficients, size of len(lev)
    P0 : reference pressure
    new_levels : the output pressure levels (Pa)
    parallel : if True, use the Numba version to parallelize interpolation step.
    """
    pressure = pres_from_hybrid(ps, hyam, hybm, P0)  # Pa
    if new_levels is None:
        pnew = 100.0 * np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1])  # mandatory levels, converted to Pa
    else:
        pnew = new_levels
    # reshape data and pressure assuming "lev" is the name of the coordinate
    zdims = [i for i in data.dims if i != 'lev']
    dstack = data.stack(z=zdims)
    pstack = pressure.stack(z=zdims)
    output = vert_remap(dstack.values, pstack.values, pnew)
    output = xr.DataArray(output, dims=("plev", "z"), coords={"plev":pnew, "z":pstack['z']})
    output = output.unstack()
    return output

def lev_to_sigma(data, ps, hyam, hybm, P0=100000., new_levels=None):
    """Interpolate data from hybrid-sigma levels to sigma levels.
    data: DataArray with a 'lev' coordinate
    ps : DataArray of surface pressure (pa), same time/space shape as data
    hyam, hybm: hybrid coefficients, size of len(lev)
    P0 : reference pressure
    new_levels " the output pressure levels (Pa)
    """

    sigma = sigma_from_hybrid(ps, hyam, hybm, P0)
    if new_levels is None:
        sigmanew = np.array([0.01,0.025,0.0375, 0.0525, 0.07, 0.085, 0.1, 0.12, 0.14, 0.16, 0.19, \
                             0.2275, 0.2675, 0.315, 0.37, 0.435, 0.51, 0.6, 0.6975, 0.7875, 0.865, \
                             0.9275, 0.97, 0.9925 ])
    else:
        sigmanew = new_levels
    # reshape data and pressure assuming "lev" is the name of the coordinate
    zdims = [i for i in data.dims if i != 'lev']
    dstack = data.stack(z=zdims)
    sigstack = sigma.stack(z=zdims)
    output = vert_remap(dstack.values, sigstack.values, sigmanew)
    output = xr.DataArray(output, dims=("sigma","z"), coords={"sigma":sigmanew, "z":sigstack['z']})
    output = output.unstack()
    return output




