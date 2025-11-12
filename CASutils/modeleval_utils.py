import xarray as xr
import numpy as np
from math import nan

def nmse(obs, mod):
    """ Calculate the normalized mean squared error metric of Williamson (1995) 
        "Skill Scores from the AMIP Simulations".
        nmse = overbar( (z_m - z_a)**2 ) / overbar( (z_a')**2 )
        where overbar is the weighted spatial average and prime is the deviation 
        from that
    """

    # make sure the model and obs have the same lons and lats
    mod['lon'] = obs.lon ; mod['lat'] = obs.lat
   
    # get the weights and weight by zero if the model or obs is nan
    w = np.cos(np.deg2rad(obs.lat))
    w = w.expand_dims({'lon':obs.lon}, axis=1)
    w = w.where( ~( np.isnan(obs) | np.isnan(mod)), 0)
    obs = obs.where( w != 0, 0)
    mod = mod.where( w != 0, 0)

    # numerator
    num = (mod - obs)**2.
    numw = num.weighted(w)
    numwm = numw.mean(['lon','lat'])
 
    # denomenator
    obsw = obs.weighted(w)
    obswm = obsw.mean(['lon','lat'])
    obsprime = obs - obswm
    obsprime2 = obsprime**2.
    obsprime2w = obsprime2.weighted(w)
    obsprime2wm = obsprime2w.mean(['lon','lat'])

    nmse = numwm / obsprime2wm

    return nmse

def compute_rpc_eade(model, model_em, obs, timeaxis='time'):
    """ Compute the ratio of predictable components"""
    r = xr.corr(model_em, obs, dim=timeaxis)
    sig2_em = model_em.std(timeaxis)**2.
    sig2_ind = model.std(timeaxis)**2.
    sig2_indm = sig2_ind.mean('M')
    rpc = r / np.sqrt( sig2_em / sig2_indm)
    num = r
    denom = np.sqrt( sig2_em / sig2_indm)
    return rpc, num, denom

def compute_rpc_hardiman(model, model_em, obs, timeaxis='time'):
    """ Compute the ratio of predictable components following Hardiman et al
        'Missing eddy feedback may explain weak signal-to-noise ratios in climat predictions"""

    r_mo = xr.corr(model_em, obs, dim=timeaxis)

    sig2signal = model_em.std(timeaxis)**2
    sig2total = model.std(timeaxis)**2
    sig2total = sig2total.mean('M')

    n = model.M.size

    r_mm = np.sqrt( (n*(sig2signal/sig2total) -1) / n-1)

    rpc = r_mo / r_mm

    return rpc, r_mo, r_mm






