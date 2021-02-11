# routines for reading in data in various forms

import xarray as xr
import pandas as pd
import numpy as np
from pandas import Timedelta as timedelta
import sys

def read_cesm_h0(filepath, datestart, dateend, var):
    """Read in a variable from CESM history files.  Adapted for 
    CESM's wierd calendar.  Setting the time axis as the average of time_bnds
    """

    dat = xr.open_mfdataset(fpath, coords="minimal", decode_times=True)
    dat = dat[var]
    try: 
        timebndavg = np.array(dat.time_bnds, 
                     dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
    except:
        timebndavg = np.array(dat.time_bounds,
                     dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
   
    dat['time'] = timebndavg
    dat = dat.sel(time=slice(datestart, dateend))

    return dat



def read_sfc_cesm(filepath, datestart, dateend):
    """Read in a time slice of a surface field from datestart to dateend.
    Adapted for CESM's wierd calendar.  Setting the time axis as the average 
    of time_bnds
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """
    
    dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = True)
    try:
        dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
    except: pass

    # setting the time axis as the verage of time bounds.
    try:
        try:
            timebndavg = np.array(dat.time_bnds, 
                     dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
        except:
            timebndavg = np.array(dat.time_bounds,
                     dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')

        dat['time'] = timebndavg
        dat = dat.sel(time=slice(datestart, dateend))
    except:
        print("warning, you're reading CESM data but there's no time_bnds")
        print("make sure you're reading in what you're expecting to")

    return dat

def read_sfc(filepath, datestart, dateend):
    """Read in a time slice of a surface field from datestart to dateend.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """

    #First try opening and doing the select assuming everything is working ok with the time axis
    try:
        dat = \
        xr.open_mfdataset\
        (filepath, coords="minimal", join="override", decode_times=True, use_cftime=True).\
        sel(time=slice(datestart, dateend))
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
            print("changing longitude --> lon, latitude --> lat")
        except: pass

    except:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = False)
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
        except: pass

        dat = xr.decode_cf(dat, use_cftime = True)
        dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    return dat
