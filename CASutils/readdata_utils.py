# routines for reading in data in various forms

import xarray as xr
import pandas as pd
import numpy as np
from pandas import Timedelta as timedelta
import sys

def read_cesm_h0(fpath, datestart, dateend, var):
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
    except:
        print("warning, you're reading CESM data but there's no time_bnds")
        print("make sure you're reading in what you're expecting to")

    dat = dat.sel(time=slice(datestart, dateend))


    return dat

def read_sfc_cesm_dailyavg(filepath, datestart, dateend):
    """Read in a time slice of a surface field from datestart to dateend.
    Adapted for CESM's wierd calendar.  Setting the time axis as the average 
    of time_bnds
    Here specifying the hour as 12:00:00 to make sure we don't take the first timestamp
    which is actually dec 31st from the year before the start
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """
   
    try:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = True)
    except:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = False)
        dat = xr.decode_cf(dat, use_cftime = True)
        #dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    try:
        dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
    except: pass

    # setting the time axis as the average of time bounds.
    try:
        try:
            timebndavg = np.array(dat.time_bnds, 
                     dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
        except:
            timebndavg = np.array(dat.time_bounds,
                     dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')

        # manually fix the average of time bounds at 02/29-03/02 for leap years
        dates = pd.DatetimeIndex(timebndavg)
        lyindices = np.argwhere( (dates.month == 2) & (dates.day == 29))
        if (len(lyindices) > 0):
            for i in lyindices:
                timebndavg[i] = str(dates.year[i].item())+"-02-28T12:00:00"


        dat['time'] = timebndavg
        dat = dat.sel(time=slice(datestart+"T12:00:00", dateend+"T12:00:00"))
    except:
        print("warning, you're reading CESM data but there's no time_bnds")
        print("make sure you're reading in what you're expecting to")
        dat = dat.sel(time=slice(datestart,dateend))



    return dat

def read_sfc_cesm_3hourly(filepath, datestart, dateend):
    """Read in a time slice of a 3 hourly surface field from datestart to dateend.
    Adapted for CESM's wierd calendar.  Setting the time axis to the average of time_bnds.
    """
    dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times=True)
    timebndavg = np.array(dat.time_bnds,
                 dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')

    dates = pd.DatetimeIndex(timebndavg)
    lyindices = np.argwhere( (dates.month == 2) & (dates.day == 29))
    if (len(lyindices) > 0):
        for i in lyindices:
            timebndavg[i] = str(dates.year[i].item())+"-02-28T22:30:00"

    dat['time'] = timebndavg
    dat = dat.sel(time=slice(datestart+"T01:30:00", dateend+"T22:30:00"))

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

def read_zonalmean(filepath, datestart, dateend):
    """Read in a time slice from datestart to dateend and calculate the zonal mean.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = path to files e.g., "/path/to/files/*.nc"
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """

    try:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override",
                 decode_times=True, use_cftime=True).\
                 sel(time=slice(datestart, dateend))

        try:
            datzm=dat.mean(dim="lon")
        except:
            # deal with problematic coordinate names
            dat=dat.rename({"longitude":"lon", "latitude":"lat"})
            datzm=dat.mean(dim="lon")

    except:
        print("Something's wierd about the time axis, decoding manually")
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override",
                   decode_times=False)
    
        try:
            datzm=dat.mean(dim="lon")
        except:
            # deal with problematic coordinate names
            dat=dat.rename({"longitude":"lon", "latitude":"lat"})
            datzm=dat.mean(dim="lon")

        datzm=xr.decode_cf(datzm, use_cftime=True)
        datzm=datzm.sel(time=slice(datestart, dateend))
        datetimeindex=datzm.indexes['time'].to_datetimeindex()
        datzm['time'] = datetimeindex

    return datzm

