import xarray as xr
import numpy as np
import pandas as pd
from CASutils import filter_utils as filt
from math import nan
from scipy.ndimage import label
from datetime import datetime as dt

def calcshfw(dat, threshold=10, timeaxis='time'):
    """ Calculate the day of year and date of the SH final warming following the method of 
    Black and McDaniel (2007): Interannual Variability in the Southern Hemisphere Circulation
    Organized by Stratospheric Final warming events.

    This method uses 5 day running means at 50hPa, 60S and finds the date at which these
    winds fall below 10m/s and don't return above 10m/s until the following autumn.

    dat = 50hPa zonal wind at 60S (without the 5 day running mean applied)

    """
    ## Smooth the data
    datsm = filt.runningmean(dat, 5, timeaxis=timeaxis, dropna=False)
    testdat = datsm.copy(deep=True)
    testdat = xr.where( testdat < threshold, 1, 0)
    testlabel, testcount = label(testdat) 
    testlabel = xr.DataArray(testlabel, coords=[datsm.time.values], dims=[timeaxis], name='below10')

    dates=[]
    dayofyear=[]
    for icount in np.arange(0,testcount,1):
        thisone = datsm.where( testlabel == (icount+1), drop=True)
        length = len(thisone)
        if (length > 50): # making sure there's more than 50 days before it returns
            dates.append(thisone.time.isel(time=0))
            dayofyear.append(thisone.time.isel(time=0).dt.dayofyear.values)

    dates = xr.concat(dates, dim=timeaxis)
    dayofyear = np.array(dayofyear)
    dayofyear = xr.DataArray(np.array(dayofyear), coords=[dates[timeaxis]], dims=[timeaxis], 
           name='dayofyear')

    # Make sure there's at least 20 days between the start of the data and the first date
    test = (pd.to_datetime(dates[0].values) - pd.to_datetime(datsm.time.isel(time=0).values)).days
    if (test < 20):
        dayofyear = dayofyear[1:len(dayofyear)]
        dates = dates[1:len(dates)]

    return dayofyear, dates

def calcnhfw(dat, timeaxis='time'):
    """Calculate the day of year and the date of the SH final warming following the method of 
       Butler and Gerber(2018) as described in Butler et al (2019) 10.1029/2019GL083346
   
       Uses 60N and 10hPa zonal mean zonal wind.
       Last date during the 1st July to 30th June period of each year on which the
       zonal winds reverse and do not return to westerly for more than 10 consecutive days.
    """
    testdat = dat.copy(deep=True)
    testdat = xr.where( testdat < 0, 1, 0 )
    testlabel, testcount = label(testdat)
    testlabel = xr.DataArray(testlabel, coords=[dat.time.values], dims=[timeaxis], name='below0')

    dates=[]
    for icount in np.arange(0,testcount,1):
        thisone = dat.where( testlabel == (icount+1), drop=True)
        length = len(thisone)
        if (length > 10): 
            dates.append(thisone.time.isel(time=0))
    dates = xr.concat(dates, dim=timeaxis)

    # Check for the latest one
    ybeg = dates[0].dt.year
    yend = dates[len(dates)-1].dt.year
    dates2=[]
    for iyear in np.arange(ybeg,yend+1,1):
        datesuse = dates.where( dates.dt.year == iyear, drop=True)
        if (len(datesuse) > 1):
            june30 = dt.strptime(str(iyear)+'/06/30',"%Y/%m/%d")
            difdays = (june30 - pd.to_datetime(datesuse.values)).days
            minind = np.argmin(difdays)
            dates2.append(datesuse[minind])
        else:
            dates2.append(datesuse)

    dates2 = xr.concat(dates2, dim=timeaxis)
    dayofyear = dates2.dt.dayofyear.values
    dayofyear = xr.DataArray(np.array(dayofyear), coords=[dates2[timeaxis]], dims=[timeaxis],
                   name='dayofyear')

    # Make sure there's at least 20 days between the start of the data and the first date
    test = (pd.to_datetime(dates2[0].values) - pd.to_datetime(dat.time.isel(time=0).values)).days
    if (test < 20):
        dayofyear = dayofyear[1:len(dayofyear)]
        dates2 = dates2[1:len(dates)]

    return dayofyear, dates2

