## routines for calculating seasonal climatology and seasonal timeseries
import xarray as xr
import numpy as np
from datetime import timedelta, datetime
import pandas as pd
from math import nan
import sys

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

dpseas = {'DJF': 90, 'MAM': 92, 'JJA': 92, 'SON': 91 }

def leap_year(year, calendar='standard'):
    """Determine if year is a leap year
    Args: 
        year (numeric)
    """
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

def get_days_per_mon(time, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    
    Args: time (CFTimeIndex): ie. ds.time.to_index()
          calendar (str): default 'standard'
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if ( (leap_year(year, calendar=calendar)) and (month == 2)):
            month_length[i] += 1
    return month_length


def season_mean(ds, var=None, season = "all", cal = "none"):
    """ calculate climatological mean by season
    Args: ds (xarray.Dataset): dataset
          var (str): variable to use
          season (str): "all", 'DJF', "MAM", "JJA", "SON"
          cal (str): "none"(default) or calendar used for weighting months by number of days
    """

    try:
        ds = ds[var]
    except:
        pass

    ## no weighting of months: 

    if (season == 'JJ'):
        ds = ds.where( (ds['time.month'] == 6) | (ds['time.month'] == 7) )
        smean = ds.mean('time')
    else:        
        if cal == "none":
            if season == "all":
                ## calculate mean for all season
                smean = ds.groupby('time.season').mean('time')
            else :
                ## calculate mean for specified season
                smean = ds.where(ds['time.season'] == season).mean('time')
    
            return smean
        ## weighted months
        else:
            ## create array of month_length (number of days in each month)
            ## assign time coords matching original ds
            month_length = xr.DataArray(get_days_per_mon(ds.time.to_index(), calendar=cal),
                                     coords=[ds.time], name='month_length')
            ## Calculate the weights by grouping by 'time.season'
            weights = month_length.groupby('time.season') / month_length.groupby('time.season').sum()
    
            if season == "all":
                ## calculate weighted mean for all season
                smean = (ds * weights).groupby('time.season').mean('time')
            else :
                ## calculate weighted mean for specified season
                smean = (ds * weights).where(ds['time.season'] == season).mean('time')

    return smean

def season_ts(ds, season, var=None):
    """ calculate timeseries of seasonal averages
    Args: ds (xarray.Dataset): dataset
          var (str): variable to calculate 
          season (str): 'DJF', 'MAM', 'JJA', 'SON'
    """

    if (season == 'JAS'):
        ds_season = ds.where(
        (ds['time.month'] == 7) | (ds['time.month'] == 8) | (ds['time.month'] == 9))

        if (var):
            ds_season = ds_season[var].rolling(min_periods=3, center=True, time=3).mean().dropna("time", how='all')
        else:
            ds_season = ds_season.rolling(min_periods=3, center=True, time=3).mean().dropna("time", how="all")
    elif (season == 'JJAS'):
        ds_season = ds.where(
        (ds['time.month'] == 6) | (ds['time.month'] == 7) | (ds['time.month'] == 8) | (ds['time.month'] ==9))
        if (var):
            ds_season = ds_season[var].rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")
        else:
            ds_season = ds_season.rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")    
    elif (season == 'OND'):
        ds_season = ds.where(
           (ds['time.month'] == 10) | (ds['time.month'] == 11) | (ds['time.month'] == 12))
        if (var):
            ds_season = ds_season[var].rolling(min_periods=3, center=True, time=3).mean().dropna("time", how="all")
        else:
            ds_season = ds_season.rolling(min_periods=3, center=True, time=3).mean().dropna("time", how="all")
    elif (season == 'MJJA'):
        ds_season = ds.where(
        (ds['time.month'] == 5) | (ds['time.month'] == 6) | (ds['time.month'] == 7) | (ds['time.month'] ==8))
        if (var):
            ds_season = ds_season[var].rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")
        else:
            ds_season = ds_season.rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")
    elif (season == 'DJFM'):
        ds_season = ds.where(
        (ds['time.month'] == 12) | (ds['time.month'] == 1) | (ds['time.month'] == 2) | (ds['time.month'] == 3))
        if (var):
            ds_season = ds_season.rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")
        else:
            ds_season = ds_season.rolling(min_periods=4, center=True, time=4).mean().dropna("time", how="all")
    else:
        ## set months outside of season to nan
        ds_season = ds.where(ds['time.season'] == season)

        # calculate 3month rolling mean (only middle months of season will have non-nan values)
        if (var):
            ds_season = ds_season[var].rolling(min_periods=3, center=True, time=3).mean().dropna("time", how='all')
        else:
            ds_season = ds_season.rolling(min_periods=3, center=True, time=3).mean().dropna("time", how="all")

    return ds_season

def group_season_daily(ds,  season, calendar='standard'):
    """ Group daily data in to seasons 
    """

    # change nans to something else so they don't also get dropped in the seasonal grouping
    ds = ds.where(~np.isnan(ds),-999)

    #move the time axis to the first 
    if (ds.dims[0] != 'time'):
        print('moving time axis to the start')
        ds = ds.transpose("time",...)


    ds_season = ds.where(ds['time.season'] == season)
    years = ds['time.year']
    ybeg = np.array(years[0])
    yend = np.array(years[len(years)-1])
    months = ds['time.month']
    mbeg = np.array(months[0])
    mend = np.array(months[len(months)-1])

    if (season == 'DJF'):
        if (mbeg < 12):
            # remove any januarys and februaries in the first year
            ds = ds.where(~((ds['time.year'] == ybeg) & ((ds['time.month'] == 1) | (ds['time.month'] == 2))))

        if (mend > 2):
            # remove december in the last year
            ds = ds.where(~((ds['time.year'] == yend) & (ds['time.month'] == 12)))

    ds_season = ds.where(ds['time.season'] == season).dropna("time", how="all")
    nyears = ds_season.time.size/dpseas[season] 

    print("nyears="+str(nyears))

    # get coords for output
    dims = ds.dims
    dimout=["year", "mon"]
    outcoords = [('year', ybeg + np.arange(nyears)), ('day', np.arange(dpseas[season]))]
    for icoord in range(1,len(dims)):
        dimout.append(dims[icoord])
        outcoords.append( (dims[icoord], np.array(ds[dims[icoord]])))

    # check you have an integer number of years
    if (nyears == int(nyears)):
        # get output dimensions (nyears, ndays, ...)
        outdims = [ int(nyears), dpseas[season] ]
        for i in ds_season.shape[1::]:
            outdims.append(i)
        # reshape array, convert to xarray and assign coords  

        ds_season = ds_season.where(ds_season != -999., nan)


        datout = np.reshape(np.array(ds_season), outdims)
        datout = xr.DataArray(datout, coords = outcoords)
    else:
        print("You don't seem to have the right number of days to have an integer number of "+season+" seasons")
        print("ndays = "+str(ds_season.time.size))
        sys.exit()

    return datout 

def date2fracofyear(time):
    year = time.dt.year.values
    month = time.dt.month.values
    day = time.dt.day.values

    cal_days = dpm['standard']

    fracofyear = [ year[i] + (np.sum(cal_days[0:month[i]]) + day[i])/365. for i in np.arange(0,len(year),1) ]
    return fracofyear



def fracofyear2date(time, caltype='standard'):
    """Convert a time series that is in terms of fractions of a year
    """
    year = time.astype(int)
    ly =np.array( [leap_year(i, calendar=caltype) for i in year])
    day = (time - year)*(365 + ly)
    d = pd.to_timedelta(day, unit="d")
    d1 = pd.to_datetime(year, format="%Y")
    date = d + d1

    #d = timedelta(days= ( time - year )*(365+leap_year(year, calendar=caltype)))
    #d = timedelta(days = np.float((time-year)*365))
    #d1 = datetime(year, 1, 1)
    #date = d + d1
    return date 


def YYYYMM2date(time, caltype='standard'):
    """ Convert a date of the form YYYYMM to a datetime64 object """
    date = pd.to_datetime(time, format='%Y%m')
    return date

def YYYYMMDD2date(date, caltype='standard'):
    time = pd.to_datetime(date, format='%Y%m%d')
    return time

def MMDD2date(date, caltype='standard'):
    time = pd.to_datetime(date, format='%m%d')
    return time

def group_daily2yearly(dat):
    years = dat['time.year']
    ybeg = np.array(years[0])
    yend = np.array(years[len(years)-1])

    datyear=[]
    for iyear in np.arange(ybeg,yend+1,1):
        yearlydat = dat.sel(time=slice(str(iyear)+'-01-01', str(iyear)+'-12-31'))
        yearlydat['time'] = np.arange(0,365,1)
        datyear.append(yearlydat)

    datyear = xr.concat(datyear, dim='year')
    return datyear

def group_monthly2yearly(dat):
    years = dat['time.year']
    ybeg = np.array(years[0])
    yend = np.array(years[len(years)-1])

    datyear=[]
    for iyear in np.arange(ybeg,yend+1,1):
        yearlydat = dat.sel(time=slice(str(iyear)+'-01-01', str(iyear)+'-12-31'))
        yearlydat['time'] = np.arange(0,12,1)
        datyear.append(yearlydat)

    datyear = xr.concat(datyear, dim='year')
    return datyear


def calcannualmean(ds, skipna=False):
    """ Calculate the annual mean weighting the months of the year appropriately if
        given the calendar type
    """

    def dothecalc(var, skipna=False):
        month_length = var.time.dt.days_in_month
        wghts = month_length.groupby('time.year') / month_length.groupby('time.year').sum()
        if (skipna):
            datsum = (var*wghts).groupby('time.year').sum(dim='time', skipna=True)
            cond = var.isnull()
            ones = xr.where(cond, 0, 1)
            onesum = (ones*wghts).groupby('time.year').sum(dim='time')
        else:
            datsum = (var*wghts).groupby('time.year').sum(dim='time', skipna=False)
            cond = var.isnull()
            ones = xr.where(cond, 1, 1)
            onesum = (ones*wghts).groupby('time.year').sum(dim='time')     

        var_am = datsum / onesum      
        return var_am

    #--Note if the NaN's are different in each variable you'll be averaging over
    #-- different times for each variable
    dset = False
    try:
        varnames = list(ds.keys())
        dset = True
    except:
        pass

    if (dset):
        for i, ivar in enumerate(varnames):
            var = ds[ivar]
            var_am = dothecalc(var, skipna=skipna)
            var_am = var_am.rename(ivar)
            if (i == 0):
                ds_am = var_am
            else:
                ds_am = xr.merge([ds_am, var_am])
    else:
        ds_am = dothecalc(ds, skipna=skipna)
        ds_am = ds_am.rename(ds.name)

    return ds_am



