## routines for calculating seasonal climatology and seasonal timeseries
import xarray as xr
import numpy as np

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

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


def season_mean(ds, var, season = "all", cal = "none"):
    """ calculate climatological mean by season
    Args: ds (xarray.Dataset): dataset
          var (str): variable to use
          season (str): "all", 'DJF', "MAM", "JJA", "SON"
          cal (str): "none"(default) or calendar used for weighting months by number of days
    """
    ## no weighting of months: 
    if cal == "none":
        if season == "all":
            ## calculate mean for all season
            smean = ds[var].groupby('time.season').mean('time')
        else :
            ## calculate mean for specified season
            smean = ds[var].where(ds['time.season'] == season).mean('time')

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
            smean = (ds[var] * weights).groupby('time.season').mean('time')
        else :
            ## calculate weighted mean for specified season
            smean = (ds[var] * weights).where(ds['time.season'] == season).mean('time')

        return smean

def season_ts(ds, var, season):
    """ calculate timeseries of seasonal averages
    Args: ds (xarray.Dataset): dataset
          var (str): variable to calculate 
          season (str): 'DJF', 'MAM', 'JJA', 'SON'
    """
    ## set months outside of season to nan
    ds_season = ds.where(ds['time.season'] == season)

    # calculate 3month rolling mean (only middle months of season will have non-nan values)
    ds_season = ds_season[var].rolling(min_periods=3, center=True, time=3).mean().dropna("time", how='all')

# IRS - this was the old way.  I think the above is more efficient (just using dropna)
#    # reduce to one value per year
#    ds_season = ds_season.groupby('time.year').mean('time')
#
#    # remove years that are nan at all grid points
#    # to get rid of seasons that overlap with the start
#    # or end of the record
#
#    for iyear in ds_season['year'] :
#        if (np.isnan(ds_season.sel(year = iyear)).all().values) :
#             ds_season = ds_season.drop_sel( year = [iyear] )
#
    return ds_season

