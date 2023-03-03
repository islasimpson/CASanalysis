# Utilities to calculate SSW dates
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import label
from math import nan
import pdb

def getseason_ndjfma(dat):
    """
    pull out the November to April seasons omitting the J-A of the first year and 
    the N-D of the last year
    """
    # only using November to April
    datestrings= xr.DataArray(dat.indexes['time'].strftime('%m-%d'), coords = dat.time.coords, name="monstr")
    ybeg = dat.time.where(datestrings == '01-01', drop=True).dt.year[0]
    yend = dat.time.where(datestrings == '12-31', drop=True).dt.year
    yend = yend[yend.size-1]
    
    datseas = dat.where( (dat.time.dt.month == 11) | 
                         (dat.time.dt.month == 12) | 
                         (dat.time.dt.month == 1) | 
                         (dat.time.dt.month == 2) | 
                         (dat.time.dt.month == 3) | 
                         (dat.time.dt.month == 4), 0 )
    
    # omitting Jan, Feb, Mar, April of year 1
    datseas = datseas.where(~( ((dat.time.dt.month == 1) & (dat.time.dt.year == ybeg)) | 
                               ((dat.time.dt.month == 2) & (dat.time.dt.year == ybeg)) | 
                               ((dat.time.dt.month == 3) & (dat.time.dt.year == ybeg)) | 
                               ((dat.time.dt.month == 4) & (dat.time.dt.year == ybeg)) ), 0 )
    
    datseas = datseas.where( ~ ( ((dat.time.dt.month == 11) & (dat.time.dt.year == yend)) | 
                                 ((dat.time.dt.month == 12) & (dat.time.dt.year == yend)) ), 0 ) 
 
    nwinters = yend-ybeg+1-1
   
    return datseas,nwinters

def ssw_cp(dat):
    """
    Obtain the SSW dates following the Charlton and Polvani criterion
    Input: dat = 10hPa, 60N, daily zonal mean zonal wind.
    """

    # get the November - April season
    datseas, nwinters = getseason_ndjfma(dat)

    # pick out times when the u_10_60 goes negative
    uneg = datseas.copy(deep=True)
    uneg = uneg.where( uneg < 0, 0)
    sswnum, count = label(uneg)

    sswnum = sswnum[:,0]
#    sswnum = sswnum[:]
    sswnum = xr.DataArray(sswnum, coords=[datseas.time.values], dims=['time'], name='sswnum')

    dset = xr.merge([datseas, sswnum])

    # first pass at SSW dates (all segments where U goes negative)
    datessw_temp = []
    for issw in np.arange(0,count,1):
        thisone = dset.where( sswnum == (issw + 1), drop=True)
        centraldate = thisone.time.isel(time=0)
        datessw_temp.append(centraldate)

    # remove dates that are in april
    datessw = []
    for issw in np.arange(0,len(datessw_temp),1):
        if (datessw_temp[issw].dt.month != 4):
            datessw.append(datessw_temp[issw])


    # remove events where there aren't 10 consecutive days of westerlies again before April 30th.
    datesswnew = []
    for issw in np.arange(0,len(datessw),1):
        startdate = str(datessw[issw].dt.year.values).zfill(4)+'-'+str(datessw[issw].dt.month.values).zfill(2)+\
          '-'+str(datessw[issw].dt.day.values).zfill(2)
        if ((datessw[issw].dt.month.values == 11) | (datessw[issw].dt.month.values == 12)):
            enddate = str(datessw[issw].dt.year.values + 1).zfill(4)+"-04-30"
        else:
            enddate = str(datessw[issw].dt.year.values).zfill(4)+"-04-30"
        test = datseas.sel(time=slice(startdate,enddate))
       
        # pick out slices that are westerly between SSW and April 30th
        test = test.where(test > 0, 0)
        westerly, countwesterly = label(test)
        westerly = westerly[:,0]
        #westerly = westerly[:]
        westerly = xr.DataArray(westerly, coords=[test.time], dims=['time'], name='westerly')
        flagok = 0
        for iwesterly in np.arange(0,countwesterly,1):
            thiswesterly = westerly.where( westerly == iwesterly + 1, drop=True)
            ndays = thiswesterly.size
            if (ndays >= 10):
                flagok = 1
        
        # drop the date if there are no periods of westerlies > 10 days before April 30th.
        if (flagok == 1):
            datesswnew.append(datessw[issw])


    # remove events that aren't separated from the last event by mre than 20 days of westerlies
    sswyears=np.zeros([len(datesswnew)])
    for issw in np.arange(0,len(datesswnew),1):
        sswyears[issw] = datesswnew[issw].dt.year.values
    
    datesswnew2=[]
    countssw=0
    for issw in np.arange(0,len(datesswnew),1):

        # definitely including the first warming
        if (issw == 0):
            datesswnew2.append(datesswnew[issw])
            countssw = countssw + 1
        else:
            startdate = str(datesswnew2[countssw-1].dt.year.values).zfill(4)+'-'+\
             str(datessw[countssw-1].dt.month.values).zfill(2)+'-'+str(datessw[countssw-1].dt.day.values).zfill(2)
            enddate = str(datesswnew[issw].dt.year.values).zfill(4)+'-'+str(datesswnew[issw].dt.month.values).zfill(2)+\
             '-'+str(datesswnew[issw].dt.day.values).zfill(2)
            sincelast = datseas.sel(time=slice(startdate, enddate))

            if (sincelast.time.size > 200):
                datesswnew2.append(datesswnew[issw])
                countssw = countssw+1
            else:
            # now testing for westerlies since the last SSW.
                sincelast = sincelast.where(sincelast > 0, 0)
                westerly, countwesterly = label(sincelast)
                westerly = westerly[:,0]
                #westerly = westerly[:]
                westerly = xr.DataArray(westerly, coords=[sincelast.time], dims=['time'], name='westerly')
                flagok = 0
                for iwesterly in np.arange(0,countwesterly,1):
                    thiswesterly = westerly.where( westerly == iwesterly + 1, drop=True)
                    ndays = thiswesterly.size
                    if (ndays >= 20):
                        flagok = 1
 
                if (flagok == 1):
                    datesswnew2.append(datesswnew[issw])
                    countssw = countssw+1

    return datesswnew2, nwinters

