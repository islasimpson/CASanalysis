###Subroutines for calculating the 2D blocking statistics of Masato et al (2013) Winter and Summer Northern Hemisphere Blocking in CMIP6 Models, J. Clim.

import importlib
import pandas as pd
import xarray as xr
import numpy as np
from numpy import nan
import sys
import warnings
import math

from scipy.ndimage import label

def getseason_pm5(dat, ystart, yend, season):
    """ Pick out seasons +/- 5 days """
    dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    mondaystr = dat.indexes['time'].strftime('%Y-%m-%d')
    nyears = yend - ystart + 1
    
    if (season == 'DJF'):
        ndays = dpm[11]+dpm[0]+dpm[1]+10
        indices = np.zeros([(nyears-1)*ndays]).astype(int)
    if (season == 'MAM'):
        ndays = dpm[2]+dpm[3]+dpm[4]+10
        indices = np.zeros([(nyears)*ndays]).astype(int)
    if (season == 'JJA'):
        ndays = dpm[5]+dpm[6]+dpm[7]+10
        indices = np.zeros([(nyears)*ndays]).astype(int)
    if (season == 'SON'):
        ndays = dpm[8]+dpm[9]+dpm[10]+10
        indices = np.zeros([(nyears)*ndays]).astype(int)
    
    if (season == 'DJF'):
        for iyear in np.arange(ystart, yend, 1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-12-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear+1)+'-02-28')[0]+5+1
            indices[ (iyear-ystart)*ndays:(iyear-ystart)*ndays+ndays ] = np.arange(index1[0], index2[0],1).astype(int)
        
    if (season == 'MAM'):
        for iyear in np.arange(ystart,yend+1,1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-03-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear)+'-05-31')[0]+5+1
            indices[ (iyear-ystart)*ndays:(iyear-ystart)*ndays+ndays ] = np.arange(index1[0], index2[0], 1).astype(int)
        
    if (season == 'JJA'):
        for iyear in np.arange(ystart,yend+1,1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-06-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear)+'-08-31')[0]+5+1
            indices[ (iyear - ystart)*ndays:(iyear-ystart)*ndays + ndays ] = np.arange(index1[0], index2[0], 1).astype(int)
        
    if (season == 'SON'):
        for iyear in np.arange(ystart,yend+1,1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-09-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear)+'-11-30')[0]+5+1
            indices[ (iyear - ystart)*ndays:(iyear-ystart)*ndays + ndays ] = np.arange(index1[0], index2[0], 1).astype(int)
        
    datout = dat.isel(time = indices)
    
    return datout   


def getseason_pm5_360day(dat, ystart, yend, season):
    """ Pick out seasons +/- 5 days """
    dpm = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
    mondaystr = dat.indexes['time'].strftime('%Y-%m-%d')
    nyears = yend - ystart + 1
    
    if (season == 'DJF'):
        ndays = dpm[11]+dpm[0]+dpm[1]+10
        indices = np.zeros([(nyears-1)*ndays]).astype(int)
    if (season == 'MAM'):
        ndays = dpm[2]+dpm[3]+dpm[4]+10
        indices = np.zeros([(nyears)*ndays]).astype(int)
    if (season == 'JJA'):
        ndays = dpm[5]+dpm[6]+dpm[7]+10
        indices = np.zeros([(nyears)*ndays]).astype(int)
    if (season == 'SON'):
        ndays = dpm[8]+dpm[9]+dpm[10]+10
        indices = np.zeros([(nyears)*ndays]).astype(int)
    
    if (season == 'DJF'):
        for iyear in np.arange(ystart, yend, 1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-12-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear+1)+'-02-30')[0]+5+1
            indices[ (iyear-ystart)*ndays:(iyear-ystart)*ndays+ndays ] = np.arange(index1[0], index2[0],1).astype(int)
        
    if (season == 'MAM'):
        for iyear in np.arange(ystart,yend+1,1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-03-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear)+'-05-30')[0]+5+1
            indices[ (iyear-ystart)*ndays:(iyear-ystart)*ndays+ndays ] = np.arange(index1[0], index2[0], 1).astype(int)
        
    if (season == 'JJA'):
        for iyear in np.arange(ystart,yend+1,1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-06-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear)+'-08-30')[0]+5+1
            indices[ (iyear - ystart)*ndays:(iyear-ystart)*ndays + ndays ] = np.arange(index1[0], index2[0], 1).astype(int)
        
    if (season == 'SON'):
        for iyear in np.arange(ystart,yend+1,1):
            index1 = np.argwhere(mondaystr == str(iyear)+'-09-01')[0]-5
            index2 = np.argwhere(mondaystr == str(iyear)+'-11-30')[0]+5+1
            indices[ (iyear - ystart)*ndays:(iyear-ystart)*ndays + ndays ] = np.arange(index1[0], index2[0], 1).astype(int)
        
    datout = dat.isel(time = indices)
    
    return datout   




def findcontiguous(bday, minlatsearch=40, maxlatsearch=70):
    """Find the local maxima in daily blocking index and grid points that are contiguous to those maxima.
    Input:
        bday[nlat,nlon] = the blocking index for each lon and lat. 
    (Optional):
        minlatsearch = the minimum latitude for searching for blocking maximum
        maxlatsearch = the maximum latitude for searching for blocking maximum
        
    Output:
        blocknumbertemp[nlat,nlon] = array where each grid point is assigned blocking number if it is
            contiguous with the current block at the current timestep. left Nan otherwise
        lonmax = the longitude of the blocking event identified at the current timestep
        latmax = the longitude of the blocking event identified at the current timestep
    """

    blocknumbertemp = np.zeros([bday.lat.size, bday.lon.size])*nan
    
    # zero out negative regions
    bzerod = bday.copy(deep=True)
    bzerod = bzerod.where(bzerod > 0, 0)

    # append 0E-90E to the end of the longitude axis to allow for blocks that span the Greenwich Meridian
    b_0_90 = bzerod.sel(lon=slice(0,90))
    b_0_90.coords['lon'] = (b_0_90.coords['lon'] + 360)

    bwrapped = xr.concat([bzerod, b_0_90], dim='lon')

    # isolate blocks within the search region
    b_minlat_maxlat = bwrapped.sel(lat=slice(minlatsearch, maxlatsearch))
    bnum_minlat_maxlat, countblocks_temp = label(b_minlat_maxlat)
    bnum_minlat_maxlat = xr.DataArray(bnum_minlat_maxlat, coords = b_minlat_maxlat.coords, name='bnum_minlat_maxlat')

    bnum_cp = bnum_minlat_maxlat.copy(deep=True)
    b_cp = b_minlat_maxlat.copy(deep=True)

    # find the location fo the blocking maxima
    lonblock = np.zeros([countblocks_temp]) ; latblock = np.zeros([countblocks_temp])
    countblocks = 0
    for iblock in np.arange(0, countblocks_temp,1):
        btest = b_cp.where(bnum_cp == (iblock+1), -9999)
        indices = btest.argmax(dim=['lat','lon'], keep_attrs = True)
        lonmax = b_cp.lon[indices['lon']]
        latmax = b_cp.lat[indices['lat']]
        if ( (lonmax >= 60) and (lonmax < (360+60)) ):
            lonblock[iblock] = lonmax.values ; latblock[iblock] = latmax.values
        else:
            lonblock[iblock] = nan ; latblock[iblock] = nan

    lonblock = lonblock[~np.isnan(lonblock)] ; latblock = latblock[~np.isnan(latblock)]
    countblocks = lonblock.size
    
    if (countblocks > 0):
    
        # second pass.  Loop over the blocks to find the corresponding number in the larger range and set the contiguous region
        #blocknumbertemp = np.zeros([bday.lat.size, bday.lon.size])*nan
        multiblockflag = np.zeros([bday.lat.size, bday.lon.size])
        for iblock in np.arange(0, countblocks,1):
            buse = bwrapped.copy(deep=True)
    
            # only search for contiguous regions within +/- 60 degrees.  Needed to avoid blocks spreading into other blocks at the highest latitudes
            # may not be necessary now with the multiblock flag below
            buse = buse.where( (buse.lon > (lonblock[iblock]-60)) & (buse.lon < (lonblock[iblock]+60)), 0)
            bnum, countblocks2 = label(buse)
            bnum = xr.DataArray(bnum, coords=bwrapped.coords, name='bnum')
            numofblock = bnum.sel(lon=lonblock[iblock], lat=latblock[iblock])

            blocknumbertemp2 = np.where(bnum == numofblock, iblock+1, nan)
            blocknumbertemp2 = xr.DataArray(blocknumbertemp2, coords=bwrapped.coords)

            # piecing together the required bits when dealing with blocks that can cross the greenwich meridian
            if (lonblock[iblock] >= 300):
                blocknumber_0on = blocknumbertemp2.sel(lon=slice(360,lonblock[iblock]+60))
                blocknumber_0on.coords['lon'] = (blocknumber_0on['lon']-360.)
                maxlon = blocknumber_0on.lon.isel(lon=blocknumber_0on.lon.size-1).values
                minlon = blocknumber_0on.lon.isel(lon=0).values
                test = np.argwhere( np.array(blocknumbertemp2.lon) == maxlon)[0]
                minindex = test + 1

                test = np.argwhere( np.array(blocknumbertemp2.lon) - 360. == minlon)[0]
                maxindex = test

                blocknumber_to360 = blocknumbertemp2.isel(lon=slice(int(minindex), int(maxindex)))
                blocknumbertemp3 = xr.concat([blocknumber_0on, blocknumber_to360], dim='lon')
        
            else:
                blocknumbertemp3 = blocknumbertemp2.isel(lon=slice(0,bday.lon.size))
            
            # flag for special cases where grid points belong to more than one block
            multiblockflag = np.where( (blocknumbertemp3 == iblock+1), multiblockflag+1, multiblockflag)
        
    
            blocknumbertemp = np.where(blocknumbertemp3 == iblock+1, iblock+1, blocknumbertemp)
    
    
        lonblock = np.where(lonblock >= 360, lonblock-360., lonblock)
    
        # fix the grid points that belong to multiple blocks
        londistance = np.zeros([countblocks, bday.lat.size, bday.lon.size])
        lonarray = np.broadcast_to(bday.lon, [bday.lat.size, bday.lon.size])
        for iblock in np.arange(0,countblocks,1):
            londistance[iblock,:,:] = np.abs(lonarray - lonblock[iblock])
    
        test = np.argmin(londistance, axis=0)
        blocknumbertemp = np.where(multiblockflag < 2, blocknumbertemp, test+1)
    
    
    return blocknumbertemp, lonblock, latblock


def calcblocking(dat, ystart, yend):
    """ Identify blocks and their characteristics at the given timestep
    
    Input:
    - dat = an array of the form [ntime, nlat, nlon] containing the geopotential height data for a season over a time period
    ( it is expected that the time axis will only contain the days for the season of interest)
    - ystart = the year in which to start the blocking calculations
    - yend = the year in which to end the blocking calculations 
     (note that for DJF, ystart and yend both refer to the year of the december)
     
    Output:
    - blocknumber[year, day, lat, lon] = contains a number assignment to each block
    - blockstart[number of blocks] = the day of the season that the block starts
    - blockend[number of blocks] = the day of the season that the block ends
    - yearblock[number of blocks] = the year in which the block is found
    - ilonmaxall[40,number of blocks] = the longitude of the maximum for each day of a block (allowing up to 40 days duration)
    - ilatmaxall[40,number of blocks] = the latitude of the maximum for each day of a block (alloing up to 40 days duration)
    
    """
    
    dphi = (( dat.lat.isel(lat=1) - dat.lat.isel(lat=0))/180.)*np.pi # assuming that lat is on a regular grid
    nlat15 = int(np.floor(15./(dat.lat.isel(lat=1) - dat.lat.isel(lat=0)))) # number of latitude spaning 15 deg
    delphi = (30./180.)*np.pi
    
    # latitude range over which b is calculated
    minlatsearch = 40 ; maxlatsearch = 70
    minlatb = 25 ; maxlatb = 75
    
    # blocking index
    b = xr.DataArray(np.zeros([dat.time.size, dat.lat.size, dat.lon.size]),
                coords=[dat.time, dat.lat, dat.lon], 
                dims=['time', 'lat','lon'], name='b')
    indices = np.argwhere( (np.array(dat.lat) >= minlatb) & (np.array(dat.lat) <= maxlatb) )[:,0]
    for latindex in indices:
        zn = (2./delphi)* ( dat.isel(lat=latindex)*dphi/2. + dat.isel(lat=slice(latindex+1,latindex+1+7)).sum("lat")*dphi )
        zs = (2./delphi)* ( dat.isel(lat=latindex)*dphi/2. + dat.isel(lat=slice(latindex-7,latindex)).sum("lat")*dphi )
        b[:,latindex,:] = zn - zs
    
    
    nblocksearch = 30 # maximum number of blocks to search for in a year
    nyears = yend-ystart+1 # number of years for which blocking is being calculated
    ndays = dat.time.size/nyears # the number of days in the season
    

    # initialize output array of blocking number
    blocknumber = xr.DataArray(np.zeros([int(nyears), int(ndays), b.lat.size, b.lon.size])*nan, 
                           coords=[np.arange(0,nyears,1), np.arange(0,ndays,1), b.lat, b.lon], 
                           dims=['year','day','lat','lon'], name='blocknumber')
    
    # initialize arrays that will contain block characteristics
    yearblock = np.empty([0])           # the year in which the block occurs
    blockstart = np.empty([0])          # the day of the season when the block started
    blockend = np.empty([0])            # the day of the season when the block ended
    ilonmaxall = np.zeros([40,1])*nan   # Array that will contain the longitude of the block maximum for each block 
                                        # (1st index = for the different days of the block, 2nd is for the different blocks)
    ilatmaxall = np.zeros([40,1])*nan   # as ilonmaxall but for the latitude of the center of the block
    
    countall = 0 # counter for block numbers
    for iyear in np.arange(ystart, yend+1, 1): # loop over years
        #print('Processing year '+str(iyear))
        byear = b.isel(time=slice(int((iyear-ystart)*ndays), int((iyear-ystart+1)*ndays))) # blocking indices for that year
        
        alreadyblocks = 0 # no blocks have been detected in this season yet
        
        for iday in np.arange(0,byear.time.size,1):
            
            bday = byear.isel(time=iday) # blocking index for the current day
            
            blocknumbertemp, lonblock, latblock = findcontiguous(bday, minlatsearch = 40, maxlatsearch = 70)
            countblocks = len(lonblock)
            newblocks = 0 # initializing counter of new blocks at current timestep
            
            if (countblocks > 0): # blocks are found on this day
                
                if (alreadyblocks > 0): # blocks have already been identified in this season (search for continuations)
                    blockbeforecontinue = np.zeros([ilonblockbefore.size])
                    blockbeforecontinue[:] = 0 # 0 of a block from before is not continuing, 1 if it is
                    blockcontinue = np.zeros([countblocks])
                    blockcontinue[:] = 0 # 0 of the current block is not a continuation, 1 if it is
                    
                    # --------Continuations -----------
                    # loop over blocks from the previous timestep to check for continuations
                    blocktaken = np.zeros([lonblock.size])
                    blocktaken[:] = 0 #set to one if the block is already designated as a continuation

                    for iblock in np.arange(0,ilonblockbefore.size,1): # loop over last timesteps blocks
                        lonblockbefore = bday.lon.isel(lon=ilonblockbefore[iblock])
                        # dealing with blocks around the Greenwich Meridian
                        if ( (lonblockbefore <= 18.) or
                         (lonblockbefore >= (360.-18.)) ): 
                            
                            if (lonblockbefore <= 18):
                                lontest = lonblock.copy()
                                result = np.argwhere(lontest >= (360.-18.))
                                if (result.size > 0): lontest[result] = lontest[result] - 360.
                                dif = np.abs(lontest - np.array(lonblockbefore))
                    
                            if (lonblockbefore >= (360-18)):
                                lontest = lonblock.copy()
                                result = np.argwhere(lontest <= 18)
                                if (result.size > 0): lontest[result] = lontest[result] + 360.
                                dif = np.abs(lontest - np.array(lonblockbefore))
                            
                        else:
                            dif = np.abs(lonblock - np.array(lonblockbefore))
                            
                        minindex = np.argmin(dif) # find the closest block
                        result = dif[minindex]
                        
                        if ((result <= 18) and (blocktaken[minindex] == 0)):
                            if (np.abs( latblock[minindex] - np.array(bday.lat.isel(lat=ilatblockbefore[iblock]))) < 14 ):
                                # check for distance from first occurrence
                                lon1 = lonmaxall[0,blocknumberbefore[iblock].astype(int)-1]
                                lat1 = latmaxall[0,blocknumberbefore[iblock].astype(int)-1]
                                lonnow = lonblock[minindex]
                                latnow = latblock[minindex]
                            
                                dif = np.abs(lonnow - lon1)
                                
                                if ( ((lonnow <= 27.) and (lon1 >= (360.-27.))) ):
                                    dif = lonnow - (lon1-360.)
                                if ( ((lon1 <= 27.) and (lonnow >= (360.-27.))) ):
                                    dif = lon1 - (lonnow -360.)
                            
                                diflat = np.abs(latnow - lat1)
                                

                                if ( (dif <= 27) and (diflat <= 20)):
                                    blockbeforecontinue[iblock] = 1 # previous block is continuing
                                    blockcontinue[minindex] = 1 # this timesteps block is a continuation
                                    getarray = np.where(blocknumbertemp == (minindex+1), blocknumberbefore[iblock], blocknumber[iyear-ystart,iday,:,:])
                                    blocknumber[iyear-ystart,iday,:,:] = getarray
                                
                                    test = np.max(np.argwhere(~np.isnan(lonmaxall[:,blocknumberbefore[iblock].astype(int)-1])))
                                    test = test+1
                                
                                    lonmaxall[test, blocknumberbefore[iblock].astype(int)-1] = lonblock[minindex]
                                    latmaxall[test, blocknumberbefore[iblock].astype(int)-1] = latblock[minindex]
                                    blocktaken[minindex]=1
                                    
                                    
                    # -------new blocks ----------                
                    newblocks = 0
                    for iblock in np.arange(0,countblocks,1):
                        if (blockcontinue[iblock] == 0):
                            blocknumber[iyear-ystart, iday, :,:] = np.where(blocknumbertemp == (iblock+1), countall+newblocks+1, 
                                                                            blocknumber[iyear-ystart,iday,:,:] )
                            blockstart = np.append(blockstart, iday)
                            yearblock = np.append(yearblock, iyear)
                            lonmaxall4append = np.zeros([40,1])*nan
                            latmaxall4append = np.zeros([40,1])*nan
                            lonmaxall4append[0,0] = lonblock[iblock]
                            latmaxall4append[0,0] = latblock[iblock]
                            lonmaxall = np.concatenate( [lonmaxall, lonmaxall4append], axis=1 )
                            latmaxall = np.concatenate( [latmaxall, latmaxall4append], axis=1 )
                            newblocks = newblocks+1
                        
                    blockend = np.append(blockend, np.zeros([newblocks])*nan)
                    
                    # -------ending blocks -----------
                    if (blocknumberbefore.size > 0):
                        result = np.argwhere(blockbeforecontinue == 0)
                        if (result.size > 0):
                            for i in np.arange(0,result.size,1):
                                blockending = blocknumberbefore[result[i]].astype(int)
                                blockend[blockending-1] = iday - 1
                    

                    # -------Set values at current timestep for use in the next ----------
                    ilonblockbefore = np.zeros([countblocks]).astype(int)
                    ilatblockbefore = np.zeros([countblocks]).astype(int)
                    blocknumberbefore = np.zeros([countblocks])
                    for i in np.arange(0, countblocks, 1):
                        ilonblockbefore[i] = np.argmin( np.array(np.abs(bday.lon - lonblock[i]))).astype(int)
                        ilatblockbefore[i] = np.argmin(np.array(np.abs(bday.lat - latblock[i]))).astype(int)
                        blocknumberbefore[i] = blocknumber[iyear-ystart, iday, ilatblockbefore[i], ilonblockbefore[i]]


                else: # first day of the season (not searching for cotinuations)
                    blocknumber[iyear-ystart,iday,:,:] = blocknumbertemp[:,:] + countall
                    blocknumberbefore = np.zeros([countblocks])
                    ilonblockbefore = np.zeros([countblocks]).astype(int)
                    ilatblockbefore = np.zeros([countblocks]).astype(int)
                    for i in np.arange(0,countblocks,1):
                        ilonblockbefore[i] = np.argmin( np.array(np.abs(bday.lon - lonblock[i]))).astype(int)
                        ilatblockbefore[i] = np.argmin( np.array(np.abs(bday.lat - latblock[i]))).astype(int)
                        blocknumberbefore[i] = blocknumber[iyear-ystart, iday, ilatblockbefore[i], ilonblockbefore[i]]
                    yearblock = np.append(yearblock, np.zeros([countblocks])+iyear)
                    blockstart = np.append(blockstart, np.zeros([countblocks]) + iday)
                    blockend = np.append(blockend, np.zeros([countblocks]))

                    
                    if (iyear == ystart):
                        lonmaxall = np.zeros([40,len(lonblock)])*nan
                        latmaxall = np.zeros([40,len(latblock)])*nan
                        lonmaxall[0,:] = lonblock
                        latmaxall[0,:] = latblock
                    else:
                        for iblock in np.arange(0,countblocks,1):
                            lonmaxall4append = np.zeros([40,1])*nan
                            latmaxall4append = np.zeros([40,1])*nan
                            lonmaxall4append[0,0] = lonblock[iblock]
                            latmaxall4append[0,0] = latblock[iblock]
                            lonmaxall = np.concatenate( [lonmaxall, lonmaxall4append], axis=1 )
                            latmaxall = np.concatenate( [latmaxall, latmaxall4append], axis=1 )
                
                    newblocks = countblocks
                    alreadyblocks = countblocks
            
            # add new blocks to the counter
            countall = countall + newblocks
            
    return blocknumber, blockstart, blockend, yearblock, lonmaxall, latmaxall
