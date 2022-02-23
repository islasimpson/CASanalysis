import importlib
import pandas as pd
import xarray as xr
import numpy as np
from numpy import nan
import sys
import warnings
import math
from glob import glob
import xesmf as xe
import matplotlib.pyplot as plt
from scipy.ndimage import label

from CASutils import blocking_utils as block

def read_data(fpath, datestart, dateend, level):
#First try opening and doing the select assuming everything is working ok with the time axis
    print(fpath)
    try:
        dat = \
        xr.open_mfdataset\
        (fpath, coords="minimal", join="override", decode_times=True, use_cftime=True).\
        sel(time=slice(datestart, dateend))
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
            print("changing longitude --> lon, latitude --> lat")
        except: pass

    except:
        dat = xr.open_mfdataset(fpath, coords="minimal", join="override", decode_times = False)
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
        except: pass
        
        dat = xr.decode_cf(dat, use_cftime = True)
        dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    dat = dat.sel(plev=level, method="nearest")
    return dat


# Input and Output paths
filepath="/project/cas02/islas/CESM2_CMIP6/obs/ERA5/day/zg/" # base path for CMIP6 data
wgtfile="/project/cas/islas/temp/cmip6/wgtfile.nc" # location for weight file for regridder
pathout="/project/cas/islas/python_savs/CESM2eval/blocking/ERA5/" # location for output files

# Start and end year of analysis
ystart=1979 # start year
yend=2014 # end year

# Reading in the data
fpath = filepath+'/*.nc'
#fpath = fpath[0]
#print(fpath)
#sys.exit()

#dat = read_data(fpath+"*.nc", str(ystart)+"-01-01",str(yend)+"-12-31",50000.)

dat = xr.open_mfdataset(fpath)
dat = dat.sel(pre=500.0)
#dat = read_data(fpath+"*.nc", str(ystart)+"-01-01",str(yend)+"-12-30",50000.)    

# Regridding to the 2deg grid

datdjf = block.getseason_pm5(dat.zg, ystart, yend, 'DJF')
#datdjf = datdjf.drop('plev')
datmam = block.getseason_pm5(dat.zg, ystart, yend, 'MAM')
#datmam = datmam.drop('plev')
datjja = block.getseason_pm5(dat.zg, ystart, yend, 'JJA')
#datjja = datjja.drop('plev')
datson = block.getseason_pm5(dat.zg, ystart, yend, 'SON')
#datson = datson.drop('plev')
  
blocknumber_djf, blockstart_djf, blockend_djf, yearblock_djf, lonmaxall_djf, latmaxall_djf = \
       block.calcblocking(datdjf, ystart, yend-1)


blocknumber_mam, blockstart_mam, blockend_mam, yearblock_mam, lonmaxall_mam, latmaxall_mam = \
       block.calcblocking(datmam, ystart, yend)

blocknumber_jja, blockstart_jja, blockend_jja, yearblock_jja, lonmaxall_jja, latmaxall_jja = \
       block.calcblocking(datjja, ystart, yend)

blocknumber_son, blockstart_son, blockend_son, yearblock_son, lonmaxall_son, latmaxall_son = \
       block.calcblocking(datson, ystart, yend)


blockstart_djf = xr.DataArray(blockstart_djf, coords=[np.arange(0,len(blockstart_djf),1)], dims=['block'], name='blockstart')
blockend_djf = xr.DataArray(blockend_djf, coords=[np.arange(0,len(blockend_djf),1)], dims=['block'], name='blockend')
yearblock_djf = xr.DataArray(yearblock_djf, coords=[np.arange(0,len(yearblock_djf),1)], dims=['block'], name='yearblock')
lonmaxall_djf = xr.DataArray(lonmaxall_djf, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_djf))], dims=['day_of_block', 'block'], 
                             name='lonmaxall')
latmaxall_djf = xr.DataArray(latmaxall_djf, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_djf))], dims=['day_of_block', 'block'], 
                             name='latmaxall')

blockstart_mam = xr.DataArray(blockstart_mam, coords=[np.arange(0,len(blockstart_mam),1)], dims=['block'], name='blockstart')
blockend_mam = xr.DataArray(blockend_mam, coords=[np.arange(0,len(blockend_mam),1)], dims=['block'], name='blockend')
yearblock_mam = xr.DataArray(yearblock_mam, coords=[np.arange(0,len(yearblock_mam),1)], dims=['block'], name='yearblock')
lonmaxall_mam = xr.DataArray(lonmaxall_mam, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_mam))], dims=['day_of_block', 'block'], 
                             name='lonmaxall')
latmaxall_mam = xr.DataArray(latmaxall_mam, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_mam))], dims=['day_of_block', 'block'], 
                             name='latmaxall')

blockstart_jja = xr.DataArray(blockstart_jja, coords=[np.arange(0,len(blockstart_jja),1)], dims=['block'], name='blockstart')
blockend_jja = xr.DataArray(blockend_jja, coords=[np.arange(0,len(blockend_jja),1)], dims=['block'], name='blockend')
yearblock_jja = xr.DataArray(yearblock_jja, coords=[np.arange(0,len(yearblock_jja),1)], dims=['block'], name='yearblock')
lonmaxall_jja = xr.DataArray(lonmaxall_jja, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_jja))], dims=['day_of_block', 'block'], 
                             name='lonmaxall')
latmaxall_jja = xr.DataArray(latmaxall_jja, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_jja))], dims=['day_of_block', 'block'], 
                             name='latmaxall')

blockstart_son = xr.DataArray(blockstart_son, coords=[np.arange(0,len(blockstart_son),1)], dims=['block'], name='blockstart')
blockend_son = xr.DataArray(blockend_son, coords=[np.arange(0,len(blockend_son),1)], dims=['block'], name='blockend')
yearblock_son = xr.DataArray(yearblock_son, coords=[np.arange(0,len(yearblock_son),1)], dims=['block'], name='yearblock')
lonmaxall_son = xr.DataArray(lonmaxall_son, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_son))], dims=['day_of_block', 'block'], 
                             name='lonmaxall')
latmaxall_son = xr.DataArray(latmaxall_son, coords=[np.arange(0,40,1), np.arange(0,len(yearblock_son))], dims=['day_of_block', 'block'], 
                             name='latmaxall')

# output the results
blocknumber_djf.to_netcdf(pathout+'ERA5_blocking_djf.nc')
blockstart_djf.to_netcdf(pathout+'ERA5_blocking_djf.nc', mode='a')
blockend_djf.to_netcdf(pathout+'ERA5_blocking_djf.nc', mode='a')
yearblock_djf.to_netcdf(pathout+'ERA5_blocking_djf.nc', mode='a')   
lonmaxall_djf.to_netcdf(pathout+'ERA5_blocking_djf.nc', mode='a')
latmaxall_djf.to_netcdf(pathout+'ERA5_blocking_djf.nc', mode='a')


blocknumber_mam.to_netcdf(pathout+'ERA5_blocking_mam.nc')
blockstart_mam.to_netcdf(pathout+'ERA5_blocking_mam.nc', mode='a')
blockend_mam.to_netcdf(pathout+'ERA5_blocking_mam.nc', mode='a')
yearblock_mam.to_netcdf(pathout+'ERA5_blocking_mam.nc', mode='a')   
lonmaxall_mam.to_netcdf(pathout+'ERA5_blocking_mam.nc', mode='a')
latmaxall_mam.to_netcdf(pathout+'ERA5_blocking_mam.nc', mode='a')

blocknumber_jja.to_netcdf(pathout+'ERA5_blocking_jja.nc')
blockstart_jja.to_netcdf(pathout+'ERA5_blocking_jja.nc', mode='a')
blockend_jja.to_netcdf(pathout+'ERA5_blocking_jja.nc', mode='a')
yearblock_jja.to_netcdf(pathout+'ERA5_blocking_jja.nc', mode='a')   
lonmaxall_jja.to_netcdf(pathout+'ERA5_blocking_jja.nc', mode='a')
latmaxall_jja.to_netcdf(pathout+'ERA5_blocking_jja.nc', mode='a')

blocknumber_son.to_netcdf(pathout+'ERA5_blocking_son.nc')
blockstart_son.to_netcdf(pathout+'ERA5_blocking_son.nc', mode='a')
blockend_son.to_netcdf(pathout+'ERA5_blocking_son.nc', mode='a')
yearblock_son.to_netcdf(pathout+'ERA5_blocking_son.nc', mode='a')   
lonmaxall_son.to_netcdf(pathout+'ERA5_blocking_son.nc', mode='a')
latmaxall_son.to_netcdf(pathout+'ERA5_blocking_son.nc', mode='a')












