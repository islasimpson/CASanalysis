import xarray as xr
import numpy as np
import glob
import sys
import dask

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal

import importlib
importlib.reload(lens)

pathout="/project/cas02/islas/SINGLEFORCING/cesm2_le/"
topdir="/project/mojave/cesm2/LENS/atm/month_1/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX','SWCF']
varnames=['AODVIS','TGCLDIWP','TGCLDLWP']

memstr = lens.lens2memnamegen_second50(50)
memstr = memstr[10:50]

for varname in varnames:
    filelist=[sorted(glob.glob(topdir+varname+'/*.BHISTsmbb.*'+imem+'*.nc'))+
              sorted(glob.glob(topdir+varname+'/*.BSSP370smbb.*'+imem+'*.nc')) 
              for imem in memstr ]
    members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time']) for i in filelist]
    for imem in np.arange(0,len(filelist),1):
        print(memstr[imem]+' '+str(members[imem].time.size))

    dat = xr.concat(members, dim='M', join='override')
    timebndavg = np.array(dat.time_bnds.isel(M=0),
    dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
    dat['time'] = timebndavg
    dat = dat[varname]
    
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, "DJF").compute()
    mam = cal.season_ts(dat, "MAM").compute()
    jja = cal.season_ts(dat, "JJA").compute()
    son = cal.season_ts(dat, "SON").compute()

    am.to_netcdf(pathout+"ALL_"+varname+"_am.nc")
    djf.to_netcdf(pathout+"ALL_"+varname+"_djf.nc")
    mam.to_netcdf(pathout+"ALL_"+varname+"_mam.nc")
    jja.to_netcdf(pathout+"ALL_"+varname+"_jja.nc")
    son.to_netcdf(pathout+"ALL_"+varname+"_son.nc")
