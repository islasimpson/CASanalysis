import xarray as xr
import numpy as np
import xesmf as xe
import pop_tools

def regridpop(ds,grid, destgrid):
    rEarth=6378.1
    popgrid = pop_tools.get_grid(grid)
    ds_grid = POP_grid_for_xesmf(popgrid)
    ny = destgrid.lat.size ; nx = destgrid.lon.size
    destgrid['mask'] = xr.DataArray(np.ones((ny,nx)), dims=['lat','lon'],
                       coords={'lat':destgrid.lat, "lon":destgrid.lon})*\
                       np.cos(np.deg2rad(destgrid.lat))*rEarth**2.
    destgrid['lat'] = destgrid['lat'].assign_attrs({'units':'degrees_north','long_name':'latitude'})
    destgrid['lon'] = destgrid['lon'].assign_attrs({'units':'degrees_east','long_name':'longitude'})
    destgrid = add_grid_bounds(destgrid)

    regridder = xe.Regridder(ds_grid, destgrid, 'conservative', periodic=True)
    try:
        ds = ds.rename({'nlat':'y', 'nlon':'x'})
    except:
        ds = ds.rename({'nj':'y', 'ni':'x'})

    datregrid = regridder(ds.load(), skipna=True, na_thres=0.5)


    return datregrid

#----------define the pop grid for xesmf
def POP_grid_for_xesmf(ds):
    """
    Given a POP grid dataset with [TLAT,TLONG,ULAT,ULONG], modify so that
    dataset is suitable for conservative regridding using xesmf.  Assumes
    mapping from T-grid so [ULAT,ULONG] give lat/lon bounds.
    Obtained from Steve Yeager (Feb 24th 2022)
    """
    ulat = ds.ULAT.rename({'nlat':'y_b','nlon':'x_b'})
    ulon = ds.ULONG.rename({'nlat':'y_b','nlon':'x_b'})
    tlat = ds.TLAT.rename({'nlat':'y','nlon':'x'})
    tlon = ds.TLONG.rename({'nlat':'y','nlon':'x'})
    ulat_b = xr.concat([ulat.isel(x_b=-1),ulat],dim='x_b')
    dy = ulat_b.isel(y_b=1) - ulat_b.isel(y_b=0)
    ulat_b = xr.concat([ulat_b.isel(y_b=0)-dy,ulat_b],dim='y_b')
    ulon_b = xr.concat([ulon.isel(x_b=-1),ulon],dim='x_b')
    ulon_b = xr.concat([ulon_b.isel(y_b=0),ulon_b],dim='y_b')
    ds = ds.assign_coords({'lon_b':ulon_b.transpose(),'lat_b':ulat_b.transpose(),'lon':tlon,'lat':tlat})
    ds = ds.drop(['ULAT','ULONG','TLONG','TLAT'])
    return ds

#----------define the output grid for xesm (also from Steve)
def add_grid_bounds(ds):
    ds = ds.cf.add_bounds(['lon','lat'])
    # Check latitude range
    lb = ds['lat_bounds']
    if ((lb<-90) | (lb>90)).any():
        saveattrs = ds['lat'].attrs
        lb = xr.where(lb>90,90,lb)
        lb = xr.where(lb<-90,-90,lb)
        ds['lat_bounds'] = lb
        ds['lat'] = ds['lat'].assign_attrs(saveattrs)
    return ds


#---------Utilities for dealing with conservative remapping
def construct_grid(lon, lat, lon_b, lat_b):
    """Construct the grid information (lons, lats, and bounds) for input to xesmf when
       doing conservative remapping"""
    if (np.array(lon).ndim == 2):
        grid = {
            'lat': lat,
            'lon': lon,
            'lat_b': lat_b,
            'lon_b': lon_b
        }
    else:
        lon_2d = np.tile(lon[np.newaxis,:], (len(lat),1))
        lat_2d = np.tile(lat, (len(lon),1)).T
        grid = {
            'lat': lat_2d,
            'lon': lon_2d,
            'lat_b': lat_b,
            'lon_b': lon_b
        }
    return grid

def get_bounds(lon,lat):
    """ Obtain the longitude and latitude bounds for the xesmf regridder
        (Needed for conservative remapping)
    """

    # Case of a 2D lon/lat array
    if (np.array(lon).ndim == 2):
        nlat = lon[:,0].size
        nlon = lon[0,:].size

        lon_b = np.zeros([nlat + 1, nlon + 1])
        lat_b = np.zeros([nlat + 1, nlon + 1])

        lon_b[1:,1:-1] = (lon[:,1:] - lon[:,0:-1])/2. + lon[:,0:-1]
        lon_b[1:,0] = lon[:,0] - (lon[:,1] - lon[:,0])/2.
        lon_b[1:,-1] = lon[:,-1] + (lon[:,-1] - lon[:,-2])/2.
        lon_b[0,:] = lon_b[1,:]
        lon_b[-1,:] = lon_b[-2,:]

        lat_b[1:-1,1:] = lat[0:-1,:] + (lat[1:,:] - lat[0:-1,:])/2.
        lat_b[0,1:] = lat[0,:] - (lat[1,:] - lat[0,:])/2.
        lat_b[-1,1:] = lat[-1,:] + (lat[-1,:] - lat[-2,:])/2.
        lat_b[:,0] = lat_b[:,1]
        lat_b[:,-1] = lat_b[:,-2]

    else:
        nlat = lat.size
        nlon = lon.size

        lon_b = np.zeros( [nlon+1] )
        lat_b = np.zeros( [nlat+1] )

        lon_b[1:-1] = (lon[1:] - lon[0:-1])/2. + lon[0:-1]
        lon_b[0] = lon[0] - (lon[1]-lon[0])/2.
        lon_b[-1] = lon[-1] + (lon[-1] - lon[-2])/2.
        lon_b = np.tile(lon_b, (nlat+1,1))

        lat_b[1:-1] = lat[0:-1] + (lat[1:] - lat[0:-1])/2.
        lat_b[0] = lat[0] - (lat[1] - lat[0])/2.
        lat_b[-1] = lat[-1] + (lat[-1] - lat[-2])/2.
        lat_b = np.tile(lat_b[:,np.newaxis], (1, nlon+1))

    return lon_b, lat_b

def regrid_conservative(dat, lonin, latin, lonout, latout, reuse_wgts=False, wgtfile='wgt.nc'):
    """ Conservatively remap """
    lonin_b, latin_b = get_bounds(np.array(lonin), np.array(latin))
    grid_in = construct_grid(np.array(lonin), np.array(latin), np.array(lonin_b), np.array(latin_b))

    lonout_b, latout_b = get_bounds(np.array(lonout), np.array(latout))
    grid_out = construct_grid(np.array(lonout), np.array(latout), np.array(lonout_b), np.array(latout_b))

    regridder = xe.Regridder(grid_in, grid_out, 'conservative', reuse_weights=reuse_wgts, filename=wgtfile)

    dat_rg = regridder(dat)
    lon1d = dat_rg.lon.isel(y=0)
    lat1d = dat_rg.lat.isel(x=0)

    dat_rg = dat_rg.drop_vars(['lon','lat'])    
    dat_rg = dat_rg.assign_coords({'x':lon1d, 'y':lat1d})
    dat_rg = dat_rg.rename({'x':'lon', 'y':'lat'})

    return dat_rg

def regrid_bilinear(dat, lonout, latout, reuse_wgts=False, wgtfile='wgt.nc'):
    """ Bilinear remap """
    grid_out = xr.Dataset({'lat': (['lat'], latout.values)}, {'lon': (['lon'], lonout.values)})
    regridder = xe.Regridder(dat, grid_out, 'bilinear', periodic=True, reuse_weights=reuse_wgts,
                     filename=wgtfile)
    dat_rg = regridder(dat)
    return dat_rg 



