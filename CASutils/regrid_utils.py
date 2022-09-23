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



