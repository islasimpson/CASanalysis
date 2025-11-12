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




def get_bounds_from_vertices(lon_vertices, lat_vertices):
    nlat, nlon, _ = lon_vertices.shape
    lon_b = np.zeros((nlat+1, nlon+1))
    lat_b = np.zeros((nlat+1, nlon+1))

    # Find the indices for the different corners
    lon_vertices_np = np.array(lon_vertices[1,1,:]) ; lat_vertices_np = np.array(lat_vertices[1,1,:])
    index_lower_left = np.argmin(np.array(lon_vertices_np) + np.array(lat_vertices_np))
    index_upper_right = np.argmax(np.array(lon_vertices_np) + np.array(lat_vertices_np))
#    index_lower_right = np.argwhere( (lat_vertices_np == lat_vertices_np[index_lower_left]) & (lon_vertices_np != lon_vertices_np[index_lower_left]) )[0][0]
#    index_upper_left = np.argwhere( (lon_vertices_np == lon_vertices_np[index_lower_left]) & (lat_vertices_np != lat_vertices_np[index_lower_left]))[0][0]
    index_lower_right = np.argwhere( (lon_vertices_np > lon_vertices_np[index_lower_left]) & 
                                     (lat_vertices_np < lat_vertices_np[index_upper_right]))[0][0] 
    index_upper_left = np.argwhere( (lon_vertices_np < lon_vertices_np[index_upper_right]) & 
                                     (lat_vertices_np > lat_vertices_np[index_lower_left]))[0][0]


    lon_b[:-1, :-1] = lon_vertices[:,:,index_lower_left] # lower left
    lon_b[1:,:-1] = lon_vertices[:,:,index_upper_left] # upper left
    lon_b[:-1,1:] = lon_vertices[:,:,index_lower_right] # lower right
    lon_b[1:,1:] = lon_vertices[:,:,index_upper_right] # upper right
    
    lat_b[:-1,:-1] = lat_vertices[:,:,index_lower_left]
    lat_b[1:,:-1] = lat_vertices[:,:,index_upper_left]
    lat_b[:-1,1:] = lat_vertices[:,:,index_lower_right]
    lat_b[1:,1:] = lat_vertices[:,:,index_upper_right]
 
#    lon_b[:-1, :-1] = lon_vertices[:, :, 0]  # lower left
#    lon_b[:-1, 1:]  = lon_vertices[:, :, 1]  # lower right
#    lon_b[1:, 1:]   = lon_vertices[:, :, 2]  # upper right
#    lon_b[1:, :-1]  = lon_vertices[:, :, 3]  # upper left
#
#    lat_b[:-1, :-1] = lat_vertices[:, :, 0]
#    lat_b[:-1, 1:]  = lat_vertices[:, :, 1]
#    lat_b[1:, 1:]   = lat_vertices[:, :, 2]
#    lat_b[1:, :-1]  = lat_vertices[:, :, 3]

   
    return lon_b, lat_b


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

def regrid_conservative_fromvertices(dat, lon_in, lat_in, lonvertices_in, latvertices_in, lon_out, lat_out, reuse_weights=False,
                                        wgtfile='wgt.nc'):

    lon_in_temp = np.unwrap(np.deg2rad(lon_in)) * 180/np.pi
    lonvertices_in_temp = np.unwrap(np.deg2rad(lonvertices_in), axis=1) * 180/np.pi

    lon_in = lon_in_temp
    lonvertices_in = lonvertices_in_temp

    """ Conservatively regrid from a curvilinear grid to a regular grid when vertices are provided """
    nlat, nlon, _ = lonvertices_in.shape
    lon_b_in, lat_b_in = get_bounds_from_vertices(lonvertices_in, latvertices_in)
    source_grid = {
       'lon': lon_in,
       'lat': lat_in,
       'lon_b': lon_b_in,
       'lat_b': lat_b_in
       }

    dlon = lon_out[1] - lon_out[0]
    dlat = lat_out[1] - lat_out[0]
    lon_bnds = np.arange(lon_out[0] - dlon/2, lon_out[-1] + 3*dlon/2, dlon)
    lat_bnds = np.arange(lat_out[0] - dlat/2, lat_out[-1] + 3*dlat/2, dlat)

    dest_grid = {
        'lon': lon_out,
        'lat': lat_out,
        'lon_b': lon_bnds,
        'lat_b': lat_bnds
        }

    regridder = xe.Regridder(source_grid, dest_grid, 'conservative',
      periodic=True, reuse_weights=False)

    dat_rg = regridder(dat)
    return dat_rg




def regrid_conservative(dat, lonin, latin, lonout, latout, reuse_wgts=False, wgtfile='wgt.nc', lonin_b=None, latin_b=None):
    """ Conservatively remap """


    if lonin_b is None:
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



