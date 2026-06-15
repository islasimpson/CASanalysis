# Functions for dealing with land model data
import xarray as xr
import numpy as np
import sys


def get_pft_aggregate_lonlat(dat, lon, lat, var, outlevel="gridcell", useactive=True):
    """Aggregate variable var at location lon and lat onto either the   
       gridcell or the landunit level"""

    # Mask to only retain the vegetated land unit if outlevel == "landunit"
    if outlevel == "landunit":
        if "pfts1d_itype_lunit" in dat:
            mask = (dat.pfts1d_itype_lunit == 1).compute()
        if "pfts1d_ityplun" in dat:
            mask = (dat.pfts1d_ityplun == 1).compute()
        dat = dat.isel(pft=mask)

    # Find the locations at lon and lat
    dist2 = (dat.pfts1d_lon - lon)**2. + (dat.pfts1d_lat - lat)**2
    index = np.argmin(dist2.values) 
    lon_near = dat.pfts1d_lon.isel(pft=index)
    lat_near = dat.pfts1d_lat.isel(pft=index)
    locmask = (
    xr.DataArray(np.isclose(dat.pfts1d_lon.compute(), lon_near), dims="pft") &
    xr.DataArray(np.isclose(dat.pfts1d_lat.compute(), lat_near), dims="pft"))
              
    dat_loc = dat.isel(pft=locmask)

    if useactive:
        active = dat_loc.pfts1d_active == 1
    else:
        active = xr.ones_like(dat_loc.pfts1d_active, dtype=bool)

    if outlevel == "landunit":
        w = dat_loc.pfts1d_wtlnd
        dat_out = (dat_loc[var]*w.where(active) ).sum('pft') / \
                  (w.where( dat_loc.pfts1d_active)).sum('pft')

    elif outlevel == "gridcell":
        if "pfts1d_wtgcell" in dat_loc:
            w = dat_loc.pfts1d_wtgcell
        elif "pfts1d_wtxy" in dat_loc:
            w = dat_loc.pfts1d_wtxy
        dat_out = (dat_loc[var]*w.where( active ) ).sum('pft')

    return dat_out




def get_veglunit_lonlat(dat, lon, lat):
    """ Pick out the vegetated landunit at location lon and lat from CLM data on the landunit level """

    # Pick out the vegetated 
    mask = (dat.land1d_ityplunit == 1).compute()
    dat = dat.isel(landunit = mask)

    # Find the landunit nearest to the desired longitude and latitude
    dist2 = ((dat.land1d_lon - lon)**2. + (dat.land1d_lat - lat)**2.)
    index = int(dist2.argmin().values)

    # Select that location
    dat_loc = dat.isel(landunit = index)
    return dat_loc

def get_gridcell_lonlat(dat, lon, lat):
    """ Pick out the land grid cell that's nearest to location lon and lat from CLM gridcell level """
    dist2 = ( (dat.lon - lon)**2. + (dat.lat - lat)**2.)
    index = int(dist2.argmin().values)

    dat_loc = dat.isel(lndgrid = index)
    return dat_loc

