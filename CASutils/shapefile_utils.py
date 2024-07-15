import xarray as xr
import numpy as np
import math
import shapefile as shp
import geopandas as gp
import regionmask
from numpy import nan
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def maskgen(shpfile, dat4mask, regionname):
    """ Generate a mask using information from a shapefile.  Mask will have 1's 
    within the desired region, nan's everywhere else
    Input: 
        shpfile = the shapefile 
        dat4mask = the data that you're planning to mask
        regionname (list) = a list of the region you want to mask.  (assuming this is specified using
         NAME_1 i.e., full name of the state or country ["Alabama", "Alaska"...])
    Output:
        mask = the mask
    """

    # setup of the grid for the mask from dat4mask
    maskcoords = xr.Dataset({'lat' : (['lat'],dat4mask['lat'].values)}, {'lon' : (['lon'],dat4mask['lon'].values)})

    mask = np.zeros([maskcoords.lat.size, maskcoords.lon.size])

    # read in shapefile
    shpcontents = gp.read_file(shpfile)

    # loop over states to mask
    for i in range(0,len(regionname),1):
        print("masking "+regionname[i]) 
        try:
            region = shpcontents[shpcontents.NAME_1 == regionname[i]]
        except:
            region = shpcontents[shpcontents.NAME_0 == regionname[i]]
        maskt = regionmask.mask_geopandas(region, maskcoords["lon"], maskcoords["lat"])
        maskt = np.where(np.isnan(maskt), 0, 1)
        mask[:,:] = mask[:,:] + maskt[:,:]

    # ensure unmasked region is set to 1, rest set to nan's
    mask = np.where(mask == 0, nan, 1)
    mask = xr.DataArray(mask, coords=maskcoords.coords)

    return mask

def plotlower48(ax, shpfile, dat4mask, color='red'):
    shp = gp.read_file(shpfile)
    maskcoords = xr.Dataset({'lat' : (['lat'],dat4mask['lat'].values)}, {'lon' : (['lon'],dat4mask['lon'].values)})
    for istat in np.arange(0,len(shp.NAME_1)):
        print(shp.NAME_1[istat])
        region = shp[shp.NAME_1 == shp.NAME_1[istat]]
        mask = regionmask.mask_geopandas(region, maskcoords["lon"], maskcoords["lat"])
        mask = np.where(np.isnan(mask),0,1)
        mask = xr.DataArray(mask, coords=maskcoords.coords)

        ax = plt.contour(mask.lon, mask.lat, np.nan_to_num(mask), levels=[0,1], colors=color, 
               linewidth=2, transform=ccrs.PlateCarree())

    return ax
