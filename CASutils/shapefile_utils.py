import xarray as xr
import numpy as np
import math
import shapefile as shp
import geopandas as gp
import regionmask
from numpy import nan

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
    maskcoords = xr.Dataset({'lat' : (['lat'],dat4mask['lat'])}, {'lon' : (['lon'],dat4mask['lon'])})

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

