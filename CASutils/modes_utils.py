import xarray as xr

from CASutils import averaging_utils as avg

def calcAMV(sst):
    """ Calculate the Trenberth and Shea AMV index.
        natl = 80W-0, 0-60N average
        glob = 0-360E, 60S-60N average
        AMV = natl/glob
    """
    natl = avg.cosweightlonlat(sst,280,360,0,60)
    glob = avg.cosweightlonlat(sst,0,360,-60,60)
    amv = natl/glob
    return amv

def calcAMOC(MOC, lat=26.5):
    """ Calculate the AMOC as measured by the maximum stream function
    in the Atlantic sector below 500m depth at 26.5 degrees North.
    This requires the input variable MOC from the ocean component.
    If another latitude is required specify lat
    """

    amoc = MOC.sum('moc_comp').sel(lat_aux_grid=lat,method='nearest').\
           isel(transport_reg=1).sel(moc_z=slice(50000.,None)).max('moc_z')

    return amoc
