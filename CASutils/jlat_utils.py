import numpy as np

def calcjetlat( uzm, minlat, maxlat):
# -----calculate the latitude of a maximum between latitude bounds 
# input: uzm = data array(nlat)
#        minlat = minimum latitude over which to search for the max
#        maxlat = maximum latitude over which to search for the max
# output: jlatv = the jet latitude
#         jmaxv = the jet maximum
# the jet is as the maximum of the quadratic fit to the grid point maximum
# and the two adjacent grid points (as Kidston and Gerber 2010)
# NaN's are skipped
  lats = uzm.sel(lat = slice(minlat,maxlat)).coords['lat']
  imax = uzm.sel(lat = slice(minlat,maxlat)).argmax(dim='lat', skipna=True, keep_attrs=True)
  if (imax == 0) or (imax == len(lats)-1):
    jlatv = np.nan
    jspeedv = np.nan
    print( "!!!! no local maximum found in calcjetlat" )

  else:
    lat4fit = lats.isel(lat=slice(int(imax)-1,int(imax)+2))
    u4fit = uzm.sel(lat = slice(minlat,maxlat)).isel(lat=slice(int(imax)-1,int(imax)+2))
    coefs = np.polyfit(lat4fit,u4fit,2)
    jlatv = -1.*coefs[1]/(2.*coefs[0])
    jspeedv = coefs[2] + coefs[1]*jlatv + coefs[0]*jlatv**2

  return jlatv, jspeedv
