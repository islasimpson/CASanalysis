import xarray as xr
import numpy as np

#------------------calculate vapor pressure from specific humidity and surface pressure
#https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
def calcvpfromhuss(huss, ps):
    """Calculate vapor pressure (in hPa) from specific humidity (in kg/kg) and 
       surface pressure (in Pa) i
    https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    """
    e = (huss*ps)/(0.622 + 0.378*huss)
    e = e/100.
    e = e.rename('vp')
    return e
#---------------------------------------------------------------------------------------

#------------------calculate saturation vapor pressure from temperature 
#                  (or vapor pressure from dew point temperature)
#  Based on Bolton (1980) The computation of equivalent potential temperature, MWR
#-------------------------------------------------------------------------------------
def calcsvp(TK):
    """calculate the saturation vapor pressure (in hPa) from T (in K)"""
    T = TK-273.15
    e = 6.112*np.exp( (17.67*T)/(T + 243.5))
    e = e.rename('svp')
    return e

# -----------------calculate the saturation specific humidity from temperature and pressure
#                  (or specific humidity from dew point temperature)
#  Based on Bolton (1980) The computation of equivalent potential temperature, MWR
def calcsq(TK, PS):
    """calculate the saturation specific humidity (in kg/kg) from T (in K) and PS (in Pa)"""
    T = TK-273.15
    e = 6.112*np.exp( (17.67*T)/(T + 243.5))
    q = (0.622*e)/(PS/100. - (0.378*e))
    q = q.rename('sq')
    return q
