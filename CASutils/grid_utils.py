import xarray as xr

def fliplon_neg2pos(dat,lonname):
   """flip longitudes if lon goes from -180 to 180"""
   negvals = dat.where(dat[lonname] < 0, drop=True)
   negvals[lonname] = negvals[lonname] + 360.
   posvals = dat.where(dat[lonname] >= 0, drop=True)
   datout = xr.concat([posvals,negvals],lonname)
   return datout 

def fliplon_pos2neg(dat,lonname):
   """flip longitudes if lon goes from 0 to 360"""
   negvals = dat.where(dat[lonname] > 180, drop=True)
   negvals[lonname] = negvals[lonname] - 360.
   posvals = dat.where(dat[lonname] <= 180, drop=True)
   datout = xr.concat([negvals,posvals], lonname)
   return datout
