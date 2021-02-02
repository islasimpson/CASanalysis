import matplotlib.pyplot as plt
import numpy as np
from CASutils import colormap_utils as mycolors

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker


def contourmap_bothcontinents_fill_nh_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, lonmid=0,labels=True, cmap="blue2red"):
    """ plot a contour map of 2D data dat with coordinates lon and lat
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
     
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-180,180,0,90], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-180, -120, -60, 0,60,120, 180], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['180W','120W','60W','0','60E','120E','180E'], fontsize=12)
        ax.set_yticks([0,30,60,90], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['0','30N','60N','90N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    return ax

def contourmap_bothcontinents_scatter_nh_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, lonmid=0,labels=True, cmap="blue2red"):
    """ plot a contour map of 2D data dat with coordinates lon and lat
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
     
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-180,180,0,90], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-180, -120, -60, 0,60,120, 180], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['180W','120W','60W','0','60E','120E','180E'], fontsize=12)
        ax.set_yticks([0,30,60,90], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['0','30N','60N','90N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap)
    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax





