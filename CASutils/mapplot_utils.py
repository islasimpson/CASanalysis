import matplotlib.pyplot as plt
import numpy as np
from CASutils import colormap_utils as mycolors

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker


def contourmap_bothcontinents_fill_nh_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
     
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    #ax.cmap.set_over(mymap(len(mymap)-1))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-180,180,0,90], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-180, -120, -60, 0,60,120, 180], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['180W','120W','60W','0','60E','120E','180E'], fontsize=fontsize-3)
        ax.set_yticks([0,30,60,90], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['0','30N','60N','90N'], fontsize=fontsize-3)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="max")

    return ax

def contourmap_bothoceans_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
     
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=210))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
   # ax.set_extent([-180,180,0,90], crs = ccrs.PlateCarree())

   # if (labels):
        #ax.set_xticks([-180, -120, -60, 0,60,120, 180], crs = ccrs.PlateCarree())
        #ax.set_xticklabels(['180W','120W','60W','0','60E','120E','180E'], fontsize=fontsize-3)
        #ax.set_yticks([0,30,60,90], crs = ccrs.PlateCarree())
        #ax.set_yticklabels(['0','30N','60N','90N'], fontsize=fontsize-3)
        #ax.xformatter = LongitudeFormatter()
        #ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="max", transform=ccrs.PlateCarree())

    return ax


def contourmap_bothcontinents_scatter_nh_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red"):
    """ plot a map plot of scatter points for the northern hemisphere with the 
        greenwich meridian at the center.
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)

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

def contourmap_northamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red"):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
     
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-170, -50, 10, 80], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    return ax





