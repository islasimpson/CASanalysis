import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from CASutils import colormap_utils as mycolors

import cartopy as cart
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

    if (cmap == "precip"):
        mymap = mycolors.precip_map(nlevs)

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

def contourmap_bothcontinents_fill_nh_pos_oplotscatter(ax, dat, lon, lat, 
    ci, cmin, cmax, cmap="blue2red", edgecolor="None"):
    """Over plot scatter points on a NH map"""

    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_map(nlevs)

    ax.scatter(lon, lat, c=dat, marker="o", edgecolors=edgecolor, 
        vmin=cmin, vmax=cmax, cmap=mymap, transform = ccrs.PlateCarree(), zorder=100)

    return ax



def contourmap_bothcontinents_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
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

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    #ax.cmap.set_over(mymap(len(mymap)-1))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-180,180,-90,90], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-180, -120, -60, 0,60,120, 180], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['180W','120W','60W','0','60E','120E','180E'], fontsize=fontsize-3)
        ax.set_yticks([-90,-60,-30,0,30,60,90], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['90S','60S','30S','0','30N','60N','90N'], fontsize=fontsize-3)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="max")

    return ax



def contourmap_centrallon_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, signifdat=None, centrallon=240,
  contourlines=False, contourlinescale=1):
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
#    clevs[np.abs(clevs) < ci/2.] = 0
#    print(clevs)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=centrallon))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=1000)
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
    
    if ( signifdat is not None ):
        lonsignif = signifdat.lon
        signifdat, lonsignif = add_cyclic_point( signifdat, coord=lonsignif)
        ax.contourf(lon, lat, signifdat, levels=[0,0.5,1], colors='lightgray', transform = ccrs.PlateCarree())

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > (ci/2.)] 
        ax.contour(lon, lat, dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())

    ax.set_global()

    return ax





def contourmap_bothoceans_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, signifdat=None, stipplesignif=False, contourlines=None, contourlinescale=1):
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
#    clevs[np.abs(clevs) < ci/2.] = 0
#    print(clevs)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=240))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)
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
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())
    
    if ( signifdat is not None ):
        lonsignif = signifdat.lon
        signifdat, lonsignif = add_cyclic_point( signifdat, coord=lonsignif)
        if (stipplesignif):
            density=3
            ax.contourf(lonsignif, lat, signifdat, levels=[0,0.5,1], colors='none', 
               hatches=[density*'.',density*'.', density*','],
               transform = ccrs.PlateCarree())
        else:
            ax.contourf(lonsignif, lat, signifdat, levels=[0,0.5,1], colors='lightgray', 
               transform = ccrs.PlateCarree())

    ax.set_global()

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci] 
        ax.contour(lon, lat, dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())

    return ax

def contourmap_bothcontinents_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, maskocean=False):
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

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)
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
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    ax.set_global()

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    return ax


def contourmap_bothcontinents_northstereo_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, maskocean=False, latlim=40., nowhite=False):
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
        mymap = mycolors.blue2red_cmap(nlevs, nowhite)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs, nowhite)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,latlim,90], ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)

    theta = np.linspace(0,2*np.pi,100)
    center, radius = [0.5,0.5],0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts*radius+center)
    ax.set_boundary(circle, transform = ax.transAxes)


    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')


    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    return ax


def contourmap_bothcontinents_southstereo_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, maskocean=False, latlim=-40):
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

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,latlim], ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)

    theta = np.linspace(0,2*np.pi,100)
    center, radius = [0.5,0.5],0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts*radius+center)
    ax.set_boundary(circle, transform = ax.transAxes)


    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    return ax





def contourmap_continentsonly_robinson_noborder_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red",nowhite=False, fontsize=15, maskocean=False, signifdat=None):
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
        mymap = mycolors.blue2red_cmap(nlevs, nowhite=nowhite)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs, nowhite=nowhite)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    ax.set_global()

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    ax.axis('off')
    ax.set_extent([-180,180,-57,90], crs = ccrs.PlateCarree())

    if ( signifdat is not None ):
        lonsignif = signifdat.lon
        signifdat, lonsignif = add_cyclic_point( signifdat, coord=lonsignif)
        ax.contourf(lon, lat, signifdat, levels=[0,0.5,1], colors='lightgray', 
            transform = ccrs.PlateCarree())


    return ax






def contourmap_bothcontinents_robinson_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, maskocean=False, oplot=False, ax=None,
 markersize=10):
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

    if (oplot == False):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
        ax.set_aspect('auto')
        ax.add_feature(cfeature.COASTLINE)

        ax.set_title(titlestr, fontsize=fontsize)

    ax.scatter(lon, lat, c=dat, s=markersize, marker="o", vmin=cmin, vmax=cmax, cmap = mymap, transform=ccrs.PlateCarree(), zorder=100)

    if (~oplot):

        ax.set_global()

        if (maskocean):
            ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

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


def contourmap_northamerica_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, markersize=20, signifdat=None):
    """ plot a map plot of scatter points for the northern america 
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

    ax.scatter(lon, lat, c=dat, marker="o", s=markersize, vmin=cmin, vmax=cmax, cmap = mymap)

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if ( signifdat is not None ):
      ax.scatter(lon, lat, c=signifdat, marker="o", s=markersize, vmin=0, vmax=2, 
                 cmap='gray_r') 
      #ax.scatter(lon, lat, c=signifdat, marker="o", color='gray')

    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax



def contourmap_northatlantic_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, contourlines=False, contourlinescale=1):
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
    ax.add_feature(cfeature.COASTLINE, zorder=1000)
    ax.set_extent([-100, 30, 10, 80], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-90,-60, -30, 0, 30], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['90W','60W','30S','0','30E'], fontsize=12)
        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend='both')

    if (contourlines):
        clevs = clevs[clevs != 0]
        clevs = clevs*contourlinescale
        cpos = ax.contour(lon,lat,dat,levels=clevs[clevs > 0], colors='black')
        cneg = ax.contour(lon,lat,dat,levels=clevs[clevs < 0], colors='black', linestyle='dotted')
        ax.clabel(cpos, cpos.levels, inline=True, fontsize=10, fmt='{:.1f} '.format)
        ax.clabel(cneg, cneg.levels, inline=True, fontsize=10, fmt='{:.1f} '.format)


    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    return ax


def contourmap_med_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, contourlines=False, contourlinescale=1,
 signifdat=None):
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
    ax.set_extent([-40, 60, 10, 80], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-90,-60, -30, 0, 30], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['90W','60W','30S','0','30E'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    if (contourlines):
        clevs = clevs[np.abs(clevs) > ci/2.]
        clevs = clevs*contourlinescale
        cpos = ax.contour(lon,lat,dat,levels=clevs[clevs > 0], colors='black')
        cneg = ax.contour(lon,lat,dat,levels=clevs[clevs < 0], colors='black', linestyle='dotted')
        ax.clabel(cpos, cpos.levels, inline=True, fontsize=10, fmt='{:.1f} '.format)
        ax.clabel(cneg, cneg.levels, inline=True, fontsize=10, fmt='{:.1f} '.format)

    if ( signifdat is not None ):
        lonsignif = signifdat.lon
        signifdat, lonsignif = add_cyclic_point( signifdat, coord=lonsignif)
        ax.contourf(lon, lat, signifdat, levels=[0,0.5,1], colors='lightgray', transform = ccrs.PlateCarree())

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    return ax





def contourmap_pacificnorthamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False,contourlines=False, contourlinescale=1):
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

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)
    ax.set_extent([100,320, -10, 80], crs = ccrs.PlateCarree(ccrs.PlateCarree(central_longitude=180)))

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()
#
    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, transform=ccrs.PlateCarree(), extend='both')

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2.]
        ax.contour(lon, lat, dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())


    return ax






def contourmap_northamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, contourlines=False, contourlinescale=1):
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

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2. ]
        ax.contour(lon,lat,dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())


    return ax


def contourmap_africa_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
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
    ax.set_extent([-25, 55, -40, 40], crs = ccrs.PlateCarree())

    #if (labels):
     #   ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
     #   ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
     #   ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
     #   ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
     #   ax.xformatter = LongitudeFormatter()
     #   ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    return ax

def contourmap_africa_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red", markersize=20, signifdat=None):
    """ plot a map plot of scatter points for africa
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
    ax.set_extent([-25, 55, -40, 40], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap, s=markersize)

    if ( signifdat is not None ):
      ax.scatter(lon, lat, c=signifdat, marker="o", s=markersize, vmin=0, vmax=2,
                 cmap='gray_r')


    #ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap)
    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax



def contourmap_southamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
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
    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    ax.add_feature(cfeature.COASTLINE, zorder=1000)
    ax.set_extent([-90, -30, -60, 20], crs = ccrs.PlateCarree())

    #if (labels):
     #   ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
     #   ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
     #   ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
     #   ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
     #   ax.xformatter = LongitudeFormatter()
     #   ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    return ax


def contourmap_middleeast_asia_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
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
    ax.set_extent([30, 125, 5, 55], crs = ccrs.PlateCarree())

    #if (labels):
     #   ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
     #   ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
     #   ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
     #   ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
     #   ax.xformatter = LongitudeFormatter()
     #   ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    return ax

def contourmap_middleeast_asia_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red"):
    """ plot a map plot of scatter points for africa
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
    ax.set_extent([30, 125, 5, 55], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap)
    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax

#def contourmap_southamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
# x1, x2, y1, y2, labels=True, cmap="blue2red"):
#    """ plot a contour map of 2D data dat with coordinates lon and lat
#        Input:
#              fig = the figure identifier
#              dat = the data to be plotted
#              lon = the longitude coordinate
#              lat = the latitude coordinate
#              ci = the contour interval
#              cmin = the minimum of the contour range
#              cmax = the maximum of the contour range
#              titlestr = the title of the map
#              x1 = position of the left edge
#              x2 = position of the right edge
#              y1 = position of the bottom edge
#              y2 = position of the top edge
#              labels = True/False (ticks and  labels are plotted if true) 
#              cmap = color map (only set up for blue2red at the moment)
#    """
#
#    # set up contour levels and color map
#    nlevs = (cmax-cmin)/ci + 1
#    clevs = np.arange(cmin, cmax+ci, ci)
#     
#    if (cmap == "blue2red"):
#        mymap = mycolors.blue2red_cmap(nlevs)
#
#    if (cmap == "precip"):
#        mymap = mycolors.precip_cmap(nlevs)
#
#    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
#    ax.set_aspect('auto')
#    ax.add_feature(cfeature.COASTLINE)
#    ax.set_extent([-85, -30, -62, 20], crs = ccrs.PlateCarree())
#
#    #if (labels):
#     #   ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#     #   ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#     #   ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#     #   ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#     #   ax.xformatter = LongitudeFormatter()
#     #   ax.yformatter = LatitudeFormatter()
#
#    ax.set_title(titlestr, fontsize=16)
#
#    dat, lon = add_cyclic_point(dat, coord=lon)
#    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)
#
#    return ax

def contourmap_southamerica_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red"):
    """ plot a map plot of scatter points for africa
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
    ax.set_extent([-85, -30, -62, 20], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap)
    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax

def contourmap_europe_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
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
    ax.set_extent([-30, 50, 35, 70], crs = ccrs.PlateCarree())

    #if (labels):
     #   ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
     #   ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
     #   ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
     #   ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
     #   ax.xformatter = LongitudeFormatter()
     #   ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    return ax

def contourmap_europe_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red"):
    """ plot a map plot of scatter points for africa
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
    ax.set_extent([-30, 50, 35, 70], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap)
    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax





def contourmap_australia_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
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
    ax.set_extent([112, 156, -45, -10], crs = ccrs.PlateCarree())

    #if (labels):
     #   ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
     #   ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
     #   ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
     #   ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
     #   ax.xformatter = LongitudeFormatter()
     #   ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    return ax

def contourmap_australia_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr, 
 x1, x2, y1, y2, labels=True, cmap="blue2red"):
    """ plot a map plot of scatter points for africa
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
    ax.set_extent([112, 156, -45, -10], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", vmin=cmin, vmax=cmax, cmap = mymap)
    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax


