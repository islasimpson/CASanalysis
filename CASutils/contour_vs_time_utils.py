import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

from CASutils import colormap_utils as mycolors

def contour_lat_time(fig, dat, time, ci, cmin, cmax, titlestr, x1, x2, y1, y2, xticks=None, yticks=None, xticklabels=None, yticklabels=None, ylabel=None):
    """ Contour lat versus time """
    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    ax.set_title(titlestr, fontsize=15)

    dat = dat.transpose(...,"time") # make the time be the last dim

    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    mymap = mycolors.blue2red_cmap(nlevs)
  
    ax.contourf(time, dat.lat, dat, levels=clevs, cmap=mymap, extend='both')

    if xticks is not None:
        ax.set_xticks(xticks)

    if xticklabels is not None:
        ax.set_xticklabels(xticklabels, fontsize=14)

    if yticks is not None:
        ax.set_yticks(yticks)

    if yticklabels is not None:
        ax.set_yticklabels(yticklabels, fontsize=14)

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14)

    return ax

