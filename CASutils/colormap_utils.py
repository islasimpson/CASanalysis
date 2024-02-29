import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap ## used to create custom colormaps
import matplotlib.colors as mcolors
import numpy as np

def blue2red_cmap(n, nowhite = False, posonly=False):
    """ combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    """

    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if (nowhite):
        nwhite=0

    if (posonly):
        colorsw = np.ones((nwhite,4))
        colors2 = plt.cm.YlOrRd(np.linspace(0,1,int(n)))
        colors = np.vstack((colorsw,colors2))
    else:
        colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg)))
        colors2 = plt.cm.YlOrRd(np.linspace(0,1, int(npos)))
        colorsw = np.ones((nwhite,4))
        colors = np.vstack((colors1, colorsw, colors2))

    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap 

def red2blue_cmap(n, nowhite = False):
    """ combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals """

    if (int(n/2) == n/2):
        #even number of contours
        nwhite=1
        nneg = n/2
        npos = n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if (nowhite):
        nwhite=0

    colors1 = plt.cm.YlOrRd_r(np.linspace(0.1,1,int(npos)))
    colors2 = plt.cm.Blues(np.linspace(0.1,1,int(nneg)))
    colorsw = np.ones((nwhite,4))

    if (nowhite):
        colors = np.vstack( (colors1, colors2))
    else:
        colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
  
    return mymap



def precip_cmap(n, nowhite=False):
    """ combine two existing color maps to create a diverging color map with white in the middle.
    browns for negative, blues for positive
    n = the number of contour intervals
    """
    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if (nowhite):
        nwhite=0


    colors1 = plt.cm.YlOrBr_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.GnBu(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))

    colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap

def precip_cmap_nowhite(n):
    """ combine two existing color maps to create a diverging color map with white in the middle.
    browns for negative, blues for positive
    n = the number of contour intervals
    """
    if (int(n/2) == n/2):
        # even number of contours
        nneg=n/2
        npos=n/2
    else:
        nneg = (n-1)/2
        npos = (n-1)/2

    colors1 = plt.cm.YlOrBr_r(np.linspace(0,0.8, int(nneg)))
    colors2 = plt.cm.GnBu(np.linspace(0.2,1, int(npos)))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap

def wkcmap(n, nwhite=2):
    """ Color map for the Wheeler and Kiladis plot
        n = the number of contour intervals
    """
    ncolors = n - nwhite
    nc1 = np.int(ncolors/5)
    nc2 = np.int(ncolors/5)
    nc3 = np.int(ncolors/5)
    nc4 = np.int(ncolors/5)
    nc5 = np.int(ncolors - (nc1+nc2+nc3+4))

    colorsw = np.ones((nwhite,4))
    colors1 = plt.cm.YlOrBr(np.linspace(0,0.4,nc1))
    colors2 = plt.cm.YlOrBr(np.linspace(0.44,0.7,nc2))
    colors3 = plt.cm.afmhot_r(np.linspace(0.6,0.8,nc3))
    colors4 = plt.cm.gist_heat_r(np.linspace(0.75,1,nc4))

    #colors = np.vstack((colorsw, colors1 + colors3 + colors4))
    colors = np.vstack((colorsw, colors1, colors2, colors3, colors4))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap



