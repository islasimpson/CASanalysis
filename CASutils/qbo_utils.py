import matplotlib.pyplot as plt
import numpy as np
from dycoreutils import colormap_utils as mycolors
import sys
from scipy.ndimage import label
from math import nan

def plotqbowinds(fig, data, time, pre, ci, cmin, cmax, titlestr, x1=None, x2=None, y1=None, y2=None, plevvar='ilev', ylim=None):
    """
    Plots a QBO time series as a function of time and log(pressure) 
    """

    data = data.transpose(plevvar,"time")

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.blue2red_cmap(nlevs)

    plt.rcParams['font.size'] = '12'

    if (x1):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    else:
        ax = fig.add_axes()

    ax.contourf(time,-1.*np.log10(pre),data, levels=clevs, cmap=mymap, extend='both')
    ax.set_ylim(-np.log10(100.),-np.log10(1))
    ax.set_yticks([-np.log10(100),-np.log10(30),-np.log10(10),
                   -np.log10(3),-np.log10(1)])
    ax.set_yticklabels(['100','30','10','3','1'])
    ax.set_ylabel('Pressure (hPa)')
    ax.set_title(titlestr, fontsize=16)


    return ax

def plotddamp(fig, data, pre, expname, x1=None, x2=None, y1=None, y2=None, color=None, oplot=False, ax=None):
    """ 
    Plot up the Dunkerton and Delisi amplitude of the QBO.
    Inputs:
        fig = the figure page
        data = the dunkerton and delisi amplitude data
        pre = the pressure axis of data
        expname = the name of the experiment (for legend)
        x1 = the bottom edge of the figure (in units of fractions of the page)
        x2 = the right edge of the figure (in units of fraction of the page)
        y1 = the bottom edge of the figure (in units of fractions of the page)
        y2 = the top edge of the figure ( in units of fractions of the page)
        oplot = if True, only over plot a line
    """

    # if overplotting, check for axis input
    if (oplot and (not ax)):
        print("This isn't going to work.  If overplotting, specify axis")
        sys.exit()

    plt.rcParams['font.size'] = '14'

    if not oplot:
        if (x1):
            ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
        else:
            ax = fig.add_axes()

        ax.set_ylim(-np.log10(100),-np.log10(3))
        ax.set_yticks([-np.log10(100),-np.log10(30),-np.log10(10),-np.log10(3)])
        ax.set_yticklabels(['100','30','10','3'])
        ax.set_ylabel('Pressure (hPa)', fontsize=16)
        ax.set_xlabel('Dunkerton and Delisi amplitude (ms$^{-1}$)',fontsize=16)
        ax.set_title('QBO amplitude', fontsize=16)


    if (color):
        ax.plot(np.array(data),-1.*np.log10(np.array(pre)),linewidth=3,label=expname, color=color)
    else:
        ax.plot(np.array(data),-1.*np.log10(np.array(pre)),linewidth=3,label=expname)

    return ax

def finde2w(dat):
    """
    Find the time indices at which the time series transitions from negative to positive
    i.e., easterly to westerly transitions
    First finding location of the minima, ensuring that they are at least 0.5 standard
    deviations below the mean.  Then finding the time at which the winds next turn
    positive
    """
    
    #---Find the times of the minima
    testdat = dat.where( dat < 0, 0)
    testlabel, testcount = label(testdat)
    nmins = testcount
    eastpeak_time = np.zeros([nmins])
    eastpeak_mag = np.zeros([nmins])
    for i in np.arange(0,nmins,1):
        imin = int(np.argmin(np.where( testlabel == i+1, dat, 1000)))
        eastpeak_time[i] = imin
        eastpeak_mag[i] = dat[imin].values

    #---Make sure the minimum wind values are more than 0.5 standard deviations
    #   below the mean.
    testval = (np.mean(dat) - 0.5*np.std(dat)).values
    eastpeak_time = eastpeak_time[ eastpeak_mag < testval ] 
    eastpeak_mag = eastpeak_mag[ eastpeak_mag < testval ]

    #---Double checking we don't have any peaks that are below 1/3 of the magnitude of the maximum
    maxmageastpeak = np.min(eastpeak_mag)
    eastpeak_time = eastpeak_time[ eastpeak_mag < maxmageastpeak/3. ]
    eastpeak_mag = eastpeak_mag[ eastpeak_mag < maxmageastpeak/3. ]

    #---Now find the transition month to westerlies
    times = np.arange(0,testdat.size,1)
    testdat = dat.where( dat > 0, 0)
    testlabel, testcount = label(testdat)
    
    #find the minimum of each positive label
    minpostime = np.zeros([testcount])
    for i in np.arange(1,testcount+1,1):
        timestest = times[( testlabel == i) ]
        minpostime[i-1] = np.min(timestest)

    transition_time = np.zeros([len(eastpeak_time)])
    for i in np.arange(0,len(eastpeak_time),1):
        timeanoms = minpostime - eastpeak_time[i]
        if (max(timeanoms) < 0):
            transition_time[i] = nan
        else:
            minpostime = minpostime[timeanoms > 0]
            timeanoms = np.where( timeanoms > 0)
            transition_time[i] = minpostime[np.argmin(np.abs(minpostime) - eastpeak_time[i])]
   
    transition_time = transition_time[~np.isnan(transition_time)]

    return transition_time






