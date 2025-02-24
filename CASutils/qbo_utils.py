import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from CASutils import colormap_utils as mycolors
import sys
from scipy.ndimage import label
from math import nan
import matplotlib.colors as colors

def plotqbowinds(fig, data, time, pre, ci, cmin, cmax, titlestr, x1=None, x2=None, y1=None, y2=None, plevvar='ilev', ylim=None, speclevs=None, ylabel=True, contourlines=False, contourlinescale=1, xlabel=None, fsize=12):
    """
    Plots a QBO time series as a function of time and log(pressure) 
    """

    data = data.transpose(plevvar,"time")

    # set up contour levels and color map
    if (speclevs):
        clevs = speclevs
        nlevs = len(clevs)
    else:
        nlevs = (cmax-cmin)/ci + 1
        clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.blue2red_cmap(nlevs)

    plt.rcParams['font.size'] = str(fsize) 

    if (x1):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    else:
        ax = fig.add_axes()

    if (speclevs):
        norm = colors.BoundaryNorm(boundaries=clevs, ncolors=256)
        ax.contourf(time,-1.*np.log10(pre),data, levels=clevs, cmap=mymap, extend='both', norm=norm)

        if (contourlines):
            clevs_lines = np.array(clevs)
            test = np.arange(0,len(clevs_lines),1) - len(clevs)/2 + 0.5
            clevs_lines = clevs_lines[ clevs_lines != 0 ]
            test = test[ test != 0 ]
            test2 = (test/2.).astype(int)
            clevs_lines = clevs_lines[ test2*2 == test ]
            ax.contour(time,-1.*np.log10(pre),data,levels=clevs_lines,colors='black')

    else:
        ax.contourf(time,-1.*np.log10(pre),data, levels=clevs, cmap=mymap, extend='both')

        if (contourlines):
            clevs_lines = clevs*contourlinescale
            clevs_lines = clevs_lines[np.abs(clevs_lines) > ci/2]
            ax.contour(time,-1.*np.log10(pre), data, levels=clevs_lines, colors='black')



    ax.set_ylim(-np.log10(100.),-np.log10(3))
    ax.set_yticks([-np.log10(100),-np.log10(30),-np.log10(10),
                   -np.log10(3)])
    if (ylabel):
        ax.set_ylabel('Pressure (hPa)')
        ax.set_yticklabels(['100','30','10','3'])
    else:
        ax.set_yticklabels([' ',' ',' ',' '])

    if (xlabel):
        ax.set_xlabel(xlabel)

    ax.set_title(titlestr, fontsize=fsize+2)


    return ax

def plotqbowinds_line(fig, data, time, titlestr, x1=None, x2=None, y1=None, y2=None, ylim=None, oplot=False, linecolor=None):
    """
    Plots a QBO time series as a function of time and log(pressure) 
    """

    if (~oplot):
        if (x1):
            ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
        else:
            ax = fig.add_axes()
  
    if (linecolor):
        ax.plot(time, data, color=linecolor, linewidth=2)
    else:
        ax.plot(time, data, linewidth=2)

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
 
    ewloc = transition_time[~np.isnan(transition_time)]
 
 
    #----Now find the westterly to easterly transition times
 
    dat = -1.*dat
    #---Find the times of the maxima
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
 
    weloc = transition_time[~np.isnan(transition_time)]


#    ----Not sure why this is needed.  Commenting out for now.
#    ---logic to deal with the situation when there's more transitions of one type than another
#    ---only works when there's more w-e transitions right now
#    print(weloc)   
#    print(ewloc)
#    
#    if (weloc[0] < ewloc[0]):
#        print('weloc first')
#        idrop=[]
#        for i in np.arange(1,len(weloc),1):
#            dif_we = weloc[i] - weloc[i-1]
#            dif_ew = weloc[i] - ewloc
#            dif_ew = dif_ew[dif_ew > 0]
#            dif_ew = np.min(dif_ew)
#            if (dif_we < dif_ew): # drop the prior w-e
#                idrop.append(i-1)
#        weloc = np.delete(weloc,idrop)
#    else:
#        print('ewloc first')
#        idrop=[]
#        for i in np.arange(0,len(weloc)-1,1):
#            dif_we = weloc[i+1] - weloc[i]
#            dif_ew = ewloc - weloc[i]
#            dif_ew = dif_ew[dif_ew > 0]
#            dif_ew = np.min(dif_ew)
#            if (dif_we < dif_ew): # drop the current w-e
#                idrop.append(i)
#        weloc = np.delete(weloc,idrop) 
    
    #print(idrop)



    return ewloc, weloc 






