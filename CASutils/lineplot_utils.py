import numpy as np
import matplotlib.pyplot as plt

def plotlinetime_j2j(fig, data, x1, x2, y1, y2, titlestr, yrange=None, yticks=None, yticklabels=None, ytitle=None, linecolor=None):
    """ plot a line plot.  Takes input from jan 1st to dec 31st and plots the line plot from 
    July to June.
    Input: fig = your figure 
           data = a 365 element array containing data to be plotted
           x1 = location of left edge of plot
           x2 = location of right edge of plot
           y1 = location of bottom edge of plot
           y2 = location of top edge of plot
           titlestr = plot title
           yrange = optional range for y axis
           yticks = optional ticks for y axis
           yticklabels = optional tick labels for y axis
           ytitle= optional title for y axis
           linecolor = optional color of line
    """

    july1 = 181
    dataplot = np.zeros([data.size])
    dataplot[0:365-july1]=data[july1:365]
    dataplot[365-july1:365]=data[0:july1]

    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])
    monticks=[0,31,62,92,123,154,185,213,244,274,304,334,365]
    monticks2=np.zeros(12)
    for i in range(0,12):
        monticks2[i] = monticks[i] + (monticks[i+1]-monticks[i])/2.

    if (yrange):
        ax.set_ylim(yrange)

    if (yticks):
        ax.set_yticks(yticks)
       
    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=14)

    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)

    ax.set_xlim([0,365])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2, minor=True)
    ax.set_xticklabels(['J','A','S','O','N','D','J','F','M','A','M','J'], minor=True, fontsize=14)
    ax.set_title(titlestr, fontsize=16)

    if (linecolor):
        ax.plot(np.arange(0,365,1),dataplot, color=linecolor, linewidth=2) 
    else:
        ax.plot(np.arange(0,365,1),dataplot, linewidth=2)
 
    return ax

def oplotlinetime_j2j(ax, data, linecolor=None):
    """ over plot a line on a plot already created using plotlinetime_j2j"""
    july1 = 181
    dataplot = np.zeros([data.size])
    dataplot[0:365-july1]=data[july1:365]
    dataplot[365-july1:365]=data[0:july1]

    if (linecolor):
        ax.plot(np.arange(0,365,1),dataplot, color=linecolor, linewidth=2)    
    else:
        ax.plot(np.arange(0,365,1),dataplot, linewidth=2)
 
    return ax
