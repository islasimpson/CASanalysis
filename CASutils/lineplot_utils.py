import numpy as np
import matplotlib.pyplot as plt

def plotlinetime_j2d_monthly(fig, data, x1, x2, y1, y2, titlestr, yrange=None, 
    yticks=None, yticklabels=None, ytitle=None, linecolor=None, points=True, label=None, linestyle='solid'):
    """ plot a line plot using monthly data from Jan to Dec
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

    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])
 
    monticks = np.arange(0,12,1)
    monticks2 = np.arange(0,12,1)+0.5
    if (yrange):
        ax.set_ylim(yrange)

    if (yticks):
        ax.set_yticks(yticks)

    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=12)

    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)

    ax.set_xlim([0,12])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2, minor=True)
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], minor=True, fontsize=12)
    ax.set_title(titlestr, fontsize=16)

    if linecolor is not None:
            ax.plot(np.arange(0,12,1)+0.5,data,color=linecolor, linewidth=2, label=label, linestyle=linestyle)
            if (points == True):
                ax.plot(np.arange(0,12,1)+0.5,data,"o",markerfacecolor=linecolor, 
                markeredgecolor="black", markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,data, linewidth=2, label=label, linestyle=linestyle)
        if (points == True):
            ax.plot(np.arange(0,12,1)+0.5,data,"o",markeredgecolor="black", 
            markersize=10, markeredgewidth=2)

    return ax

def oplotlinetime_j2d_monthly(ax, data, linecolor=None, linewidth=1, points=True, label=None):
    """ overplot a line on a January - December monthly line plot"""

    if (linecolor[0]):
        ax.plot(np.arange(0,12,1)+0.5,data,color=linecolor, linewidth=linewidth, label=label)
        if (points == True):
            ax.plot(np.arange(0,12,1)+0.5,data,"o",markerfacecolor=linecolor, 
             markeredgecolor="black", markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,data, linewidth=linewidth, label=label)
        if (points == True):
            ax.plot(np.arange(0,12,1)+0.5,data,"o",markeredgecolor="black", 
            markersize=10, markeredgewidth=2)

    return ax

def plotlinetime_j2j_monthly(fig, datatemp, x1, x2, y1, y2, titlestr, yrange=None, 
    yticks=None, yticklabels=None, ytitle=None, linecolor=None, points=True, label=None, linestyle='solid'):
    """ plot a line plot using monthly data from Jan to Dec
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

    data = datatemp.copy(deep=True)
    data[0:6] = np.array(datatemp[6:12])
    data[6:12] = np.array(datatemp[0:6])

    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])
 
    monticks = np.arange(0,12,1)
    monticks2 = np.arange(0,12,1)+0.5
    if (yrange):
        ax.set_ylim(yrange)

    if (yticks):
        ax.set_yticks(yticks)

    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=12)

    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)

    ax.set_xlim([0,12])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2, minor=True)
    ax.set_xticklabels(['J','A','S','O','N','D','J','F','M','A','M','J'], minor=True, fontsize=12)
    ax.set_title(titlestr, fontsize=16)

    if linecolor is not None:
            ax.plot(np.arange(0,12,1)+0.5,data,color=linecolor, linewidth=2, label=label, linestyle=linestyle)
            if (points == True):
                ax.plot(np.arange(0,12,1)+0.5,data,"o",markerfacecolor=linecolor, 
                markeredgecolor="black", markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,data, linewidth=2, label=label, linestyle=linestyle)
        if (points == True):
            ax.plot(np.arange(0,12,1)+0.5,data,"o",markeredgecolor="black", 
            markersize=10, markeredgewidth=2)

    return ax

def oplotlinetime_j2j_monthly(ax, datatemp, linecolor=None, linewidth=1, points=True, label=None):
    """ overplot a line on a January - December monthly line plot"""

    data = datatemp.copy(deep=True)
    data[0:6] = np.array(datatemp[6:12])
    data[6:12] = np.array(datatemp[0:6])

    if (linecolor[0]):
        ax.plot(np.arange(0,12,1)+0.5,data,color=linecolor, linewidth=linewidth, label=label)
        if (points == True):
            ax.plot(np.arange(0,12,1)+0.5,data,"o",markerfacecolor=linecolor, 
             markeredgecolor="black", markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,data, linewidth=linewidth, label=label)
        if (points == True):
            ax.plot(np.arange(0,12,1)+0.5,data,"o",markeredgecolor="black", 
            markersize=10, markeredgewidth=2)

    return ax






def plotlinetime_j2d_monthly_cdf(fig, data, x1, x2, y1, y2, titlestr, yrange=None, 
    yticks=None, yticklabels=None, ytitle=None, linecolor=None):
    """
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

    dataplot = np.zeros([12])
    for i in np.arange(0,12,1):
        dataplot[i] = np.sum(data[0:i+1])

    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])
 
    monticks = np.arange(0,12,1)
    monticks2 = np.arange(0,12,1)+0.5
    if (yrange):
        ax.set_ylim(yrange)

    if (yticks):
        ax.set_yticks(yticks)

    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=14)

    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)

    ax.set_xlim([0,12])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2, minor=True)
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], minor=True, fontsize=14)
    ax.set_title(titlestr, fontsize=16)

    if (linecolor[0]):
        ax.plot(np.arange(0,12,1)+0.5,dataplot,color=linecolor, linewidth=2)
        ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markerfacecolor=linecolor, markeredgecolor=linecolor, markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,dataplot, linewidth=2)
        ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markeredgecolor="black", markersize=10, markeredgewidth=2)

    return ax


def plotlinetime_j2d_monthly_accummean(fig, data, x1, x2, y1, y2, titlestr, yrange=None, 
    yticks=None, yticklabels=None, ytitle=None, linecolor=None, points=True):
    """
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

    dataplot = np.zeros([12])
    for i in np.arange(0,12,1):
        dataplot[i] = np.mean(data[0:i+1])

    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])
 
    monticks = np.arange(0,12,1)
    monticks2 = np.arange(0,12,1)+0.5
    if (yrange):
        ax.set_ylim(yrange)

    if (yticks):
        ax.set_yticks(yticks)

    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=14)

    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)

    ax.set_xlim([0,12])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2, minor=True)
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], minor=True, fontsize=14)
    ax.set_title(titlestr, fontsize=16)

    if (linecolor[0]):
        ax.plot(np.arange(0,12,1)+0.5,dataplot,color=linecolor, linewidth=2)
        if (points):
            ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markerfacecolor=linecolor, 
              markeredgecolor=linecolor, markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,dataplot, linewidth=2)
        if (points):
            ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markeredgecolor="black", markersize=10, 
             markeredgewidth=2)

    return ax




def oplotlinetime_j2d_cdf(ax, data, linecolor=None):
    """ overplot a line on a January - December monthly line plot"""

    dataplot = np.zeros([12])
    for i in np.arange(0,12,1):
        dataplot[i] = np.sum(data[0:i+1])

    if (linecolor):
        ax.plot(np.arange(0,12,1)+0.5,dataplot,color=linecolor, linewidth=2)
        ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markerfacecolor=linecolor, markeredgecolor=linecolor, markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,dataplot, linewidth=2)
        ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markeredgecolor="black", markersize=10, markeredgewidth=2)

    return ax

def oplotlinetime_j2d_accummean(ax, data, linecolor=None, points=True):
    """ overplot a line on a January - December monthly line plot"""

    dataplot = np.zeros([12])
    for i in np.arange(0,12,1):
        dataplot[i] = np.mean(data[0:i+1])

    if (linecolor[0]):
        ax.plot(np.arange(0,12,1)+0.5,dataplot,color=linecolor, linewidth=2)
        if (points):
            ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markerfacecolor=linecolor, 
             markeredgecolor=linecolor, markersize=10)
    else:
        ax.plot(np.arange(0,12,1)+0.5,dataplot, linewidth=2)
        if (points):
            ax.plot(np.arange(0,12,1)+0.5,dataplot,"o",markeredgecolor="black", 
            markersize=10, markeredgewidth=2)

    return ax




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

def oplotfill_j2j(ax, dat1, dat2, fillcolor='lightgray'):
    """ Over plot fill on a line plot already created using plotlinetime_j2j"""
    july1 = 181
    dat1plot = np.zeros([dat1.size])
    dat1plot[0:365-july1] = dat1[july1:365]
    dat1plot[365-july1:365] = dat1[0:july1]

    dat2plot = np.zeros([dat2.size])
    dat2plot[0:365-july1] = dat2[july1:365]
    dat2plot[365-july1:365] = dat2[0:july1]

    ax.fill_between(np.arange(0,365,1),dat1plot, dat2plot, color=fillcolor)
    return ax


def oplotrange_j2j(ax, minline, maxline, color=None):
    """overplot a range on a plot already created using plotlinetime_j2j"""
    july1 = 181
    minplot = np.zeros([minline.size])
    maxplot = np.zeros([maxline.size])
    minplot[0:365-july1] = minline[july1:365]
    minplot[365-july1:365] = minline[0:july1]
    maxplot[0:365-july1] = maxline[july1:365]
    maxplot[365-july1:365] = maxline[0:july1]

    if (color):
        ax.fill_between(np.arange(0,365,1), minplot, maxplot, color=color, alpha=0.5)
    else:
        ax.fill_between(np.arange(0,365,1), minplot, maxplot, color='salmon', alpha=0.5)

    return ax






 
