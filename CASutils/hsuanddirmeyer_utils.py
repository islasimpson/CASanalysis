import numpy as np
from scipy import optimize
import xarray as xr
from math import nan

"""
Model functions for linear, 2 segment piecewise linear, and 3 segment piecewise linear

x = independent variable
x0 = x value at critical point 0
x1 = x value at critical point 1
y0 = y value at critical point 0 (or y intercept)
k1 = slope of segment 1
k2 = slope of segment 2
k3 = slope of segment 3

Lilian Zhu, with additions from Isla Simpson (2025)

"""

#--- Wrapper to call the fitting functions and output the data with additional metadata
def fit_hd2022(x, y, curve_type=None):
    curve_type, p, bic = curve_fit(np.array(x), np.array(y), curve_type=curve_type)
    
    curve_type_out = xr.DataArray(curve_type, name='curve_type')
    bic_out = xr.DataArray(bic, name='bic') 

    if (curve_type == '111'):
        names=['x0','x1','y0','k1','k2','k3']
        long_names=['x value at critical point 0',
                    'x value at critical point 1',
                    'y value at critical point 0 (or y intercept)',
                    'slope of segment 1',
                    'slope of segment 2',
                    'slope of segment 3']
        p_xr = xr.Dataset({
           names[i]: xr.DataArray(p[i], attrs={"long_name":long_names[i]}) 
           for i in range(len(p))
           })

    if ((curve_type == '010') | (curve_type == '001')):
        names=['y0','k']
        long_names=['y intercept', 'slope']
        p_xr = xr.Dataset({
           names[i]: xr.DataArray(p[i], attrs={"long_name":long_names[i]})
           for i in range(len(p))
           })

    if ((curve_type == '110') | (curve_type == '011')):
        names=['x0','y0','k1','k2']
        long_names=['x value at critical point',
                    'y value at critical point',
                    'slope of segment 1',
                    'slope of segment 2']
        p_xr = xr.Dataset({
           names[i]: xr.DataArray(p[i], attrs={"long_name": long_names[i]})
           for i in range(len(p))
           })

    fit_output = xr.merge([curve_type_out, bic_out, p_xr])
    return fit_output

#--- Function for overplotting the fit over some range
def oplot_fit(ax, xmin, xmax, fitdat, color='red', linewidth=2, linestyle='solid'):
    """ Over plotting the fit over the range xmin to xmax on axis ax using
        fitdata which is the output to the above fit_hd2022 function
    """

    if (fitdat.curve_type == '111'):

        x0 = np.array(fitdat.x0)
        x1 = np.array(fitdat.x1)
        y0 = np.array(fitdat.y0)
        k1 = np.array(fitdat.k1)
        k2 = np.array(fitdat.k2)
        k3 = np.array(fitdat.k3)
        y1 = y0 + k2*(x1 - x0)

        #- first segment
        xvals = np.array([xmin, x0])
        ax.plot(xvals, k1*xvals + y0, color=color, linewidth=linewidth, linestyle=linestyle)
        #- second segment
        xvals = np.array([x0, x1])
        ax.plot(xvals, y0 + k2*(xvals - x0), color=color, linewidth=linewidth, linestyle=linestyle)
        #- third segment
        xvals = np.array([x1, xmax])
        ax.plot(xvals, k3*(xvals - x1) + y1, color=color, linewidth=linewidth, linestyle=linestyle)

    if ((fitdat.curve_type == '010') | (fitdat.curve_type == '001')):
        y0 = np.array(fitdat.y0)
        k = np.array(fitdat.k)
        xvals = np.array([xmin, xmax])
        ax.plot(xvals, k*xvals + y0, color=color, linewidth=linewidth, linestyle=linestyle)

    if ((fitdat.curve_type == '110') | (fitdat.curve_type == '011')):
        y0 = np.array(fitdat.y0)
        x0 = np.array(fitdat.x0)
        k1 = np.array(fitdat.k1)
        k2 = np.array(fitdat.k2)

        #- first segment
        xvals = np.array([xmin, x0])
        ax.plot(xvals, k1*(xvals-x0) + y0, color=color, linewidth=linewidth, linestyle=linestyle)
        xvals = np.array([x0, xmax])
        ax.plot(xvals, y0 + k2*(xvals - x0), color=color, linewidth=linewidth, linestyle=linestyle)


    return ax

def get_percentage_in_transitional(x, fitdat):
    """ compute the percentage of data points within the transitional regime
        based on the input x values (x) and the Hsu and Dirmeyer fit parameters (fitdat) """
    if (fitdat.curve_type == '111'):
        x0 = np.array(fitdat.x0)
        x1 = np.array(fitdat.x1)
        transx = x.where( (x >= x0) & (x <= x1), drop=True)
        transper = (len(transx) / len(x))*100.
    if ( (fitdat.curve_type == '001') | (fitdat.curve_type == '100') ):
        transper = 0
    if (fitdat.curve_type == '010'):
        transper = 100
    if (fitdat.curve_type == '011'):
        x0 = np.array(fitdat.x0)
        transx = x.where( (x <= x0), drop=True)
        transper = (len(transx) / len(x))*100.
    if (fitdat.curve_type == '110'):
        x0 = np.array(fitdat.x0)
        transx = x.where( (x >= x0), drop=True)
        transper = (len(transx) / len(x))*100.
    return transper

def get_critical_point(fitdat):
    """ find the critical point """

    if (fitdat.curve_type == '111'):
        cp = fitdat.x1.item()
    if (fitdat.curve_type == '011'):
        cp = fitdat.x0.item()
    if ((fitdat.curve_type != '111') & (fitdat.curve_type != '011')):
        cp = nan

    return cp



def get_transitional_slope(fitdat):
    """ Output the slope of the transitional regime based on the Hsu and Dirmeyer fit parameters (fitdat) """
    if (fitdat.curve_type == '111'):
        b = fitdat.k2.values.item()    
    if (fitdat.curve_type == '010'):
        b = fitdat.k.values.item()
    if (fitdat.curve_type == '001'):
        b = nan
    if ( (fitdat.curve_type == '011') | (fitdat.curve_type == '110') ):
        b = np.max([fitdat.k1, fitdat.k2])
    return b


#--- linear model
def single_linear(x,y0,k1):
    return y0+k1*x

#--- from Hsu and Dirmeyer 2 segment
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

#--- from Hsu and Dirmeyer 3 segment (modified by removing y1)
def piecewise3sg_linear(x, x0, x1, y0, k1, k2, k3):
    condlist = [x < x0, (x >= x0) & (x < x1), x >= x1]
    funclist = [lambda x:k1*(x-x0) + y0, lambda x:k2*(x-x0) + y0, lambda x:k2*(x1-x0) + y0 + k3*(x-x1)]
    return np.piecewise(x, condlist, funclist)

"""
Curve fitting functions for 001, 010, 110, 011, and 111
"""
def curve_001(x,y):
    """
    order: y0, k1

    initial guesses (p0):
        y0 = median of y
        k1 = 0

    bounds:
        min(y) < y0 < max(y)
        0 < k1 < 0.01
    """
    p0 = [np.median(y),0]  # y0 and k1=0
    bounds = ([np.min(y),0],  # bounds y0 to range of y, 0<k1<0.01
              [np.max(y),0.01])

    p,e = optimize.curve_fit(single_linear,x,y,p0=p0,bounds=bounds)
    yfit = single_linear(x,*p)

    return p,yfit

def curve_010(x,y):
    """
    order: y0, k1

    initial guesses (p0):
        y0 = median of y
        k1 = 50

    bounds:
        min(y) < y0 < max(y)
        0.05 < k1 < 1000
    """
    p0 = [np.median(y),50]  # y0 and k1 = 50
    bounds = ([np.min(y),0.01],  # bounds y0 to range of y, 0.05<k1<1000
              [np.max(y),1000])

    p,e = optimize.curve_fit(single_linear,x,y,p0=p0,bounds=bounds)
    yfit = single_linear(x,*p)

    return p,yfit

def curve_110(x,y):
    """
    order: x0,y0,k1,k2

    initial guesses (p0):
        x0 -> midrange of x
        y0 -> median of y
        k1 = 0, k2 = 50

    bounds:
        min(x) < x0 < max(x)
        min(y) < y0 < max(y)
        -0.0001 < k1 < 0.0001  (near zero slope)
        0.05 < k2 < 1000 (positive slope)
    """
    p0 = [(np.max(x)+np.min(x))/2,np.median(y),0,50]  # x0, y0, k1=0, k2=50

    bounds = ([np.min(x),np.min(y),-0.01,0.05],  # bounds x0 and y0 to range of data
              [np.max(x),np.max(y),0.01,1000])  # -0.001<k1<0.001 and 0.001<k2<1000

    p,e = optimize.curve_fit(piecewise_linear,x,y,p0=p0,bounds=bounds)
    yfit = piecewise_linear(x, *p)

    return p,yfit

def curve_011(x,y):
    """
    order: x0,y0,k1,k2

    initial guesses (p0):
        x0 -> midrange of x
        y0 -> median of y
        k1 = 50, k2 = 0

    bounds:
        min(x) < x0 < max(x)
        min(y) < y0 < max(y)
        0.05 < k1 < 1000 (positive slope)
        -0.0001 < k2 < 0.0001  (near zero slope)
    """
    p0 = [(np.max(x)+np.min(x))/2,np.median(y),50,0]  # x0, y0, k1=1, k2=0

#    bounds = ([np.min(x),np.min(y),0.05,-0.01],  # bounds x0 and y0 to range of data
#              [np.max(x),np.max(y),1000,0.01])  # 0<k1<1000 and -0.001<k2<0.001
    bounds = ([np.min(x),np.min(y),0.01,-0.01],  # bounds x0 and y0 to range of data
              [np.max(x),np.max(y),1000,0.01])  # 0<k1<1000 and -0.001<k2<0.001




    p,e = optimize.curve_fit(piecewise_linear,x,y,p0=p0,bounds=bounds)
    yfit = piecewise_linear(x, *p)

    return p,yfit

def curve_111(x,y):
    """
    order: x0,x1,y0,k1,k2,k3

    initial guesses (p0):
        x0 and x1 -> midrange of x
        y0 -> median of y
        k1 = 0, k2 = 0.1, k3 = 0

    bounds:
        min(x) < x0,x1 < max(x)
        min(y) < y0 < max(y)
        -0.001 < k1,k3 < 0.001  (near zero slope)
        0.05 < k2 < 1000 (positive slope)
    """
    p0 = [(np.max(np.array(x)+0.001)+np.min(np.array(x)))/2,(np.max(np.array(x)+0.001)+np.min(np.array(x)))/2,
    np.median(np.array(y)),0,0.1,0]

#    bounds = ([np.min(np.array(x)),np.min(np.array(x)),np.min(np.array(y)),-0.001,0.05,-0.001],
#              [np.max(np.array(x)+0.001),np.max(np.array(x)+0.001),np.max(np.array(y))+0.001,0.001,1000.0,0.001])

    bounds = ([np.min(np.array(x)),np.min(np.array(x)),np.min(np.array(y)),-0.01,0.01,-0.01],
              [np.max(np.array(x)+0.001),np.max(np.array(x)+0.001),np.max(np.array(y))+0.001,0.01,1000.0,0.01])

    # curve fitting
    p,e = optimize.curve_fit(piecewise3sg_linear,x,y,p0=p0,bounds=bounds)
    yfit = piecewise3sg_linear(x, *p)

    return p,yfit

#--- compute BIC
def compute_BIC(candidate,x,y,yfit):
    if candidate in ['001','010']:
        k = 2
    elif candidate in ['110','011']:
        k = 4
    elif candidate == '111':
        k = 6
    else:
        raise ValueError(f"Unknown candidate code '{candidate}'. "
                         "Expected one of {'001', '010', '110', '011', '111'}.")

    n = len(x)
    RSS = np.sum(np.square(y-yfit))
    BIC = n*np.log(RSS/n) + k*np.log(n)

    return BIC


"""
Fits all 5 models to input data and selects lowest BIC Model

params:
    x: array-like
        Independent variable values
    y: array-like
        Depedent variable values

Returns tuple of (p,yfit) values where:
    p is a tuple of parameters corresponding to model with lowest BIC
    yfit is an array of predicted y-values from best fitting model
"""
def curve_fit(x,y,curve_type=None):
    if curve_type is not None:
        if curve_type == '001':
            p,yfit = curve_001(x,y)
        elif curve_type == '010':
            p,yfit = curve_010(x,y)
        elif curve_type == '110':
            p,yfit = curve_110(x,y)
        elif curve_type == '011':
            p,yfit = curve_011(x,y)
        elif curve_type == '111':
            p,yfit = curve_111(x,y)

        bic = compute_BIC(curve_type,x,y,yfit)
        curve_type = xr.DataArray(curve_type, name='curve_type')
        return curve_type, p, bic

    else:
        candidates = ['001', '010', '110', '011', '111']
        BIC = []; model_fits = []
        for i,candi in enumerate(candidates):
            if candi == '001':
                p,yfit = curve_001(x,y)
            elif candi == '010':
                p,yfit = curve_010(x,y)
            elif candi == '110':
                p,yfit = curve_110(x,y)
            elif candi == '011':
                p,yfit = curve_011(x,y)
            elif candi == '111':
                p,yfit = curve_111(x,y)

            BICi = compute_BIC(candi,x,y,yfit)
            BIC.append(BICi)
            model_fits.append((candi,p,BICi))

        print(BIC)
        min_BIC_idx = np.argmin(BIC)

        # Now check whether there's a simpler model that has a BIC less than 10 
        # more than the chosen model
        # if min_BIC_idx == 0 or 1, do nothing
        # if min_BIC_idx == 4 chck if 2 or 3 have less than 10 more than 4 and choose the minimum 
        if (min_BIC_idx == 4):
            dif2 = BIC[2]-BIC[4]
            dif3 = BIC[3]-BIC[4]
            if ( (dif2 < 10) | (dif3 < 10) ):
                min_BIC_idx = np.argmin([9999,9999,BIC[2],BIC[3]])
        # if min_BIC_idx == 2 or 3 then check if the BIC for 0 or 1 is less than 10 more than BIC of 2 or 3
        # If so, choose whatever is the minimum of 0 or 1
        # Note that if the original BIC was 4 and 2 or 3 were less than 10 away from it, and 1 and 2 are then less than 10 away from 2 or 3
        # then a BIC that's more than 10 (but less than 20) away from the original could be chosen.
        if (min_BIC_idx == 2):
            dif0 = BIC[0] - BIC[2]
            dif1 = BIC[1] - BIC[2]
            if ( (dif0 < 10) | (dif1 < 10) ):
                min_BIC_idx = np.argmin([BIC[0],BIC[1]])
        if (min_BIC_idx == 3):
            dif0 = BIC[0] - BIC[3]
            dif1 = BIC[0] - BIC[3]
            if ( (dif0 < 10) | (dif1 < 10) ):
                min_BIC_idx = np.argmin([BIC[0],BIC[1]]) 
        


        return model_fits[min_BIC_idx]


#----functions for generating synthetic data
def gen_001_or_010(y0,k1,xmin=0,xmax=4,n=100,noise_std=0.03):
    x = np.random.uniform(xmin,xmax,n)
    y = y0 + k1*x + np.random.normal(0.0, noise_std, n)
    return x,y

def gen_110_or_011(x0,y0,k1,k2,xmin=0,xmax=4,n=100,noise_std=0.03):
    x = np.random.uniform(xmin,xmax,n)
    y = np.where(x < x0,
                 y0 + k1 * (x-x0),
                 y0 + k2 * (x-x0)
                ) + np.random.normal(0.0, noise_std, n)
    return x,y

def gen_111(x0,y0,x1,y1,k1,k2,k3,xmin=0,xmax=4,n=100,noise_std=0.03):
    x = np.random.uniform(xmin,xmax,n)
    y = np.where(
        x < x0,
        y0 + k1 * (x - x0),
        np.where(
            x < x1,
            y0 + k2 * (x - x0),
            y0 + k2 * (x1 - x0) + k3 * (x - x1)
        )
    ) + np.random.normal(0.0, noise_std, n)
    return x,y




