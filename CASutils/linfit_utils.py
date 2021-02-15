import xarray as xr
import numpy as np
import scipy.optimize as opt 
import sys
import statsmodels.api as sm

# various linear regression functions
# linfit_lonlat = calculating the trend for a lat,lon array
# linfit_xy = using ordinary least squares fit
# multi_ols = multiple linear regression using OLS
# tls = total least squares regresion
# bhm = using a Bayesian Hierarchical Model.

def linfit_lonlat(darray,dimname):
    """Calculate the linear trend coefficients for an [:, lat, lon] array"""
    darray_2d=darray.stack(allpoints=['lat','lon'])
    result=np.polyfit(darray_2d[dimname], darray_2d, 1)
    b=result[0,:] ; a=result[1,:]
    a_xr=xr.DataArray(a, coords=[darray_2d.allpoints], name='a')
    b_xr=xr.DataArray(b, coords=[darray_2d.allpoints], name='b')
    a_xr=a_xr.unstack('allpoints')
    b_xr=b_xr.unstack('allpoints')

    return a_xr, b_xr

def linfit_xy(x,y, sigma=None):
    """Calculate a weighted least squares linear fit between 1d arrays x and y
    Input: x (the predictor)
           y (the predictand)
           sigma (optional) the standard deviation of the y values
    Output: a and b and the residuals
    """

    if sigma is not None:
        w = 1./(sigma)
        coefs = np.polyfit(x,y,1,w=w)
    else:
        coefs = np.polyfit(x,y,1)

    a = coefs[1]
    b = coefs[0]
    return a, b

def multi_ols(x, y, sigma=None):
    """ Calculate a weighted least squares regression with multiple X variables
    Input: x (numpy array of size nobs x npredictors)
           y (the predictand)
           sigma (optional) the standard deviation of the y values
    Outbut: a and bs for each predictor
    """
    x = np.array(x)
    x = sm.add_constant(x)
    
    y = np.array(y)
    
    if sigma is not None: 
        sigma = np.array(sigma)
        w = 1./(sigma**2)
        model = sm.WLS(y, x, weights = w)
    else:
        model = sm.WLS(y, x)
        
    result = model.fit()
    return result.params

def tls(xin ,yin, sigmaxin, sigmayin):
    """Calculate the total least squares regression with errors in both x and y.
       by minimizing the chi-sq ( y - a -bx)/(sigmay^2 + b^2*sigmax^2)
    """
    x = np.array(xin)
    y = np.array(yin)
    sigmax = np.array(sigmaxin)
    sigmay = np.array(sigmayin)

    # normalize
    xm = np.mean(x) ; ym =np.mean(y)
    xsd = np.std(x) ; ysd = np.std(y)
    x = (x[:] - xm)/xsd
    y = (y[:] - ym)/ysd
    sigmay = sigmay[:]/ysd
    sigmax = sigmax[:]/xsd

    # compute first guess for B slope using standard 1-D least squares fit
    # where the errors in x are ignored 
    a, b = linfit_xy(x,y,sigma=sigmay)
    firstfit = (b,a)

    # Find a and b that minimize chisq 
    def chisq_func(coefs):
        denom = (sigmay[:]**2 + (coefs[1]**2)*(sigmax[:]**2))
        num = ( y[:] - (coefs[0] + coefs[1]*x[:]))**2
        chisq = np.sum(num/denom)
        return chisq

    coefs = np.array([a, b])
    coefsin=coefs
    result = opt.minimize(chisq_func, coefs)

    coefs=result.x
    a = coefs[0]
    b = coefs[1]

    # un-normalize
    bout = b*(ysd/xsd)
    aout = ym + a*ysd - b*(ysd/xsd)*xm 
    return aout, bout

def bhm(xin, yin, sigxin, sigyin, rxyin, ntrue=1000, nburn = 30, nparams=30, iseed=None):
    """
    Bayesian linear regression using Gibbs sampling.
    Assumes correlated errors in x and y
    -------
    Input:
    xin = x observations 
    yin = y observations
    sigx = sigma for the x values
    sigy = sigma for the y values
    rxy = the correlation between x and y errors
    Optional input:
    ntrue = the number of samples for xt and xt (true x and y values)
    nburn = the number of iterations to throw out at the beginning
    nparams = the number of iterations for sampling parameters
    iseed = a specified random number seed (for reproducibility)
    Output:
    aout = ntrue values of the a parameter of y_t=a+bx_t 
    bout = ntrue values of the b parameter of y_t=a+bx_t
    deldelout = ntrue values of sigma^2_delta 
    muxout = the mean of the true x values
    delxdelxout = the variance of the true x values
    """

    # make sure input's are numpy arrays
    x = np.array(xin)
    y = np.array(yin)
    sigx = np.array(sigxin)
    sigy = np.array(sigyin)
    rxy = np.array(rxyin)

    # normalize inputs
    xm = np.mean(x) ; ym =np.mean(y)
    xsd = np.std(x) ; ysd = np.std(y)
    x = (x[:] - xm)/xsd
    y = (y[:] - ym)/ysd
    sigy = sigy[:]/ysd
    sigx = sigx[:]/xsd
    
    # sigma squareds
    sigx2 = sigx**2.
    sigy2 = sigy**2.

    nobs = x.size

    if (x.size != y.size):
        print("x and y don't have the same dimensions")
        sys.exit("blr")

    # set seed if not specified
    if (iseed):
        np.random.seed(iseed)
    
    # set seed for each sampler
    #iseed_xt = iseed+1
    #iseed_yt = iseed+2
    #iseed_paramorder = iseed+3
    #iseed_a = iseed+4
    #iseed_b = iseed+5
    #iseed_sigsig = iseed+6
    #iseed_mux = iseed+8
    #iseed_delxdelx = iseed+9

    # initialize output arrays
    av = np.zeros(ntrue)
    bv = np.zeros(ntrue)
    deldel =np.zeros(ntrue)
    mux = np.zeros(ntrue)
    delxdelx = np.zeros(ntrue)
   # varyt = np.zeros(ntrue) # variance of y_true

    # initialize unknown parameters with first guesses
    
    xt = x ; yt = y # set true vals to observed vals
    iav, ibv = linfit_xy(x,y, sigma=sigy) # ols params for a amnd b
    ydetrend = y[:] - (iav + ibv*x[:])
    ideldel = np.var(ydetrend)
    imux = np.mean(x)
    # initialize idelxdelx with the excess variance beyond the variance
    # attached to the x observations.  Ensure there is more variance in
    # the x values than the uncertainty in each individual x value. 
    idelxdelx = np.var(x) - np.mean(sigx2)
    if (idelxdelx < 0):
        print("Variance in x values is less than the mean error variance on x values")
        sys.exit("x variance issue")

    # set ceilings to prevent inverse gamma functions from blowing up
    deldelceiling = 10.
    delxdelxceiling = 10.

    # iterate through xt, yt and all 5 unknown parameters
    for i in range(0,ntrue+nburn,1):
    #for i in range(0,1,1):
       
        # sampling y_t
        mean_y = np.zeros(y.size)
        var_y = np.zeros(y.size)
        var_y[:] = ( (1./ideldel) + (1./(sigy2[:]*(1. - rxy[:]**2.))) )**(-1.)
        mean_y[:] = ( ((ibv*xt[:] + iav)/ideldel) +
                    ( (y[:] - (rxy[:]*(np.sqrt(sigy2[:]/np.sqrt(sigx2[:]))*(x[:]-xt[:]))))/
                    (sigy2[:]*(1-rxy[:]**2))))*var_y[:]
        yt = np.sqrt(var_y[:])*np.random.normal(0,1,nobs) + mean_y[:]
 
        # sampling x_t
        mean_x = np.zeros(x.size)
        var_x = np.zeros(x.size)
        var_x[:] = ( (ibv**2./ideldel) + (1./(sigx2[:]*(1-rxy[:]**2.))) + (1./idelxdelx))**(-1.) 
        mean_x[:] = ( (ibv*(yt[:]-iav)/ideldel) + 
                      ((x[:] -rxy[:]*(np.sqrt(sigx2[:])/np.sqrt(sigy2[:]))*(y[:]-yt[:]))/
                      (sigx2*(1.-rxy[:]**2.)) ) + 
                      (imux/idelxdelx))*var_x[:]
        xt = np.sqrt(var_x[:])*np.random.normal(0,1,nobs) + mean_x[:]

        # iterative loop for sampling the 5 unknowns
        for j in range(0,nparams,1):
            
            # sampling a
            amean = (np.sum(y[:]) - ibv*np.sum(x[:]))/float(nobs)
            astdev = np.sqrt( ideldel /float(nobs) )
            iav = astdev * np.random.normal(0,1,1) + amean

            # sampling b
            bmean = (np.sum(xt[:]*yt[:]) - np.sum(xt[:]*iav))/np.sum(xt[:]**2.)
            bstdev = np.sqrt( ideldel / np.sum(xt[:]**2.) )
            ibv = bstdev * np.random.normal(0,1,1) + bmean

            # sampling delta^2
            igalpha = float(nobs)/2.
            igbeta = np.sum( (yt[:] - (iav + ibv*xt[:]))**2.) / 2.
            gammaval = np.random.gamma( igalpha, 1./igbeta, 1 )
            ideldel = (1./gammaval)
            if (ideldel >  deldelceiling):
                 ideldel = deldelceiling

            # sampling mu_x
            mux_mean = np.sum(xt[:])/float(nobs)
            mux_stdev = np.sqrt( (idelxdelx)/float(nobs))
            imux = mux_stdev * np.random.normal(0,1,1) + mux_mean

            # sampling delx^2
            igalpha = float(nobs)/2.
            igbeta = np.sum( (xt[:] - imux)**2.)/2.
            gammaval = np.random.gamma( igalpha, 1./igbeta, 1)
            idelxdelx = (1./gammaval)
            if (idelxdelx > delxdelxceiling): 
                idelxdelx = delxdelxceiling


        #print(str(np.var(yt))+'='+str(idelxdelx*(ibv**2) + ideldel))

        if (i >= nburn):
            av[i-nburn] = iav
            bv[i-nburn] = ibv
            deldel[i-nburn] = ideldel
            mux[i-nburn] = imux
            delxdelx[i-nburn] = idelxdelx
            ytrenorm = yt[:]*ysd
           # varyt[i-nburn] = np.var(ytrenorm)
    
    #re-normalize
    aout = ym + av[:]*ysd - bv[:]*(ysd/xsd)*xm
    bout = bv[:]*(ysd/xsd)
    delxdelxout = delxdelx[:]*(xsd**2.)
    muxout = xsd*mux[:] + xm
    deldelout = (ysd**2.)*deldel[:]

    return aout, bout, deldelout, muxout, delxdelxout



def compute_slope(y):
    """
    function to compute slopes along a dimension of an xarray DataArray
    combine to apply_ufunc to implement
    
    Example usage: 
    xr.apply_ufunc(linfit.compute_slope, da, vectorize=True, input_core_dims=[['time']])
    """
    x = np.arange(len(y))
    return np.polyfit(x, y, 1)[0]





































