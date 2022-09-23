import xarray as xr
import numpy as np
from CASutils import linfit_utils as linfit
from CASutils import bootstrap_utils as boot
import sys

def dotheconstraint(xem, yem, x1mem, y1mem, obsx, sigxem=None, sigyem=None,sigx1mem=None,
                    sigy1mem=None, rxyem=None, rxy1mem=None, seed=None, nboots=1000, method='OLS', 
                    outputsamples=False):
    """ Performing the constraint 
    Inputs:
        xem = the ensemble mean predictor
        yem = the ensemble mean predictand
        x1mem = a single member predictor
        y1mem = a single member predictand
        obsx = the observed predictors
        sigxem = the standard deviation on the ensemble mean predictor
        sigyem = the standard deviation on the ensemble mean predictand
        sigx1mem = the standard deviation on the single member predictor
        sigy1mem = the standard deviation on the single member predictand
        rxyem = the ensemble mean correlation between predictor and predictand
        rxy1mem = the single member correlation between predictor and predictand
        seed = a random number seed which may be specified for reproducibility
        nboots = the number of bootstrap samples for each part. Default 1000.
        method = 'OLS', 'TLS' or 'BHM'
        outputsamples = if True, the samples used to build up the constrained distribution are output.
    Outputs: 
        datout = contains the mean, 95% and 66% confidence intervals and the fraction 
         of samples greater than the ensemble mean predictand for both the forced response plus
         internal variability and the forced response in isolation.
    """

    # check that all the needed things are there for the method
    if sigyem is None:
        print("You need to specify the sigma_y's for the ensemble mean for any constraint")
        sys.exit()
    if sigx1mem is None:
        print("You need to specify the sigma_x for 1 member for any constraint")
        sys.exit()
    if sigy1mem is None:
        print("You need to specify the sigma_y for 1 member for any constraint")
        sys.exit()

    if ((method == "TLS") or (method == "BHM")):
        if sigxem is None:
            print("You need to specify the sigma_x for the ensemble mean for TLS or BHM")
            sys.exit() 

    if ((method != "OLS") and (method !="TLS") and (method != "BHM")):
        print("You have chosen an invalid method.  Choose OLS, TLS or BHM")
        sys.exit()

    if (method == "BHM"):
        if rxyem is None:
            print("You need to specify rxy for the ensemble mean to use the BHM")
            sys.exit()
        if rxy1mem is None:
            print("You need to specify rxy for 1 member to use the BHM")
            sys.exit()

    # calculate the regression coefficients using the ensemble mean
    if (method == 'OLS'):
        print("Constraining using OLS")
        a, b = linfit.linfit_xy(xem, yem, sigma=sigyem)
    if (method == 'TLS'):
        print("Constraining using TLS") 
        a, b = linfit.tls(xem, yem, sigxem, sigyem)
    if (method == 'BHM'):
        print("Constraining using the BHM")
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3)

    if (method == "BHM"): 
        # an array of standard deviations for the forced noise
        sigforced = np.sqrt(del2[:])
        # standard deviation for the internal variability noise component
        sigyiv=np.sqrt( (sigy1mem**2)*(1.-rxy1mem**2.))
    else:
        # calculate the single member residuals from the linear regression fit
        # their standard deviation and the standard deviation of the noise term
        # both with (sigwithiv) and without (sigforced) internal variability
        eps = y1mem[:] - (a + b*x1mem[:])
        sigeps = np.std(eps)
        sigwithiv = np.sqrt(sigeps**2 - (b**2)*sigx1mem)
        sigyiv = sigy1mem
        sigforced = np.sqrt(sigwithiv**2. - sigyiv**2.)

    # sampling the uncertainty on the observed predictor
    # 250 values for each observational value
    nobs=obsx.size
    obspdf = np.zeros([nobs*250])
    obstrue = np.zeros([nobs*250])
    obssample = np.random.normal(0,sigx1mem,250)
    for iobs in range(0,obsx.size,1):
        obspdf[iobs*250:(iobs+1)*250] = obsx[iobs] + obssample[:]
        if (method == "BHM"):
            obstrue[iobs*250:(iobs+1)*250] = obsx[iobs]


    # combine all the sampling
    # OLS and TLS
    if (method != "BHM"):

        # sample the noise terms and regression coefficients
        randomvals = np.random.normal(0,1,nboots)
        # forced + internal
        noise_withiv = randomvals*np.array(sigwithiv)
        # forced
        noise_forced = randomvals*np.array(sigforced)
        if (method == "OLS"):
            sigxin=None
        if (method == "TLS"):
            sigxin=sigxem
        aboots, bboots = boot.boot_regcoefs(xem, yem, sigx=sigxin, sigy=sigyem)

        # first, regression coefficient uncertainty with observational predictor uncertainty
        y = np.zeros([nobs*250*nboots])
        for iboot in range(0,nboots,1):
            y[iboot*nobs*250:(iboot+1)*nobs*250] = aboots[iboot] + bboots[iboot]*obspdf[:]

        # now adding on the noise terms
        yplusiv = np.zeros([nobs*250*nboots*nboots])
        yforced = np.zeros([nobs*250*nboots*nboots])
        for iboot in range(0,nboots,1):
            yplusiv[iboot*(nobs*250*nboots):(iboot+1)*(nobs*250*nboots)] = y[:] + noise_withiv[iboot]
            yforced[iboot*(nobs*250*nboots):(iboot+1)*(nobs*250*nboots)] = y[:] + noise_forced[iboot]
 
    else:
        #sample the noise terms
        # only do internal variability here because forced noise is 
        # dependend on delxdelx which is paired with alpha and beta's
        randomvals = np.random.normal(0,1,nboots*nboots)
        randomvals2 = np.random.normal(0,1,nboots*nboots)
        noiseiv = np.array(sigyiv)*randomvals2[:]
 
        yplusiv = np.zeros([nobs*250*nboots*nboots])
        yforced = np.zeros([nobs*250*nboots*nboots])
        for iboot in range(0,nboots,1):
            y=np.zeros([nobs*250])
            y[:] = aboots[iboot] + bboots[iboot]*obspdf[:]
            noise_forced = randomvals[:]*sigforced[iboot]
            
            yplusivt = np.zeros([nobs*250*nboots])
            yforcedt = np.zeros([nobs*250*nboots])
            for inoise in range(0,nboots,1):
                yforcedt[inoise*nboots:(inoise+1)*nboots]=\
                y[inoise]+noise_forced[inoise*nboots:(inoise+1)*nboots]
        
                yplusivt[inoise*nboots:(inoise+1)*nboots]=\
                yforcedt[inoise*nboots:(inoise+1)*nboots]+\
                noiseiv[inoise] + \
                np.array(rxy1mem)*(np.array(sigy1mem)/np.array(sigx1mem))*(obspdf[:]-obstrue[:])
               
            yforced[iboot*nobs*250*nboots:(iboot+1)*nobs*250*nboots]=yforcedt[:]
            yplusiv[iboot*nobs*250*nboots:(iboot+1)*nobs*250*nboots]=yplusivt[:] 

    meanwithiv = np.mean(yplusiv)
    meanforced = np.mean(yforced)

    print("starting to calculate the percentiles - this could take a while")
    min95withiv = np.percentile(yplusiv, 2.5)
    max95withiv = np.percentile(yplusiv, 97.5)
    min66withiv = np.percentile(yplusiv, 17)
    max66withiv = np.percentile(yplusiv, 83)

    min95forced = np.percentile(yforced, 2.5)
    max95forced = np.percentile(yforced, 97.5)
    min66forced = np.percentile(yforced, 17)
    max66forced = np.percentile(yforced, 83)

    print("calculating the percentage greater than the ensemble mean - this may also take a while")
    ymean = np.mean(yem)
    numberabovey = np.sum(yplusiv > np.array(ymean))
    gtymean_withiv = (numberabovey/float(nobs*250*nboots*nboots))*100.
    numberabovey = np.sum(yforced > np.array(ymean))
    gtymean_forced = (numberabovey/float(nobs*250*nboots*nboots))*100.


    datout={ "meanwithiv":meanwithiv, "meanforced":meanforced,
             "min95withiv":min95withiv, "max95withiv":max95withiv,
             "min66withiv":min66withiv, "max66withiv":max66withiv,
             "min95forced":min95forced, "max95forced":max95forced,
             "min66forced":min66forced, "max66forced":max66forced,
             "gtymean_withiv":gtymean_withiv, "gtymean_forced":gtymean_forced}

    if (outputsamples):
        return datout, yplusiv
    else: 
        return datout




