import xarray as xr
import numpy as np
import sys
from math import nan
from CASutils import linfit_utils as linfit

def bootdif2means(dat1, dat2):
    """Obtain the significance of the difference between two means using bootstrapping """
    dat1boot = bootgen(dat1)
    dat2boot = bootgen(dat2)

    dat1bootm = dat1boot.mean('isample')
    dat2bootm = dat2boot.mean('isample')
 
    diff = dat2bootm - dat1bootm
    diffmin = diff.quantile(0.025, dim='iboot')
    diffmax = diff.quantile(0.975, dim='iboot')

    dimsignif = diffmin.dims
    coordsignif = diffmin.coords

    signif = diffmin*0+1
    signif = signif.where( (diffmin < 0) & (diffmax > 0), nan)

    return signif


def bootgen(darray, nsamples=None, seed=None, nboots=1000):
    """Generate nboots bootstrap samples from darray with nsamples within each bootstrap.
    Sampling is done from the left most dimension
    If nsamples = None then nsamples = the length of the left most dimension.

    Input:
        darray = the data array on which you want to do the resampling (on the left most dimension)

    Optional input:
        nsamples = the number of members to go into each bootstrap samples.
        seed = an optional number to put in for the random number seed.  Required in 
               cases where reproducibility is needed.
        nboots = the number of bootstrap samples consisting of nsamples each

    """

    ### exit if darray is a dataset.
    if (str(type(darray)) == "xarray.core.dataset.Dataset"):
        print("this function doesn't accept datasets, convert to data array")
        sys.exit()

    ###if it's an xarray dataset, set up the dimensions
    try:
        dims = darray.dims
        if nsamples is None:
            nsamples = darray[dims[0]].size
 
        nmemin = darray[dims[0]].size

        dimboot = [nsamples*nboots]
        dimboot2d = [nboots, nsamples]
        bootcoords = [('iboot', np.arange(0,nboots,1)), ('isample', np.arange(0,nsamples,1))]
        for icoord in range(1,len(dims)):
            dimboot.append(darray[dims[icoord]].size)
            dimboot2d.append(darray[dims[icoord]].size)
           # bootcoords.append( (dims[icoord], darray[dims[icoord]] ))
            bootcoords.append( (dims[icoord], np.array(darray[dims[icoord]])) )
        print("you are using an xarray dataarray")
   
    except:
        if nsamples is None:
            nsamples = darray.shape[0]
        nmemin = darray.shape[0]

        dimboot = [nsamples*nboots]
        dimboot2d = [nboots, nsamples]
        for icoord in range(1,len(darray.shape)):
            dimboot.append(darray.shape[icoord])
            dimboot2d.append(darray.shape[icoord])

        print("you are using a numpy array")

    ### generate random number for bootstrapping
    if (seed):
        np.random.seed(seed)

    ### do the resampling
    ranu = np.random.uniform(0,nmemin,nboots*nsamples)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.array(darray[ranu])
    bootdat = bootdat.reshape(dimboot2d)

    #print(bootdat)

    try:
        bootdat = xr.DataArray(bootdat, coords=bootcoords)
    except:
        pass

    return bootdat

def bootgen_multimem(darray, nmems, seed=None, nboots=1000):
    """ Generate nboots bootstrap samples from darray with nmems members for each sample
        calculates the mean over members 

    Input: darray = an xarray data array with the sampling being performed on the first dimension
           nmems = the number of members in each bootstrap sample
           nboots = the number of bootstrap samples (optional)

    Output: bootdatxr = an xarray data array containing the bootstrap samples
            with dimensions (nboots, all the other coords of darray except the first)

    Option: a seed if you want to specify the seed for the random number generator 
    """

    # set up the dimensions and coordinates of the bootstrap array
    dims = darray.dims
    dimboot = [nmems*nboots]
    dimboot2d = [nmems, nboots]
    bootcoords = [("iboot", np.arange(0,nboots,1))]

    for icoord in range(1,len(dims)):
        dimboot.append(darray[dims[icoord]].size)
        dimboot2d.append(darray[dims[icoord]].size)
        bootcoords.append( (dims[icoord], darray[dims[icoord]] ))

    # generate random numbers for bootstrapping
    if (seed):
        np.random.seed(seed)

    nmemin = darray[dims[0]].size
    ranu = np.random.uniform(0,nmemin,nboots*nmems)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros(dimboot)
    bootdat = np.array(darray[ranu])
    bootdat = bootdat.reshape(dimboot2d)
    bootdatm = np.mean(bootdat, axis=0)
    bootdatxr = xr.DataArray(bootdatm, coords=bootcoords)


    return bootdatxr



def bootgenchunk_multimem(darray, nyears, nmems, nboots, seed=None):
    """Generate nboot samples with nmems members containing chunks of length nyears"""

    ### exit if darray is a dataset.
    if (str(type(darray)) == "xarray.core.dataset.Dataset"):
        print("this function doesn't accept datasets, convert to data array")
        sys.exit()

    ###if it's an xarray dataset, set up the dimensions
    try:
        dims = darray.dims
 
        nmemin = darray[dims[0]].size

        dimboot = [nyears*nmems*nboots]
        dimboot2d = [nboots, nmems, nyears]
#        bootcoords = [('iboot', np.arange(0,nboots,1)), ('imem', np.arange(0,nmems,1)), ('isample', np.arange(0,nyears,1))]
        bootcoords={'iboot': np.arange(0,nboots,1), 'imem':np.arange(0,nmems,1), 'isample':np.arange(0,nyears,1)}
        for icoord in range(1,len(dims)):
            dimboot.append(darray[dims[icoord]].size)
            dimboot2d.append(darray[dims[icoord]].size)
            #bootcoords.append( (dims[icoord], darray[dims[icoord]] ))
            bootcoords[dims[icoord]] = darray[dims[icoord]]

#        print("you are using an xarray dataarray")
   
    except:
        nmemin = darray.shape[0]

        dimboot = [nyears*nboots]
        dimboot2d = [nboots, nyears]
        for icoord in range(1,len(darray.shape)):
            dimboot.append(darray.shape[icoord])
            dimboot2d.append(darray.shape[icoord])

        print("you are using a numpy array")

    ### generate random number for bootstrapping
    if (seed):
        np.random.seed(seed)

    ### do the resampling
    ranu = np.random.uniform(0,nmemin-nyears,nboots*nmems)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros(dimboot2d)
    for iboot in np.arange(0,nboots,1):
        for imem in np.arange(0,nmems,1):
            bootdat[iboot,imem,...] = darray[ranu[iboot*nmems + imem]:ranu[iboot*nmems+imem]+nyears]
   

#    bootdat = np.array(darray[ranu])
#    bootdat = bootdat.reshape(dimboot2d)

#    print(bootdat)
#    print(bootcoords)


    try:
        bootdat = xr.DataArray(bootdat, coords=bootcoords)
    except:
        pass

    return bootdat


def boot_regcoefs(a1,a2,sigx=None,sigy=None,nboots=1000):
    """ Output bootstrap samples of regression coefficients

    Input:
        a1 = first array
        a2 = second array
    Optional input:
        nboots = the number of bootstrap samples used to generate the ci
        sigx = the standard deviation on the predictor points
        sigy = the standard deviation on the predictand points

    Output:
        acoefs = nboots samples of the coefficient a 
        bcoefs = nboots samples of the coefficient b
    
    where y = a + bx
 
    Different regression methods are used 
    depending on the availability of sigx or sigy
    if no sigx then ordinary least squares regression
    if sigx and sigy then total least squares regression
    """

    if (a1.size != a2.size):
        print("The two arrays must have the same size")
        sys.exit()

    samplesize = a1.size
    ranu = np.random.uniform(0,samplesize,nboots*samplesize)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros([samplesize,nboots])
    bootdat1 = np.array(a1[ranu])
    bootdat2 = np.array(a2[ranu])
    bootdat1 = bootdat1.reshape([samplesize,nboots])
    bootdat2 = bootdat2.reshape([samplesize,nboots])

    if sigx is not None:
        bootdatsigx = np.array(sigx[ranu])
        bootdatsigx = bootdatsigx.reshape([samplesize,nboots])
    if sigy is not None:
        bootdatsigy = np.array(sigy[ranu])
        bootdatsigy = bootdatsigy.reshape([samplesize,nboots])

    acoef = np.zeros(nboots) ; bcoef=np.zeros(nboots)

    if sigx is not None:
        for iboot in range(0,nboots,1):
            acoef[iboot], bcoef[iboot] = linfit.tls(bootdat1[:,iboot],bootdat2[:,iboot],
                              bootdatsigx[:,iboot],bootdatsigy[:,iboot])

    else:
        for iboot in range(0,nboots,1):
            acoef[iboot], bcoef[iboot] = linfit.linfit_xy(bootdat1[:,iboot],bootdat2[:,iboot])

    return acoef, bcoef

def boot_corcoefs(a1, a2, nboots=1000):
    """ Output bootstrap samples of correlation coefficients """
    if (a1.size != a2.size):
        print("The two arrays must have the same size")
        sys.exit()

    samplesize = a1.size
    ranu = np.random.uniform(0,samplesize,nboots*samplesize)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros([samplesize,nboots])
    bootdat1 = np.array(a1[ranu])
    bootdat2 = np.array(a2[ranu])
    bootdat1 = bootdat1.reshape([samplesize,nboots])
    bootdat2 = bootdat2.reshape([samplesize,nboots])

    bootdat1 = xr.DataArray(bootdat1, dims=['model','boot'])
    bootdat2 = xr.DataArray(bootdat2, dims=['model','boot'])

    rvals = xr.corr(bootdat1, bootdat2, dim='model')

    return rvals

def boot_corsignif_multidim(a1, a2, dim, nboots=1000, seed=None, signan=True):
    """ Output bootstrap significance results for correlation of two arrays over dimension dim """

    nsamples = a1[dim].size

    if (seed):
        np.random.seed(seed)

    ranu = np.random.uniform(0,nsamples, nboots*nsamples)
    ranu = np.floor(ranu).astype(int)

    #Figure out the coordinates
    dims_1 = a1.dims
    dimboot_1 = [nboots]
    bootcoords_1 = [ ("iboot", np.arange(0,nboots,1))]
    for icoord in range(0,len(dims_1)):
        dimboot_1.append(a1[dims_1[icoord]].size)
        bootcoords_1.append( (dims_1[icoord], a1[dims_1[icoord]].values ))

    dims_2 = a2.dims
    dimboot_2 = [nboots]
    bootcoords_2 = [ ("iboot", np.arange(0,nboots,1))]
    for icoord in range(0,len(dims_2)):
        dimboot_2.append(a2[dims_2[icoord]].size)
        bootcoords_2.append( (dims_2[icoord], a2[dims_2[icoord]].values ))

    kwargs = {dim: ranu}
    a1boot = a1.isel(**kwargs)
    a2boot = a2.isel(**kwargs)

    a1boot_reshape = np.reshape(np.array(a1boot), dimboot_1)
    a2boot_reshape = np.reshape(np.array(a2boot), dimboot_2)


    a1boot_xr = xr.DataArray(a1boot_reshape, coords=bootcoords_1)
    a2boot_xr = xr.DataArray(a2boot_reshape, coords=bootcoords_2)

    rvals = xr.corr(a1boot_xr, a2boot_xr, dim=dim)
    min95 = rvals.quantile(0.025, dim='iboot')
    max95 = rvals.quantile(0.975, dim='iboot')

    signif = min95*0 + 1
    if (signan): # significant values are nan
        signif = signif.where( (min95 < 0) & (max95 > 0), nan)
    else:
        signif = signif.where( (min95 > 0) | (max95 < 0), nan)

    return signif 








  




