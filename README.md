# CASanalysis

## Climate analysis python tools

Here we are gathering various python functions and jupyter notebook examples that can be used for common climate analysis tasks, such as reading in large datasets, calculating regional averages, modes, statistical analysis etc.

### Installing Utilities

Assuming you have cloned this github repostory to $DIR, in order to install the python functions located in ./CASutils which form the basis of this package you will need to do the following

```bash
cd $DIR
pip install -e . --user
```


## Code examples

Located in ./examples

* **output_lens1_seasons_prect.ipynb** = reading in 1979-2014 for LENS1 from mojave and calculating yearly seasonal averages.
* **output_lens2_seasons_prect.ipynb** = reading in 1979-2014 for LENS2 from mojave and calculating yeraly seasonal averages.
* **readisddata.ipynb** = reading ISD station data from mojave.  parses inventory files and csv data files and chooses stations based on data availability for a required period.  Outputs to netcdf.
* **grabera5t2m_from_rda.ipynb** = getting daily mean 2m temperature from the ERA5 archive on the CISL RDA.  Reads in hourly data on a 0.25 deg grid, calculates daily mean and regrids to the CESM grid.
* **grabera5zmfluxes.ipynb** = computing daily zonal mean fluxes using full column U, V, T and OMEGA for ERA5 from the CISL RDA.
* **grabera5strd.ipynb** = computing daily averaged radiative fluxes from ERA5 on the CAM grid i.e., this one deals with ERA5 forecast fields.
* **fluxnetread_allstations.ipynb** = reading in fluxnet data. Retaining stations over a specified longitude latitude range and with a specified record length. Removing Feb 29th and outputting to netcdf.
* **outputera5deseas.ipynb** = reading in ERA5 daily data, removing leap years, deseasonalizing by removing harmonics from the seasonal cycle and then the seasonal mean
* **plotvpd_cmip6.ipynb** = an example using lots of plotting techniques: jointpdfs, different colored points etc.  Shows saturation vapor pressure versus actual vapor pressure over the US southwest for CMIP6 models and ERA5.
* **GEOHEAT_jlats.ipynb** = an example reading in some zonal wind data, calculating the DJF zonal mean and the location and speed of the jet maxima in each hemisphere.
