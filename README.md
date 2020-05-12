# ExtremeEddies

## Tracking overview

Eddy tracking algorithms, including a distinctive third eddy type [anticyclonic modewater eddy (ACME)] identification procedure in addition to the standard cyclones and anticylones, as well as additional eddy dynamics metrics including rotational velocity, translational velocity, and nonlinearity.

With a focus on the biogeochemical consequences of nonlinear eddies.

### N.B.

**Baseline structure for the tracking is from Eric Oliver's implementation of the Chelton et al. (2011) method - please see his original code at https://github.com/ecjoliver/eddyTracking and the methodology in Chelton et al. Progress in Oceanography, 2011. Additional code is written by Jamie Atkins to create a modified tracking procedure with an additional eddy type and additional eddy dynamics metrics for future analysis.**

**The author makes no claim of Eric Oliver's work to be their own, and instead intends for it to be a useful set of complementary functions for novel dynamics/climatological/biogeochemical eddy analysis.**

## Contents

| File | Description |
| ---- | ----------- |
| paramsADD.py | Parameters (adjustable) to be used in subsequent detection, tracking etc. Including original Eric Oliver eddyTracking code but with additional novel code written by Jamie Atkins. |
| eddy_functionsADD.py | Functions to be used throughout detection, tracking process. Including original Eric Oliver eddyTracking functions but with additional novel functions added by Jamie Atkins. |
| eddy_detectionADD.py | Detection procedure for eddy features. Including original Eric Oliver eddyTracking code but with additional novel ('part 2') code written by Jamie Atkins. |
| eddy_trackingADD.py | Tracking procedure for eddies, using detected features form eddy_detectionADD.py. Including original Eric Oliver eddyTracking code but with additional novel ('part 2') code written by Jamie Atkins. |
| anomFunctions.py | Supporting functions for calculating SST/SSS anomalies in detection/tracking procedure. |

## How to use

N.B.
1) SSH data must be netcdf file and the tracking requires individual files for individual days, e.g. 365 individual .nc files for one year of data. If data is currently one file then it can be split easily using the Climate Data Operators (CDO) command line suite 'splitsel' function and then each file can be renamed. In the source code posted in this repository the file naming follows the format of day0.nc, day1.nc, day2, ... , day364.nc etc.

2) For use on different datasets, parameters need to be adjusted. E.g. the `NAME` variable in params.py as well as any other variables (timesteps, resolution etc.). Also functions in eddy_functionsADD.py will also need to be adjusted in order to properly loaded in, see funcitons `load_eta()`, `load_additional()`

3) For additional variables/metrics calculations etc. in eddy_detection.py part 2, if all variables (SSH, SST, SSS, MLD etc.) come in separate files rather than as one netcdf as in the case of this source code then the loading of the variables in part 2 must be adjusted accordingly.

4) Rossby radius data (rossrad.dat) file required is available from Eric Oliver eddyTracking repository; or by `wget http://www-po.coas.oregonstate.edu/research/po/research/rossby_radius/rossrad.dat`
