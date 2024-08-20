burstatmos
========

Xspec model for fitting the atmosphere models of [Suleimanov et al. (2012)](http://adsabs.harvard.edu/abs/2012A%26A...545A.120S) to the burst data. 

### Three versions included:
 - burstatmos.f90
    * The preferred version. The pre-calculated atmosphere data is hardcoded in this version. The data consist of over 13500 lines, which can make the handling of this code a bit difficult. Compiling this can take a long time (tens of minutes), but running this is code is very fast, fitting one burst takes roughly one minute. 
 - burstatmos_read.f90
    * The atmosphere data are stored in files `ells.dat`,` teffs.dat`, and `fluxes.dat` instead of hardcoding them into the model. Running this version is much slower than the one above, but compiling takes only seconds and handling and debugging is easier. The paths to the datafiles need to be set manually in the code at the moment. 
 - burstatmos_calculate.f90
    * This is the "raw" version. This version reads directly the atmosphere datafiles of [Suleimanov et al. (2012)](http://adsabs.harvard.edu/abs/2012A%26A...545A.120S) (in folder J_A+A_545_A120, available at [Vizier](http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/545/A120)) and calculates and interpolates the needed ells, teffs, and fluxes. Running this code is very slow and it should not be used for the fitting, but it can be (and has been) used for calculating the data in files `ells.dat`,` teffs.dat`, and `fluxes.dat` and the data hardcoded in `burstatmos.f90`. 


### Installing:
  - Copy the fortran code and the required file `lmodel.dat` to the local model directory, e.g. /home/UTU/jmjkuu/heasoft-6.21/local/. The default local model directory can be set in file ~/.xspec/Xspec.init (see also [Customizing Xspec](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/node33.html))
  - Then launch Xspec and build the model with [initpackage](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node101.html#initpackage)
    ```sh
    XSPEC12> initpackage burstatmos lmodel.dat /home/UTU/jmjkuu/heasoft-6.21/local/
    ```
    This can take up to an hour if using the model with hardcoded data, but it needs to be done only once.
  - The model can then be loaded with [lmod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node102.html#lmod)
    ```sh
    XSPEC12> lmod burstatmos /home/UTU/jmjkuu/heasoft-6.21/local/
    ```
    The model loading can be done automatically on start-up by adding this command to the file ~/.xspec/xspec.rc (see [Xspec FAQ](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/faq.html#Q8)).
    
  - See also this [Xspec model guide](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node294.html).

  
### Parameters:
This model has five free parameters:
   - 1. g_rad/g
   - 2. Hydrogen mass fraction X
   - 3. Mass (in solar)
   - 4. Radius (in km)
   - 5. Distance (in kpc)

The initial values and hard and soft limits for each parameter are defined in `lmodel.dat` (see [Xspec model guide](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node294.html)).
