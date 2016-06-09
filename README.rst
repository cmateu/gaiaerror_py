Gaia Error Py-Wrapper
======

This is a simple Python wrapper for Merce Romero-Gomez's `Gaia-Error code <https://github.com/mromerog/Gaia-errors>`__

**DESCRIPTION**

This version of Merce Romero-Gomez's code includes the evaluation of the Drimmel extinction map.


Codes and files provided

- **gaiaerr.py**: Simulates Gaia errors 
- **gaiaerr_photerr.py**: Simulates Gaia errors for all quantities except distance. For the distance, (gaussian) errors are simulated with a fixed relative error. 
- **example_data/myfile.ne.dat**: Example input file


**INSTALL**

- First run the compile script inside the gaia_errors_color_tmission folder
- Edit the gaiaerr.py and gaiaerr_photerr.py files and substitute the appropriate path in the gerr_path variable at the very beginning of the file

**REQUIREMENTS**

The numpy and scipy libraries are required

**QUICK GUIDE**

The codes take input files (ascii) containing  X Y Z VX VY VZ Mv V-I, spatial coordinates must be in kpc and velocities in km/s. Mv and V-I are the absolute magnitude and V-I intrinsic colour.

You can use the example input file provided to run the Gaia error simulation as::

	gaiaerr.py myfile.ne.dat

The output file will be myfile.ge.dat.

To simulate 25% distance errors plus Gaia errors for the remaining quantities, run the gaiaerr_photerr.py code as::

	gaiaerr_photerr.py myfile.ne.dat 0.25

The output file will be myfile.pe.dat.

By default both codes simulate nominal end-of-mission Gaia errors, i.e. the errors expected at the end of the nominal mission lifetime of 5 yr. To simulate errors at a different mission time, use the -tm option followed by the desired mission time in years.

The following example will provide errors simulated for 1.5 yrs of operation::

	gaiaerr.py myfile.ne.dat -tm 1.5

You can also simulate errors for extended mission lifetimes, e.g. if the Gaia mission is extended for 2 more years, simulate the errors as::

	gaiaerr.py myfile.ne.dat -tm 7.

It is assumed the errors scale as t\*\*3/2 for up to 10 years of operation. For larger times a t\*\*1/2 scaling factor is assumed (A.G.A. Brown, private communication) 

**OUTPUT**

For both codes the output file names are the same as the input files with the **.ge.dat** extension for  **gaiaerr.py** and **.pe.dat** for  **gaiaerr_photerr.py**.

Both codes produce output files with the following structure:

- Av, V, Gmag, Grvs : V-band extinction, V-band, G and Grvs apparent magnitudes
- xX,xY,xZ,xVX,xVY,xVZ : error-free cartesian galactocentric coordinates (kpc) and velocities (km/s)
- xl_deg,xb_deg,xRhel_kpc,xmulcosb_masyr,xmub_masyr,xvrad : error-free observables (l,b, heliocentric distance in kpc, proper motions in mas/yr, radial velocities in km/s)
- gX,gY,gZ,gVX,gVY,gVZ : error-convolved cartesian galactocentric coordinates and velocities.
- gl_deg,gb_deg,gRhel_kpc,gmulcosb_masyr,gmub_masyr,xvrad : error-convolved observables
- relerr_D: relative error in distance
- sig_mub: mean error in proper motion (gaussian standard deviation)
- sig_vrad: mean error in radial velocity (gaussian standard deviation)
- VI: observed V-I colour
- relerr_mub: relative error in proper motion
- relerr_vrad: relative error in radial velocity

**AUTHORS**

The Gaia-Errors fortran library and package has been developed by Merce Romero-Gomez (University of Barcelona, UB) and is maintained at https://github.com/mromerog/Gaia-errors.
The python wrapped has been developed by Cecilia Mateu (CIDA).




